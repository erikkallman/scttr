/* Copyright (C) 2015 Erik Källman */
/* This file is part of the scttr program. */

/* scttr is free software: you can redistribute it and/or modify */
/* it under the terms of the GNU Lesser General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or */
/* (at your option) any later version. */

/* scttr is distributed in the hope that it will be useful, */
/* but without any warranty; without even the implied warranty of */
/* merchantability or fitness for a particular purpose. See the */
/* GNU Lesser General Public License for more details. */

/* You should have received a copy of the GNU Lesser General Public License */
/* along with scttr, found in the "license" subdirectory of the root */
/* directory of the scttr program. */
/* If not, see <http://www.gnu.org/licenses/>. */
/**
   * @file main.c
   * @author Erik Källman
   * @date November 2015
   * @brief The main() function handles all user i/o, as well as higher
   * level function calls to parse the provided input, use that data to
   * calculate a spectrum, and write that spectrum to a data file, among
   * other things. it handles memory deallocation once the program enters
   * termination.
   * @note side effects: intializes an input node (@p inp_node) and metadata (@p metadata)
   * struct for the provided input file.
   */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>
#include <omp.h>
#include "std_char_ops.h"
#include "get_nums.h"
#include "calc_spec.h"
#include "formats.h"
#include "scttr_io.h"
#include "scttr_cfg.h"
#include "spectrum.h"
#include "spectrum_s.h"
#include "inp_node_s.h"
#include "cpu_opt.h"
#include "cache_opt.h"

#define BUF_SIZE 256

struct ccfg *cache_cfg;

int
main (int argc, char *argv[])
{

  printf("\n\n scttr calculation initiated.\n\n");
  int j, k, l, m;                    /* iteration variables */
  int n_t;                      /* number of numbers read from input flag */
  int n_er = 0;                     /* number of numbers eigenstate energy values provided */
  int n_res;
  int len_op = 0;
  int len_fn = 0;
  int out_set = 0;
  int so_enrg = 0;

  double *res = malloc(2 * sizeof(double));
  /* arrays for storing input file name data */
  char *input_sbuff = malloc(BUF_SIZE);
  char *inp_fn = NULL;
  char *outpath = NULL;
  char *inpath = NULL;
  char *inp_sfx = NULL;
  char *num_buf;

  char curr_dir[BUF_SIZE] = {0};
  char out_buf[BUF_SIZE] = {0};

  /* use first element to specify number of values stored in matrix */
  /* threshold values for the three or six different state types */
  double *state_t = malloc(4 * sizeof(double));

  /* initial/final and intermediate state energy ranges */
  double *state_er = malloc(7 * sizeof(double));
  double *fwhm_inp = malloc(2 * sizeof(double));

  struct metadata *md = init_md();
  struct inp_node *inp;

  /* set default values for the parameter arrays */
  fwhm_inp[0] = 0.5;
  fwhm_inp[1] = 0.5;

  n_t = 2;
  state_t[0] = n_t;
  state_t[1] = 0.01; /* by default, keep 99% of the total
                       intensity */
  state_t[2] = state_t[1] * 0.0001;
  res[0] = 0.05;
  res[1] = 0.05;

  state_er[0] = 0;

  struct stat st = {0};

  /* process the input arguments */
  if (argc == 1) {
    printf("No command line arguments provided. To read the help documentation, provide the argument \"-h\". Program terminating.\n");
    fflush(stdout);
  }

  while (argc > 1 && (argv[1][0] == '-')) {

    switch(argv[1][1]) {

    case 'e' :

      n_er = 1;

      for (j = 3, k = 0; argv[1][j] != '\0'; j++) {

        input_sbuff[k++] = argv[1][j];

        if ((argv[1][j] == ',') || (argv[1][j + 1] == '\0')) {
          num_buf = malloc(k + 1);

          for (l = 0; l < k; l++) {
            num_buf[l] = input_sbuff[l];
          }
          num_buf[k] = '\0';

          state_er[n_er] = atof(num_buf);
          free(num_buf);
          num_buf        = NULL;
          n_er++;

          k = 0;
        }
      }
      state_er[0] = n_er - 1;

      break;

    case 'f' :

      /* the user specified thresholds for each state type (ground, initial,
         final) to be used for screening the states when calculating the map */
      m = 0;
      for (j = 3, k = 0; argv[1][j] != '\0'; j++) {

        input_sbuff[k++] = argv[1][j];
        if ((argv[1][j] == ',') || (argv[1][j + 1] == '\0')) {
          num_buf = malloc(k + 1);

          for (l = 0; l < k; l++) {
            num_buf[l] = input_sbuff[l];
          }
          num_buf[k] = '\0';

          fwhm_inp[m++] = atof(num_buf);
          free(num_buf);
          num_buf        = NULL;

          k = 0;
        }
      }
      break;

    case 'F' :

      so_enrg = 1;
      break;

    case 'h' :
      printf("See the \"Usage\" section of the documentation provided in the doc directory of the program.\n");
      exit(1);
      break;

    case 'i' :

      k = 0;
      /* extract the file name */
      for (j=3;argv[1][j] != '\0'; j++) {
        if (argv[1][j] == '/') {
          k = j + 1;
        }
      }

      inp_fn = malloc(j - k + 1);

      for (l = 0, j = k; argv[1][j] != '\0'; j++) {

        if (argv[1][j] == '.') {
          j--;
          break;
        }

        inp_fn[l] = argv[1][j];
        l++;
      }

      for (j = 3; argv[1][j] != '\0'; j++) {
        input_sbuff[j - 3] = argv[1][j];
      }

      len_fn = j - 3;
      inpath = malloc(len_fn + 1);

      for (j = 0; j < len_fn; j++) {
        inpath[j] = input_sbuff[j];
      }
      inpath[len_fn] = '\0';

      /* set the inp_sfx  */
      j = len_fn;
      while(inpath[--j] != '.') {};

      inp_sfx = malloc(len_fn - j + 1);

      for (k = 0; j <= len_fn; j++) {
        inp_sfx[k] = inpath[j];
        k++;
      }

      break;

    case 'o' :
      /* printf("processing input file: "); */
      out_set = 1;

      /* extract the output path */
      for (j = 3; argv[1][j] != '\0'; j++) {
        input_sbuff[j - 3] = argv[1][j];
      }

      len_op = j - 3;
      outpath = malloc(len_op + 1);

      for (j = 0; j < len_op; j++) {
        outpath[j] = input_sbuff[j];
      }

      outpath[len_op] = '\0';
      printf("  - outpath specified  = %s\n",outpath );
      break;

    case 'r' :
      n_res = 0;

      for (j = 3,k=0; argv[1][j] != '\0'; j++) {

        input_sbuff[k++] = argv[1][j];
        if ((argv[1][j] == ',') || (argv[1][j + 1] == '\0')) {

          num_buf = malloc(k + 1);

          for (l = 0; l < k; l++) {
            num_buf[l] = input_sbuff[l];
          }
          num_buf[k] = '\0';

          res[n_res] = atof(num_buf);
          free(num_buf);
          num_buf = NULL;
          n_res++;

          k = 0;
        }
      }
      break;

    case 't' :
      n_t = 1;

      for (j = 3, k = 0; argv[1][j] != '\0'; j++) {

        input_sbuff[k++] = argv[1][j];

        if ((argv[1][j] == ',') || (argv[1][j + 1] == '\0')) {
          num_buf        = malloc(k + 1);

          for (l = 0; l < k; l++) {
            num_buf[l] = input_sbuff[l];
          }
          num_buf[k] = '\0';

          state_t[n_t] = atof(num_buf);
          free(num_buf);
          num_buf      = NULL;

          if ((n_t) == 3) {
            break;
          }
          n_t++;
          k = 0;
        }
      }

      state_t[1] = 1 - state_t[1] / 100;
      if (n_t == 1) {
        state_t[2] = state_t[1] * state_t[2];
      }

      state_t[0] = n_t;

      break;

    default :
      fprintf(stderr, "main.c, function main: Unknown flag %c\n", argv[1][1]);
      printf( "program terminating due to the previous error.\n");
      exit(EXIT_FAILURE);

    }
    argv++;
    argv++;
    argc--;
    argc--;
  }

  set_ccnuma_affinity();
  cache_cfg = set_ccfg(0.9); /* use 90% of cache space */

  if (!out_set ) {
    if (getcwd(curr_dir,sizeof(curr_dir)) != NULL) {
      printf("  - output path not specified in input, it defaults to: %s\n",
             curr_dir);
    } else {
      fprintf(stderr, "scttr, main.c: unable to obtain the current directory for use as default output path.  \n");
      printf("program terminating due to the previous error.\n");
      exit(1);
    }
    for (j = 0; curr_dir[j] != '\0'; j++) {};

    outpath = malloc(j);
    for (k = 0; k <= j; k++) {
      outpath[k] = curr_dir[k];
    }
    len_op = j - 1;
  }

  /* check if the output directory exists. if not, create it */
  if (stat(outpath, &st) == -1) {
    printf("      output path not existing. creating directory.\n");

    out_buf[0] = outpath[0];
    /* expect the output path to start with a backslashe */
    for (j = 1; j < len_op; j++) {

      out_buf[j] = outpath[j];
      if (out_buf[j] == '/') {
        out_buf[j + 1] = '\0';

        /* some directories in the outpath might exists, so check
         for that*/
        if (stat(out_buf, &st) == -1) {
          mkdir(out_buf, 0777);
        }
      }
    }
  }

  stat(inpath, &st);
  md -> sz_inp = (int)st.st_size;
  md -> so_enrg = so_enrg;

  md -> outpath = outpath;
  md -> inpath = inpath;
  md -> inp_fn = inp_fn;
  md -> inp_sfx = inp_sfx;

  md -> state_er  = state_er;
  md -> state_t = state_t;
  md -> res = res;
  md -> fwhm = fwhm_inp;
  md -> intf_mode = 0; /* only constructive interference implemented for now */

  l = j;
  inp = init_inp(md);

  if (len_fn == 0) {
   fprintf(stderr, "\n\Error: calc_spec.c, main: you didnt provide the path to an input file.");
    printf("program terminating due to the previous error.\n");
    exit(EXIT_FAILURE);
  }
  else {
    printf("  - parsing input data ..");
    fflush(stdout);

    if (parse_input(inp)) {
      fprintf(stderr, "calc_spec.c, main: unable to parse the input data contained in %s.\n"
              , inpath);
      printf("program terminating due to the previous error.\n");
      exit(EXIT_FAILURE);
    } else {
      printf("    done.\n");
    }
  }

  printf("\n executing scttr with the following..\n\n");
  printf("  - data contained in the input file:\n    %s\n\n"
         , inpath);
  printf("  - intensity and ground state boltzmann weight threshold values respectively:\n    ");

  for (j = 0; j < n_t; j++) {
    printf("%le, ", state_t[j + 1]);
  }

  printf("\n\n");

  printf("  - state energy intervals values:\n    ");

  for (j = 1; j < n_er; j++) {
    printf("%le, ", state_er[j]);
  }

  printf("\n\n");
  printf(" execution progress:\n\n");

  set_spec(inp);

  set_trs_red(inp, 2);
  calc_spec(inp, 2);

  write_spec(inp, get_spec(inp,2));
  write_plotscript(inp, get_spec(inp,2));

  free(inp_sfx);
  free_inp(inp);

  free(input_sbuff);
  free(state_er);
  free(state_t);
  printf("\n scttr successfully executed.\n");
  printf(" program terminating.\n\n");

  return 0;
}
