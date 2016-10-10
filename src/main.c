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
#include <ctype.h>
#include "sci_const.h"
#include "timing.h"
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
#include "glob_time.h"
#define BUF_SIZE 256

double serial_t;

struct ccfg *cache_cfg;

char **
parse_inpfile (char * fname)
{
  int j, k, l;
  int foundspace = 0;

  FILE *inpfile;

  char **argv_file;
  argv_file = malloc(BUF_SIZE * sizeof(char *));

  char *argv_buf = NULL;
  size_t len = 0;
  ssize_t line_sz;
  int lsz;

  for (j = 0; j < BUF_SIZE; j++) {
    argv_file[j] = malloc(BUF_SIZE * sizeof(char));
  }

  if((inpfile = fopen(fname, "r")) == NULL) {
    fprintf(stderr, "\n\nmain.c, function parse_inpfile: unable to open the input file %s.\n"
            , fname);
    printf("program terminating due to the previous error.\n\n");
    exit(1);
  }

  j = 1;
  k = 0;

  while ((line_sz = getline(&argv_buf, &len, inpfile)) != -1) {
    lsz = (int)line_sz - 1;

    /* printf("Retrieved line of length %d :\n", lsz); */
    /* printf("containing: %s\n", argv_buf); */

    if ((lsz > 0) && (argv_buf[0] != '#')) {

      for (k = 0; k < lsz; k++) {
        if (argv_buf[k] == ' ') {
          foundspace = 1;
        }
      }

      argv_file[j][0] = '-';

      if (foundspace) {

        /* the line contains a space, meaning it is of type "arg params" */
        for (l = 1, k = 0; k < lsz; k++) {

          if (argv_buf[k] == ' ') {
            /* next line */
            argv_file[j][l] = '\0';
            j++;
            l = 0;
          }
          else {
            argv_file[j][l++] = argv_buf[k];
          }
        }
        argv_file[j][l] = '\0';
        foundspace = 0;
      }
      else {
        for (k = 0; k < lsz; k++) {
          argv_file[j][k+1] = argv_buf[k];
        }
        argv_file[j][k+1] = '\0';
      }
      j++;
    }
  }

  sprintf(argv_file[0], "%d", j);

  if (fclose(inpfile) != 0) {
    fprintf(stderr, "\n\nmain.c, function parse_inpfile: unable to close the file:\n%s\n",
            fname);
    printf("program terminating due to the previous error.\n\n");
    exit(1);
  }
  free(argv_buf);

  return argv_file;
}

int
main (int argc, char *argv[])
{
  serial_t = omp_get_wtime(); /* start mark of the serial runtime
                                        of the program */
  int j, k, l, m, n;                    /* iteration variables */
  int n_t;                      /* number of numbers read from input flag */
  int n_er = 0;                     /* number of numbers eigenstate energy values provided */
  m = 0;
  int n_res;
  int len_op = 0;
  int len_fn = 0;
  int out_set = 0;
  int so_enrg = 0;
  int intf_mode = 0;
  int elastic = 0; /* force non-elastic transition mode to be actiated */
  int verbosity = 0;
  int lorz = 0;

  double *res = malloc(2 * sizeof(double));

  char argv_buf[BUF_SIZE];
  /* arrays for storing input file name data */
  char *input_sbuff = malloc(BUF_SIZE);
  char *inp_fn = NULL;
  char *outpath = NULL;
  char *inpath = NULL;
  char *inp_sfx = NULL;
  char *cache_fpstr;
  char *num_buf;
  char *end;

  char *ltime = malloc(20);
  char curr_dir[BUF_SIZE] = {0};
  char out_buf[BUF_SIZE] = {0};

  /* use first element to specify number of values stored in matrix */
  /* threshold values for the three or six different state types */
  double *state_t = malloc(4 * sizeof(double));

  /* initial/final and intermediate state energy ranges */
  double *state_er = malloc(7 * sizeof(double));

  double *gx_inp = NULL;
  double *gy_inp = NULL;
  double *lx_inp = NULL;
  double *ly_inp = NULL;
  double *fwhm_inp = NULL;
  struct metadata *md = init_md();
  struct spectrum *tmp_spec;
  struct inp_node *inp;

  /* set default values for the parameter arrays */
  n_t = 2;
  state_t[0] = n_t;
  state_t[1] = 0.99; /* by default, keep 99% of the total
                       intensity */
  state_t[2] = 0.01;
  res[0] = 0.05;
  res[1] = 0.05;

  state_er[0] = 0;

  struct stat st = {0};

  printf("\n\nscttr calculation initiated (%s).\n\n", get_loctime(ltime));

  /* process the input arguments */
  if (argc == 1) {
    printf("No command line arguments provided. To read the help documentation, provide the argument \"-h\". Program terminating.\n");
    fflush(stdout);
  }

  /* if only three arguments were provided, assume that the input parameters
     can be found in the path of the input file. */
  if (argc == 3) {
    printf("program provided with the following input file: %s\n\n",argv[2]);
    argv = parse_inpfile(argv[2]);
    argc = strtol(argv[0], &end, 10);

    /* write the input file as a header to the log file */
    for (j = 1; j < argc; j++) {
      printf("%s\n", argv[j]);
    }
    printf("\n\n");
  }

  while (argc > 1 && (argv[1][0] == '-')) {

    n = 0;
    while(argv[1][n] != '\0') {
      argv_buf[n] = argv[1][n];
      n++;
    }
    argv_buf[n++] = '\0';

    if(strstr(argv_buf,"v") || strstr(argv_buf,"verbosity")) {
      num_buf = malloc(2);
      /* num_buf[1] = '\0'; */
      num_buf[0] = argv[2][0];
      verbosity = (int)atof(num_buf);

      free(num_buf);
      num_buf = NULL;

      argv++;
      argc--;
    }
    else if(strstr(argv_buf,"e")) {
      n_er = 1;

      for (j = 0, k = 0; argv[2][j] != '\0'; j++) {

        input_sbuff[k++] = argv[2][j];

        if (((argv[2][j] == ',')
             || (argv[2][j] == ';'))
            || (argv[2][j + 1] == '\0')) {
          num_buf = malloc(k + 1);

          for (l = 0; l < k; l++) {
            num_buf[l] = input_sbuff[l];
          }
          num_buf[k] = '\0';

          state_er[n_er] = atof(num_buf);
          free(num_buf);
          num_buf = NULL;
          n_er++;

          k = 0;
        }
      }
      state_er[0] = n_er - 1;
      argv++;
      argc--;
    }
    else if(strstr(argv_buf, "F")) {
      so_enrg = 1;
    }
    else if(strstr(argv_buf,"gx")) {

      m = 0;
      for (j = 0; argv[2][j] != '\0'; j++) {
        if (((argv[2][j] == ',')
             || (argv[2][j] == ';'))
            || (argv[2][j + 1] == '\0')) {
          m++;
        }
      }
      if (m > 1) {
        gx_inp = malloc((m + 2) * sizeof(double));
        gx_inp[0] = m;
        m = 1;
        for (j = 0, k = 0; argv[2][j] != '\0'; j++) {

          input_sbuff[k++] = argv[2][j];
          if (((argv[2][j] == ',')
               || (argv[2][j] == ';'))
              || (argv[2][j + 1] == '\0')) {
            num_buf = malloc(k + 1);
            for (l = 0; l < k; l++) {
              num_buf[l] = input_sbuff[l];
            }
            num_buf[k] = '\0';

            gx_inp[m++] = atof(num_buf);
            free(num_buf);
            num_buf = NULL;

            k = 0;
          }
        }
      }
      argv++;
      argc--;
    }
    else if(strstr(argv_buf,"gy")) {

      /* the user specified thresholds for each state type (ground, initial,
         final) to be used for screening the states when calculating the map */
      m = 0;
      for (j = 0; argv[2][j] != '\0'; j++) {
        if (((argv[2][j] == ',')
             || (argv[2][j] == ';'))
            || (argv[2][j + 1] == '\0')) {
          m++;
        }
      }
      if (m > 1) {
        gy_inp = malloc((m + 2) * sizeof(double));
        gy_inp[0] = m;
        m = 1;
        for (j = 0, k = 0; argv[2][j] != '\0'; j++) {

          input_sbuff[k++] = argv[2][j];
          if (((argv[2][j] == ',')
               || (argv[2][j] == ';'))
              || (argv[2][j + 1] == '\0')) {
            num_buf = malloc(k + 1);
            for (l = 0; l < k; l++) {
              num_buf[l] = input_sbuff[l];
            }
            num_buf[k] = '\0';

            gy_inp[m++] = atof(num_buf);
            free(num_buf);
            num_buf = NULL;

            k = 0;
          }
        }
      }
      argv++;
      argc--;
    }
    else if(strstr(argv_buf, "h")) {
      printf("See the \"Usage\" section of the documentation provided in the doc directory of the program.\n");
      exit(1);
    }
    else if(strstr(argv_buf, "i")) {

      k = 0;

      /* extract the file name */
      for (j = 0; argv[2][j] != '\0'; j++) {
        if (argv[2][j] == '/') {
          k = j + 1;
        }
      }

      inp_fn = malloc(j - k + 1);

      for (l = 0, j = k; argv[2][j] != '\0'; j++) {

        if (argv[2][j] == '.') {
          inp_fn[l] = '\0';
          j--;
          break;
        }
        inp_fn[l] = argv[2][j];
        l++;
      }

      for (j = 0; argv[2][j] != '\0'; j++) {
        input_sbuff[j] = argv[2][j];
      }

      /* set the path to the input file */
      len_fn = j;
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

      printf("  - input path specified  = %s\n", inpath);

      argv++;
      argc--;
    }
    else if(strstr(argv_buf,"lx")) {

      m = 0;
      for (j = 0; argv[2][j] != '\0'; j++) {
        if (((argv[2][j] == ',')
             || (argv[2][j] == ';'))
            || (argv[2][j + 1] == '\0')) {
          m++;
        }
      }
      if (m > 1) {
        lx_inp = malloc((m + 2) * sizeof(double));
        lx_inp[0] = m;
        m = 1;
        for (j = 0, k = 0; argv[2][j] != '\0'; j++) {

          input_sbuff[k++] = argv[2][j];
          if (((argv[2][j] == ',')
               || (argv[2][j] == ';'))
              || (argv[2][j + 1] == '\0')) {
            num_buf = malloc(k + 1);
            for (l = 0; l < k; l++) {
              num_buf[l] = input_sbuff[l];
            }
            num_buf[k] = '\0';

            lx_inp[m++] = atof(num_buf);
            free(num_buf);
            num_buf = NULL;

            k = 0;
          }
        }
      }
      lorz = 1;
      argv++;
      argc--;
    }
    else if(strstr(argv_buf,"ly")) {

      m = 0;
      for (j = 0; argv[2][j] != '\0'; j++) {
        if (((argv[2][j] == ',')
             || (argv[2][j] == ';'))
            || (argv[2][j + 1] == '\0')) {
          m++;
        }
      }
      if (m > 1) {
        ly_inp = malloc((m + 2) * sizeof(double));
        ly_inp[0] = m;
        m = 1;
        for (j = 0, k = 0; argv[2][j] != '\0'; j++) {

          input_sbuff[k++] = argv[2][j];
          if (((argv[2][j] == ',')
               || (argv[2][j] == ';'))
              || (argv[2][j + 1] == '\0')) {
            num_buf = malloc(k + 1);
            for (l = 0; l < k; l++) {
              num_buf[l] = input_sbuff[l];
            }
            num_buf[k] = '\0';

            ly_inp[m++] = atof(num_buf);
            free(num_buf);
            num_buf        = NULL;

            k = 0;
          }
        }
      }
      lorz = 1;
      argv++;
      argc--;
    }
    else if(strstr(argv_buf,"o")) {
      /* printf("processing input file: "); */
      out_set = 1;

      /* extract the output path */
      for (j = 0; argv[2][j] != '\0'; j++) {
        input_sbuff[j] = argv[2][j];
      }

      len_op = j;
      outpath = malloc(len_op + 1);

      for (j = 0; j < len_op; j++) {
        outpath[j] = input_sbuff[j];
      }

      outpath[len_op] = '\0';
      printf("  - output path specified  = %s\n", outpath);

      argv++;
      argc--;
    }
    else if(strstr(argv_buf,"r")) {
      n_res = 0;

      for (j = 0, k = 0; argv[2][j] != '\0'; j++) {

        input_sbuff[k++] = argv[2][j];
        if (((argv[2][j] == ',')
             || (argv[2][j] == ';'))
            || (argv[2][j + 1] == '\0')) {

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

      argv++;
      argc--;
    }
    else if(strstr(argv_buf,"t")) {
      n_t = 1;

      for (j = 0, k = 0; argv[2][j] != '\0'; j++) {

        input_sbuff[k++] = argv[2][j];

        if (((argv[2][j] == ',')
             || (argv[2][j] == ';'))
            || (argv[2][j + 1] == '\0')) {
          num_buf = malloc(k + 1);

          for (l = 0; l < k; l++) {
            num_buf[l] = input_sbuff[l];
          }
          num_buf[k] = '\0';

          state_t[n_t] = atof(num_buf);

          free(num_buf);
          num_buf = NULL;

          if ((n_t) == 3) {
            break;
          }
          n_t++;
          k = 0;
        }
      }
      state_t[0] = n_t;
      argv++;
      argc--;
    }
    else if(strstr(argv_buf,"I")) {
      intf_mode = 1;
    }
    else if(strstr(argv_buf,"E")) {
      elastic = 1;
    }
    else{
      fprintf(stderr, "main.c, function main: Unknown flag %s. See the \"Usage\" section of the documentation provided in the doc directory of the program.\n"
              , argv_buf);
      printf( "program terminating due to the previous error.\n");
      exit(EXIT_FAILURE);
    }
    argv++;
    argc--;
  }

  if (fwhm_inp == NULL) { /* set default broadenings */
    fwhm_inp = malloc(4 * sizeof(double));
    fwhm_inp[0] = 1;
    fwhm_inp[1] = 0;
    fwhm_inp[2] = 0.5;
    fwhm_inp[3] = 0.5;
    fwhm_inp[3] = 0.5;
  }

  set_ccnuma_affinity();
  cache_cfg = set_ccfg(0.5); /* use 90% of cache space */

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

  /* process the input */
  stat(inpath, &st);

  gx_inp[m] = (state_er[3] + 9) / AUTOEV;
  gy_inp[m] = (state_er[6] + 9) / AUTOEV;

  if (lorz) {
    lx_inp[m] = gx_inp[m];
    ly_inp[m] = gy_inp[m];
  }

  if (state_t[1] > 1) {
    state_t[1] = state_t[1] / 100;
  }

  md -> sz_inp = (int)st.st_size;
  md -> so_enrg = so_enrg;

  md -> outpath = outpath;
  md -> inpath = inpath;
  md -> inp_fn = inp_fn;
  md -> inp_sfx = inp_sfx;

  md -> state_er  = state_er;
  md -> state_t = state_t;
  md -> res = res;
  md -> gx = gx_inp;
  md -> gy = gy_inp;
  md -> lx = lx_inp;
  md -> ly = ly_inp;
  md -> lorz = lorz;
  md -> intf_mode = intf_mode; /* only constructive interference implemented for now */
  md -> v = verbosity;

  l = j;
  inp = init_inp(md);

  if ((elastic == 1)
      && ((md -> state_er[1] == md -> state_er[5])
      && (md -> state_er[2] == md -> state_er[6]))){
    inp -> el = 1;
  }
  else {
    inp -> el = 0;
  }

  if (len_fn == 0) {
   fprintf(stderr, "\n\Error: calc_spec.c, main: you didnt provide the path to an input file.");
    printf("program terminating due to the previous error.\n");
    exit(EXIT_FAILURE);
  }
  else {
    printf("  - parsing input data (%s)..", get_loctime(ltime));
    fflush(stdout);

    if (parse_input(inp)) {
      fprintf(stderr, "calc_spec.c, main: unable to parse the input data contained in %s.\n"
              , inpath);
      printf("program terminating due to the previous error.\n");
      exit(EXIT_FAILURE);
    } else {
      fflush(stdout);
      printf("\n    done (%s).\n", get_loctime(ltime));
    }
  }

  printf("\n executing scttr with the following..\n\n");
  printf("  - data contained in the input file:\n    %s\n\n"
         , inpath);
  printf("  - intensity and ground state boltzmann weight threshold values respectively:\n    ");

  for (j = 1; j < n_t; j++) {
    printf("%le, ", state_t[j]);
  }

  printf("\n\n");

  printf("  - state energy intervals values:\n    ");

  for (j = 1; j < n_er; j++) {
    printf("%le, ", state_er[j]);
  }

  printf("\n\n");
  printf(" execution progress:\n\n");
  printf("  - screening transitions, resulting in the removal of (%s)..\n", get_loctime(ltime));

  set_spec(inp);
  tmp_spec = get_spec(inp, 2);
  printf("      .. %d transitions (%f%% of the total intensity) on ground state boltzmann weight.\n", tmp_spec -> n_sst_bw, (tmp_spec -> iscr_bw / tmp_spec -> itot)*100);
  printf("      .. %d transitions (%f%% of the total intensity) on the total transition intensity.\n", tmp_spec -> n_sst_bw, (tmp_spec -> iscr_int / (tmp_spec -> itot - tmp_spec -> iscr_bw))*100);
  printf("    done (%s).\n", get_loctime(ltime));

  printf("  - forming the reduced transition matrix (%s)..", get_loctime(ltime));
  fflush(stdout);

  set_trs_red(inp, 2);

  /* fflush(stdout); */
  /* strs2str(inp, get_spec(inp,1)); */
  /* fflush(stdout); */
  /* printf("\n\n" ); */
  /* strs2str(inp, get_spec(inp,2)); */
  /* fflush(stdout); */
  /* fprintf(stderr, "\n\n=======Valgrind eject point=======\n\n"); */
  /* fflush(stdout); */
  /* exit(1); */
  printf(" done (%s).\n", get_loctime(ltime));

  calc_spec(inp, 2);

  /* trs2str(get_spec(inp,2)); */

  write_spec(inp, get_spec(inp,2));
  write_plotscript(inp, get_spec(inp,2));

  if ((verbosity == 1) || (verbosity == 3)) {
    write_sticks(inp, get_spec(inp,2), md);
    cache_fpstr = concs(2, inp -> md -> outpath,
                        "cache_info.txt");
    write_timings(inp);
    cache2file(cache_fpstr);
    free(cache_fpstr);
  }

  if ((verbosity == 2) || (verbosity == 3)) {
    strs2str(inp, get_spec(inp,1));
    printf("\n\n" );
    strs2str(inp, get_spec(inp,2));
  }

  free_inp(inp);
  free(input_sbuff);

  printf("\nscttr calculation successfully executed.\n");
  printf(" program terminating (%s).\n\n", get_loctime(ltime));
  free(ltime);
  return 0;
}
