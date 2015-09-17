#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>
#include "spec_info.h"
#include "structs.h"
#include "state.h"
#include "std_char_ops.h"
#include "get_numsl.h"
#include "smap.h"
#include "formats.h"
#include "parse_input.h" /* for input_data array, as well as other
                            calculation-specific variables */
#include "smap_cfg.h"

#define BUF_SIZE 256
int sz_inp;
int etype;

int
main (int argc, char * argv[]) {

  printf( "\n\n smap calculation initiated.\n\n" );
  int j,k,l,m;                    /* iteration variables */
  int n_t;                      /* number of numbers read from input flag */
  int n_er = 0;                     /* number of numbers eigenstate energy values provided */
  int n_res;
  int len_op   = 0;
  int len_fn   = 0;
  int out_set  = 0;
  etype = 0;

  double * res = malloc(2*sizeof(double));
  /* arrays for storing input file name data */
  char *   input_sbuff = malloc(BUF_SIZE);
  char *   inp_fn = NULL;
  char *   outpath = NULL;
  char *   inpath = NULL;
  char *   inp_sfx = NULL;
  char *   num_buf;

  char curr_dir[BUF_SIZE] = {0};
  char out_buf[BUF_SIZE] = {0};

  /* use first element to specify number of values stored in matrix */
  /* threshold values for the three or six different state types */
  double * state_t = malloc(4*sizeof(double));

  /* initial/final and intermediate state energy ranges */
  double * state_er = malloc(7*sizeof(double));
  double * fwhm_inp = malloc(2*sizeof(double));

  metadata md = init_md();
  spec_info s;

  /* set default values for the parameter arrays */
  fwhm_inp[0] = 0.5;
  fwhm_inp[1] = 0.5;

  state_t[0] = 3;
  for (j=1; j<4; j++) {
    state_t[j] = 0.001;
  }
  n_t = 3;
  state_t[0] = n_t;

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

    switch(argv[1][1]){

    case 'h' :
      printf("Welcome to the yet-to-be implemented help message.\n");
      exit(1);
      break;

    case 'i' :
      /* printf( "processing input file: " ); */

      k = 0;
      /* extract the file name */
      for (j=3;argv[1][j] != '\0'; j++) {
        if (argv[1][j]    == '/') {
          k                = j+1;
        }
      }

      inp_fn = malloc(j-k+1);

      for (l = 0,j=k; argv[1][j] != '\0'; j++){

        if (argv[1][j] == '.') {
          j--;
          break;
        }

        inp_fn[l] = argv[1][j];
        l++;
      }

      for (j = 3; argv[1][j] != '\0'; j++) {
        input_sbuff[j-3] = argv[1][j];
        /* printf( "%c", input_sbuff[j-3]); */
      }

      /* printf( "\n" ); */
      len_fn     = j-3;
      inpath = malloc(len_fn+1);

      for (j          = 0; j<len_fn; j++) {
        inpath[j] = input_sbuff[j];
        /* printf( "%c", inpath[j]); */
      }
      /* printf( "\n" ); */
      inpath[len_fn] = '\0';

      /* set the inp_sfx  */
      j = len_fn;
      while(inpath[--j] != '.'){};

      inp_sfx = malloc(len_fn-j+1);

      for (k=0; j<=len_fn; j++) {
        inp_sfx[k] = inpath[j];
        k++;
      }

      break;

    case 'o' :
      /* printf( "processing input file: " ); */
      out_set = 1;

      /* extract the output path */
      for (j = 3; argv[1][j] != '\0'; j++) {
        input_sbuff[j-3] = argv[1][j];
        /* printf( "%c", input_sbuff[j-3]); */
      }

      len_op = j-3;
      outpath = malloc(len_op+1);

      for (j          = 0; j<len_op; j++) {
        outpath[j] = input_sbuff[j];
      }
      /* printf( "\n" ); */
      outpath[len_op] = '\0';
      printf( "  - outpath specified  = %s\n",outpath );
      break;

    case 'r' :
      n_res = 0;

      for (j = 3,k=0; argv[1][j] != '\0'; j++) {

        input_sbuff[k++] = argv[1][j];
        /* if ((argv[1][j] == ',') || ((argv[1][j+1] == '\0')||(argv[1][j+1] == '-'))) { */

        if ((argv[1][j] == ',') || (argv[1][j+1] == '\0')) {
          /* k--; /\* we dont want the comma or null sent to atof *\/ */
          num_buf = malloc(k+1);

          for (l = 0; l<k; l++) {
            num_buf[l] = input_sbuff[l];
            /* printf( "%c", num_buf[l]); */
          }
          num_buf[k] = '\0';

          /* printf( "\n" ); */

          res[n_res] = atof(num_buf);
          free(num_buf);
          num_buf = NULL;
          n_res++;

          k = 0;
        }
      }

      break;

    case 'e' :

      n_er = 1;

      for (j = 3,k=0; argv[1][j] != '\0'; j++) {

        input_sbuff[k++] = argv[1][j];
        /* if ((argv[1][j] == ',') || ((argv[1][j+1] == '\0')||(argv[1][j+1] == '-'))) { */

        if ((argv[1][j] == ',') || (argv[1][j+1] == '\0')) {
          /* k--; /\* we dont want the comma or null sent to atof *\/ */
          num_buf = malloc(k+1);

          for (l = 0; l<k; l++) {
            num_buf[l] = input_sbuff[l];
            /* printf( "%c", num_buf[l]); */
          }
          num_buf[k] = '\0';

          /* printf( "\n" ); */

          state_er[n_er] = atof(num_buf);
          free(num_buf);
          num_buf        = NULL;
          n_er++;

          k = 0;
        }
      }
      state_er[0] = n_er-1;

      break;

    case 'f' :

      /* the user specified thresholds for each state type (ground, initial,
         final) to be used for screening the states when calculating the map */
      m = 0;
      for (j = 3,k=0; argv[1][j] != '\0'; j++) {

        input_sbuff[k++] = argv[1][j];
        /* if ((argv[1][j] == ',') || ((argv[1][j+1] == '\0')||(argv[1][j+1] == '-'))) { */

        if ((argv[1][j] == ',') || (argv[1][j+1] == '\0')) {
          /* k--; /\* we dont want the comma or null sent to atof *\/ */
          num_buf = malloc(k+1);

          for (l = 0; l<k; l++) {
            num_buf[l] = input_sbuff[l];
            /* printf( "%c", num_buf[l]); */
          }
          num_buf[k] = '\0';

          /* printf( "\n" ); */

          fwhm_inp[m++] = atof(num_buf);
          free(num_buf);
          num_buf        = NULL;
          /* if ((n_er > 4)) { */
          /*   /\* switch on the dipole+quadrupole flag  *\/ */
          /* break; */
          /* } */
          /* if ((n_er > 8)) { */
          /* break; */
          /* } */

          k = 0;
        }
      }
      break;

    case 't' :
      n_t = 1;

      for (j = 3,k=0; argv[1][j] != '\0'; j++) {

        input_sbuff[k++] = argv[1][j];
        /* if ((argv[1][j] == ',') || ((argv[1][j+1] == '\0')||(argv[1][j+1] == '-'))) { */

        if ((argv[1][j] == ',') || (argv[1][j+1] == '\0')) {
          /* k--; /\* we dont want the comma or null sent to atof *\/ */
          num_buf        = malloc(k+1);

          for (l = 0; l<k; l++) {
            num_buf[l] = input_sbuff[l];
            /* printf( "%c", num_buf[l]); */
          }
          num_buf[k] = '\0';

          /* printf( "\n" ); */

          state_t[n_t] = atof(num_buf);
          free(num_buf);
          num_buf      = NULL;

          if ((n_t++) > 3) {
            break;
          }

          k = 0;
        }
      }

      state_t[0] = n_t;
      break;

    case 'F' :

      etype = 1;
      break;

    default :
      printf("Can not process the flag \"-%c\". To read the help documentation, \
 call tau with the flag \"-h\". Program terminating.\n",(char)argv [1][1]);
      exit(1);
    }
    argv++;
    argv++;
    argc--;
    argc--;
  }

  if (!out_set ) {
    if (getcwd(curr_dir,sizeof(curr_dir)) != NULL) {
      printf( "  - output path not specified in input, it defaults to: %s\n", curr_dir);
    } else {
      fprintf(stderr, "smap, main.c: unable to obtain the current directory for use as default output path.  \n");
      printf( "program terminating due to the previous error.\n");
      exit(1);
    }
    for (j=0; curr_dir[j] != '\0'; j++){};

    outpath = malloc(j);
    for (k=0; k<=j; k++) {
      outpath[k] = curr_dir[k];
    }
    len_op = j-1;
  }

  /* check if the output directory exists. if not, create it */
  if (stat(outpath,&st) == -1) {
    printf( "      output path not existing. creating directory.\n");

    out_buf[0] = outpath[0];
    /* expect the output path to start with a backslashe */
    for (j=1; j<len_op; j++) {

      out_buf[j] = outpath[j];
      if (out_buf[j] == '/') {
        out_buf[j+1] = '\0';

        /* some directories in the outpath might exists, so check
         for that*/
        if (stat(out_buf,&st) == -1) {
          mkdir(out_buf,0777);
        }
      }
    }
  }

  /* set the input file size */
  stat(inpath,&st);
  sz_inp = (int)st.st_size;
  l = j;

  md -> outpath = outpath;
  md -> inpath = inpath;
  md -> inp_fn = inp_fn;
  md -> inp_sfx = inp_sfx;

  md -> state_er  = state_er;
  md -> state_t = state_t;
  md -> res = res;
  md -> fwhm = fwhm_inp;

  s = init_sinfo(md);

  if (len_fn == 0) {
   fprintf(stderr, "\n\Error: smap.c, main: you didnt provide the path to an \
input file.");
    printf( "program terminating due to the previous error.\n");
    exit(EXIT_FAILURE);
  }
  else {
    /* extract the needed data from the input */
    printf( "  - extracting input data ..");
    fflush(stdout);

    if (parse_input(md)) {
      fprintf(stderr, "smap.c, main: unable to parse the input data \
contained in %s.\n",inpath);
      printf( "program terminating due to the previous error.\n");
      exit(EXIT_FAILURE);
    } else {
      printf( "    done.\n" );
    }
  }

  printf( "\n executing smap with the following..\n\n" );
  printf( "  - data contained in the input file:\n    %s\n\n", inpath);
  printf( "  - threshold values:\n    " );

  for (j=1; j<n_t; j++) {
    printf( "%le, ", state_t[j]);
  }
  printf( "\n\n" );

  printf( "  - state energy intervals values:\n    " );

  for (j=1; j<n_er; j++) {
    printf( "%le, ", state_er[j]);
  }

  printf( "\n\n" );
  printf( " execution progress:\n\n");

  calc_smap_m(md);
  write_log(md);

  for (j=0; j<6; j++) {
    free(parsed_input[j]);
  }
  free(parsed_input);
  free(inp_sfx);
  free_md(md);
  free(input_sbuff);
  free(state_er);
  free(state_t);
  printf( "\n smap successfully executed.\n" );
  printf( " program terminating.\n\n" );

  return 0;
}
