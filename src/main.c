#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>
#include "state.h"
#include "std_char_ops.h"
#include "get_numsl.h"
#include "smap.h"
#include "input_formats.h"
#include "parse_input.h" /* for input_data array, as well as other
                            calculation-specific variables */
#include "smap_cfg.h"

#define BUF_SIZE 256
int sz_inp;

int
main (int argc, char * argv[]) {

  printf( "\n\n smap calculation initiated.\n\n" );
  int j,k,l;                    /* iteration variables */
  int n_t;                      /* number of numbers read from input flag */
  int n_er;                     /* number of numbers eigenstate energy values provided */
  int n_res;
  int len_op   = 0;
  int len_fn   = 0;
  int dbg_flag = 0;
  int out_set  = 0;

  int len_lf,len_df,len_pf,len_bf,len_tf;

  double * res         = malloc(2*sizeof(double));
  /* arrays for storing input file name data */
  char *   input_sbuff = malloc(BUF_SIZE);
  char *   fn_infile;
  char *   outpath;
  char *   fn_relpath;

  char method[5]       = "test";
  char *   num_buf;

  char format[BUF_SIZE] = {0};
  char curr_dir[BUF_SIZE] = {0};
  char out_buf[BUF_SIZE] = {0};

  /* strings of all output files from the program */
  char * dat_fpstr;
  char * plot_fpstr;
  char * log_fpstr;
  char * bin_fpstr;
  char * tmp_fpstr;

  char datformat[5] = ".dat";
  char plotformat[4] = ".gp";
  char logformat[5] = ".txt";
  char binformat[5] = ".bin";
  char tmpformat[5] = ".tmp";

  len_df = 4;
  len_pf = 3;
  len_lf = 4;
  len_bf = 4;
  len_tf = 4;

  /* double state_t[3] = {0.005, 0.001, 0.001}; */

  /* use first element to specify number of values stored in matrix */
  /* threshold values for the three or six different state types */
  double * state_t = malloc(7*sizeof(double));
  /* initial/final and intermediate state energy ranges */
  double * state_er = malloc(9*sizeof(double));

  state_t[0] = 0;
  state_er[0] = 0;

  struct stat st = {0};
    /* char fn_relpath[BUF_SIZE] = {0}; */

  /* process the input arguments */
  if (argc == 1) {
    printf("No command line arguments provided. To read the help documentation, provide the argument \"-h\". Program terminating.\n");
    fflush(stdout);
  }

  while (argc > 1 && (argv[1][0] == '-')) {

    switch(argv[1][1]){

    case 'd' :
      dbg_flag = 1;
      break;

    case 'h' :
      printf("some help was requested. \n");
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

      fn_infile = malloc(j-k+1);

      for (l = 0,j=k; argv[1][j] != '\0'; j++){

        if (argv[1][j] == '.') {
          j--;
          break;
        }

        fn_infile[l] = argv[1][j];
        l++;
      }

      for (j = 3; argv[1][j] != '\0'; j++) {
        input_sbuff[j-3] = argv[1][j];
        /* printf( "%c", input_sbuff[j-3]); */
      }

      /* printf( "\n" ); */
      len_fn     = j-3;
      fn_relpath = malloc(len_fn+1);

      for (j          = 0; j<len_fn; j++) {
        fn_relpath[j] = input_sbuff[j];
        /* printf( "%c", fn_relpath[j]); */
      }
      /* printf( "\n" ); */
      fn_relpath[len_fn] = '\0';

      /* set the format  */
      j = len_fn;
      while(fn_relpath[--j] != '.'){};

      for (k=0; j<=len_fn; j++) {
        format[k] = fn_relpath[j];
        k++;
      }
      format[k] = '\0';

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
      /* the user specified a method to be used for calculating the scattering map */
      for (j = 3,k=0; argv[1][j] != '\0'; j++) {

        input_sbuff[k++] = argv[1][j];
        /* if ((argv[1][j] == ',') || ((argv[1][j+1] == '\0')||(argv[1][j+1] == '-'))) { */

        if ((argv[1][j] == ',') || (argv[1][j+1] == '\0')) {
          /* k--; /\* we dont want the comma or null sent to atof *\/ */
          num_buf                                 = malloc(k+1);

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
          /* if ((n_res > 2)) { */
          /*   /\* switch on the dipole+quadrupole flag  *\/ */
          /*   break; */
          /* } */
          /* if ((n_res > 8)) { */
          /* break; */
          /* } */

          k = 0;
        }
      }

      break;

    case 'e' :

      n_er = 1;
      /* the user specified a method to be used for calculating the scattering map */
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
      state_er[0] = n_er-1;

      break;

    case 't' :
      n_t = 1;
      /* the user specified a method to be used for calculating the scattering map */
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
  stat(fn_relpath,&st);
  sz_inp = (int)st.st_size;

  /* the output path string length is known, build the output file strings */
  dat_fpstr = malloc(len_op+len_df+len_fn);
  plot_fpstr = malloc(len_op+len_pf+len_fn);
  log_fpstr = malloc(len_op+len_lf+len_fn);
  bin_fpstr = malloc(len_op+len_bf+len_fn);
  tmp_fpstr = malloc(len_op+len_tf+len_fn);

  for (j=0; outpath[j] != '\0'; j++) {
    dat_fpstr[j] = outpath[j];
    plot_fpstr[j] = outpath[j];
    log_fpstr[j] = outpath[j];
    bin_fpstr[j] = outpath[j];
    tmp_fpstr[j] = outpath[j];
  }

  /* dat_fpstr[j] = '/'; */
  /* plot_fpstr[j] = '/'; */
  /* log_fpstr[j] = '/'; */
  /* bin_fpstr[j] = '/'; */
  /* tmp_fpstr[j] = '/'; */
  /* j++; */

  /* for (k=0; j<len_op+len_infile; j++) { */
  for (k=0; fn_infile[k] != '\0'; j++,k++) {
    dat_fpstr[j] = fn_infile[k];
    plot_fpstr[j] = fn_infile[k];
    log_fpstr[j] = fn_infile[k];
    bin_fpstr[j] = fn_infile[k];
    tmp_fpstr[j] = fn_infile[k];
  }

  l = j;

  for (j=l, k=0; datformat[k] != '\0'; j++) {
    dat_fpstr[j] = datformat[k++];
  }

  for (j=l, k=0; datformat[k] != '\0'; j++) {
    plot_fpstr[j] = plotformat[k++];
  }

  for (j=l, k=0; logformat[k] != '\0'; j++) {
    log_fpstr[j] = logformat[k++];
  }

  for (j=l, k=0; binformat[k] != '\0'; j++) {
    bin_fpstr[j] = binformat[k++];
  }

  for (j=l, k=0; tmpformat[k] != '\0'; j++) {
    tmp_fpstr[j] = tmpformat[k++];
  }

  dat_fpstr[j] = '\0';
  plot_fpstr[j] = '\0';
  log_fpstr[j] = '\0';
  bin_fpstr[j] = '\0';
  tmp_fpstr[j] = '\0';

  /* printf( "%s\n",dat_fpstr ); */
  /* printf( "%s\n", log_fpstr); */
  /* printf( "%s\n",tmp_fpstr ); */
  /* printf( "%s\n",bin_fpstr ); */
  /* printf( "%s\n",plot_fpstr ); */

  if (len_fn == 0) {
   fprintf(stderr, "\n\Error: smap.c, main: you didnt provide the path to an\
 input file.");
    printf( "program terminating due to the previous error.\n");
    exit(EXIT_FAILURE);
  } else {
    /* extract the needed data from the input */
    printf( "  - extracting input data ..");
    fflush(stdout);

    if (parse_input(state_er, fn_relpath, tmp_fpstr, format, bin_fpstr, len_fn+1)) {
      fprintf(stderr, "smap.c, main: unable to parse the input data \
contained in %s.\n",fn_relpath);
      printf( "program terminating due to the previous error.\n");
      exit(EXIT_FAILURE);
    } else {
      printf( "    done.\n" );
    }
  }

  printf( "\n executing smap with the following..\n\n" );
  printf( "  - data contained in the input file:\n    %s\n\n", fn_relpath);
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

  calc_smap_m(fn_infile,dat_fpstr,plot_fpstr,state_er, state_t, res);
  write_log(state_er, state_t, res, fn_relpath, log_fpstr, 5);

  for (j=0; j<6; j++) {
    free(parsed_input[j]);
  }
  free(parsed_input);

  free(outpath);
  free(fn_relpath);
  free(fn_infile);
  free(input_sbuff);
  free(state_er);
  free(state_t);
  printf( "\n smap successfully executed.\n" );
  printf( " program terminating.\n\n" );

  return 0;
}
