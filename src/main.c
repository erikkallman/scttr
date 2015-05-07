#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include "smap.h"
#include "input_formats.h"
#include "parse_input.h" /* for input_data array, as well as other
                            calculation-specific variables */
#include "smap_cfg.h"
#define BUF_SIZE 256

int
main (int argc, char * argv[]) {
  int j,k,l; /* iteration variables */
  int n_t; /* number of numbers read from input flag */
  int n_er; /* number of numbers eigenstate energy values provided */
  int n_mom;
  int len_fn = 0;
  int len_mn = 0;
  int dbg_flag;
  int * mom = malloc(4*sizeof(int));
  /* arrays for storing input file name data */
  char * input_sbuff = malloc(BUF_SIZE);
  char * fn_infile;
  char method[5] = "test";
  char * num_buf;


  /* double state_t[3] = {0.005, 0.001, 0.001}; */

  /* use first element to specify number of values stored in matrix */
  /* threshold values for the three or six different state types */
  double * state_t = malloc(7*sizeof(double));
  /* initial/final and intermediate state energy ranges */
  double * state_er = malloc(9*sizeof(double));

  state_t[0] = 0;
  state_er[0] = 0;

  /* char fn_infile[BUF_SIZE] = {0}; */

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

      for (j=3; argv[1][j] != '\0'; j++) {
        input_sbuff[j-3] = argv[1][j];
        /* printf( "%c", input_sbuff[j-3]); */
      }

      /* printf( "\n" ); */
      len_fn = j-3;
      fn_infile = malloc(len_fn+1);

      for (j=0; j<len_fn; j++) {
        fn_infile[j] = input_sbuff[j];
        /* printf( "%c", fn_infile[j]); */
      }
      /* printf( "\n" ); */
      fn_infile[len_fn] = '\0';

      break;

    /* case 'm' : */
    /*   /\* the user specified a method to be used for calculating the scattering map *\/ */
    /*   for (j=3; argv[1][j] != '\0'; j++) { */
    /*     input_sbuff[j-3] = argv[1][j]; */
    /*   } */
    /*   argv[1] + j; */
    /*   len_mn = j-3; */
    /*   method = malloc(len_mn+1); */

    /*   for (j=0; j<len_mn; j++) { */
    /*     method[j] = input_sbuff[j]; */
    /*   } */
    /*   method[len_mn] = '\0'; */
    /*   printf( "got the method: %s\n", method); */

    /*   break; */

    case 'm' :
      n_mom = 1;
      /* the user specified a method to be used for calculating the scattering map */
      for (j=3,k=0; argv[1][j] != '\0'; j++) {

        input_sbuff[k++] = argv[1][j];
        /* if ((argv[1][j] == ',') || ((argv[1][j+1] == '\0')||(argv[1][j+1] == '-'))) { */

        if ((argv[1][j] == ',') || (argv[1][j+1] == '\0')) {
          /* k--; /\* we dont want the comma or null sent to atof *\/ */
          num_buf = malloc(k+1);

          for (l=0; l<k; l++) {
            num_buf[l] = input_sbuff[l];
            /* printf( "%c", num_buf[l]); */
          }
          num_buf[k] = '\0';

          /* printf( "\n" ); */

          mom[n_mom] = atoi(num_buf);
          free(num_buf);
          num_buf = NULL;
          n_mom++;
          if ((n_mom > 2)) {
            /* switch on the dipole+quadrupole flag  */
            break;
          }
          if ((n_mom > 8)) {
            break;
          }

          k = 0;
        }
      }
      mom[0] = n_mom-1;

      break;

    case 'r' :

      n_er = 1;
      /* the user specified a method to be used for calculating the scattering map */
      for (j=3,k=0; argv[1][j] != '\0'; j++) {

        input_sbuff[k++] = argv[1][j];
        /* if ((argv[1][j] == ',') || ((argv[1][j+1] == '\0')||(argv[1][j+1] == '-'))) { */

        if ((argv[1][j] == ',') || (argv[1][j+1] == '\0')) {
          /* k--; /\* we dont want the comma or null sent to atof *\/ */
          num_buf = malloc(k+1);

          for (l=0; l<k; l++) {
            num_buf[l] = input_sbuff[l];
            /* printf( "%c", num_buf[l]); */
          }
          num_buf[k] = '\0';

          /* printf( "\n" ); */

          state_er[n_er] = atof(num_buf);
          free(num_buf);
          num_buf = NULL;
          n_er++;
          if ((n_er > 4)) {
            /* switch on the dipole+quadrupole flag  */
            break;
          }
          if ((n_er > 8)) {
            break;
          }

          k = 0;
        }
      }
      state_er[0] = n_er-1;

      break;

    case 't' :
      n_t = 1;
      /* the user specified a method to be used for calculating the scattering map */
      for (j=3,k=0; argv[1][j] != '\0'; j++) {

        input_sbuff[k++] = argv[1][j];
        /* if ((argv[1][j] == ',') || ((argv[1][j+1] == '\0')||(argv[1][j+1] == '-'))) { */

        if ((argv[1][j] == ',') || (argv[1][j+1] == '\0')) {
          /* k--; /\* we dont want the comma or null sent to atof *\/ */
          num_buf = malloc(k+1);

          for (l=0; l<k; l++) {
            num_buf[l] = input_sbuff[l];
            /* printf( "%c", num_buf[l]); */
          }
          num_buf[k] = '\0';

          /* printf( "\n" ); */

          state_t[n_t] = atof(num_buf);
          free(num_buf);
          num_buf = NULL;

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

  if (len_fn == 0) {
    fprintf(stderr, "\n\Error: smap.c, main: you didnt provide the path to an\
 input file.");
    printf( "program terminating due to the previous error.\n");
    exit(EXIT_FAILURE);
  } else {
    /* extract the needed data from the input */
    if (parse_input(state_er, mom, fn_infile, len_fn+1)) {
      fprintf(stderr, "smap.c, main: unable to parse the input data \
contained in %s.\n",fn_infile);
      printf( "program terminating due to the previous error.\n");
      exit(EXIT_FAILURE);
    }
  }

  printf( "\nexecuting smap with the following..\n\n" );
  printf( "  - data contained in the input file:\n    %s\n\n", fn_infile);
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

  printf( "  - transition moments:\n    " );

  for (j=1; j<n_mom; j++) {
    printf( "%d, ", mom[j]);
  }
  printf( "\n\n" );

  fprintf(stderr, "\n\n=======Valgrind eject point=======\n\n");
  exit(1);
  printf( "  - method:\n    %s\n", method);

  printf( "\n\n" );
  printf( "execution progress:\n\n");

  if (dbg_flag == 1) {
    calc_smap_dbg(method, fn_infile,\
                  screen_states(fn_infile, state_t, state_er));
  }

  calc_smap_m(method, fn_infile, \
              screen_states(fn_infile, state_t, state_er));


  free(fn_infile);
  free(input_sbuff);
  free(state_er);
  free(state_t);
  free(mom);
  printf( "\nsmap successfully executed.\n" );
  printf( "program terminating.\n\n" );

  return 0;
}
