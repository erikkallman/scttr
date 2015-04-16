#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include "smap.h"
#include "input_formats.h"
#include "parse_input.h" /* for input_data array, as well as other
                            calculation-specific variables */
#include "rmap_cfg.h"
#define BUF_SIZE 256

int
main (int argc, char * argv[]) {
  int j,k,l; /* iteration variables */
  int len_fn;
  int dbg_flag;
  /* arrays for storing input file name data */
  char * input_sbuff = malloc(BUF_SIZE);
  char * fn_infile;
  char * method;
  char * num_buf;

  /* threshold values for the three different state types */
  double state_t[3] = {0.005, 0.00001, 0.000001};

  /* number of threshold values processed from command line */
  int n_t = 0;

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
      printf( "processing input file: " );

      for (j=3; argv[1][j] != '\0'; j++) {
        input_sbuff[j-3] = argv[1][j];
        printf( "%c", input_sbuff[j-3]);
      }

      printf( "\n" );
      len_fn = j-3;
      fn_infile = malloc(len_fn+1);

      for (j=0; j<len_fn; j++) {
        fn_infile[j] = input_sbuff[j];
      }

      fn_infile[len_fn] = '\0';

      /* extract the needed data from the input */
      if (parse_input(fn_infile, len_fn+1)) {
        fprintf(stderr, "rmap.c, main: unable to parse the input data \
contained in %s.\n",fn_infile);
        printf( "program terminating due to the previous error.\n");
        exit(EXIT_FAILURE);
      }

      break;

    case 'm' :
      /* the user specified a method to be used for calculating the scattering map */
      for (j=3; argv[1][j] != '\0'; j++) {
        input_sbuff[j-3] = argv[1][j];
      }
      argv[1] + j;

      method = malloc(len_fn+1);

      for (j=0; j<len_fn; j++) {
        method[j] = input_sbuff[j];
      }
      method[len_fn] = '\0';
      printf( "got the method: %s\n", method);

      break;

    case 't' :
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
          n_t++;
          k = 0;
        }
      }

      break;

    default :
      printf("Cannot process the flag \"-%c\". To read the help documentation, \
 call tau with the flag \"-h\". Program terminating.\n",(char)argv [1][1]);
      exit(1);
    }
    argv++;
    argv++;
    argc--;
    argc--;
  }

  printf( "\nexecuting rmap with the following..\n\n" );
  printf( "  - data contained in the input file:\n    %s\n\n", fn_infile);
  printf( "  - threshold values:\n    " );

  for (j=0; j<3; j++) {
    printf( "%le, ", state_t[j]);
  }
  printf( "\n\n" );
  printf( "  - method:\n    %s\n", method);

  printf( "\n\n" );
  printf( "execution progress:\n\n");


  calc_smap_m(method, fn_infile, \
  screen_states(fn_infile, 3, state_t[0], state_t[1], state_t[2]));

  if (dbg_flag) {
    calc_smap_dbg(method, fn_infile,\
    screen_states(fn_infile, 3, state_t[0], state_t[1], state_t[2]));
  }

  free(method);
  free(fn_infile);
  free(input_sbuff);
  printf( "\nrmap successfully executed.\n" );
  printf( "program terminating.\n\n" );

  return 0;
}
