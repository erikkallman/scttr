#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include "input_formats.h"
#include "parse_input.h" /* for input_data array, as well as other
                            calculation-specific variables */
#include "rmap_cfg.h"
#define BUF_SIZE 256


int
main (int argc, char * argv[]) {
  int j,k,l; /* iteration variables */
  int len_fn;
  /* arrays for storing input file name data */
  char * fn_infile_buff = malloc(BUF_SIZE);
  char * fn_infile;
  /* char fn_infile[BUF_SIZE] = {0}; */

  /* process the input arguments */
  if (argc == 1) {
    printf("No command line arguments provided. To read the help documentation, provide the argument \"-h\". Program terminating.\n");
    fflush(stdout);
  }

  while (argc > 1 && (argv[1][0] == '-')) {

    switch(argv[1][1]){

    case 'h' :
      printf("some help was requested. \n");
      break;

    case 'i' :
      printf( "processing input file.. \n" );

      for (j=3; argv[1][j] != '\0'; j++) {
        fn_infile_buff[j-3] = argv[1][j];
      }

      len_fn = j-3;
      fn_infile = malloc(len_fn+1);

      for (j=0; j<len_fn; j++) {
        fn_infile[j] = fn_infile_buff[j];
      }

      fn_infile[len_fn] = '\0';

      /* extract the needed data from the input */
      printf( "%s\n", fn_infile);
      parse_input(fn_infile, len_fn+1);

      break;

    default :
      printf("Cannot process the flag \"-%c\". To read the help documentation, \
 call tau with the flag \"-h\". Program terminating.\n",(char)argv [1][1]);
      exit(1);
    }
    argc--;
    argv++;
  }

  /* for (j=0; j<=3; j++) { */
  /*     free(input_data[j]);/\* memory allocated in parse_input.c *\/ */
  /* } */
  /* free(input_data); */

  /* for (j=0; j<3; j++) { */
  /*   free(state_indices[j]); */
  /* } */
  /* free(state_indices); */

  /* free(state_indices); */
  free(fn_infile);
  free(fn_infile_buff);
  printf( "rmap successfully executed.\n" );
  return 0;
}
