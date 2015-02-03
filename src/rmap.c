#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "parse_input.h"
#include "input_formats.h"
#include "rmap_cfg.h"
#define BUF_SIZE 256

int
main (int argc, char * argv[]) {
  int j,k,l; /* iteration variables */
  int fn_len;

  /* arrays for storing input file name data */
  char fn_infile[BUF_SIZE] = {0};
  char format[BUF_SIZE] = {0};
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

      /* open the job file and read its data */
      for (j=3; argv[1][j] != '\0'; j++) {
        fn_infile[j-3] = argv[1][j];
      }
      fn_len = j-3;

      /* loop over the file name and extract the ending */
      for (j=fn_len; j>0; j--) {

        format[j] = fn_infile[j-3];
        if (format[j] = '.') {

          if (strcmp(format,MOLCAS_FORMAT)){
            printf( "found molcas.\n" );
          }
          else {
              printf( "format not found in input_formats.h\n" );
            }
        }
        break;                  /* for loop */
      }

      break;

    default :
      printf("Cannot process the flag \"-%c\". To read the help documentation, \
 call tau with the flag \"-h\". Program terminating.\n",(char)argv [1][1]);
      exit(1);
    }
    argc--;
    argv++;
  }
  if (parse_input(fn_infile, j-3)) {
    fprintf(stderr, "unable to process input the input in %s",fn_infile);
    exit(1);
  }

  return 0;
}
