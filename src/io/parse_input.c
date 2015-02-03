#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "rmap.h"
#include "parse_input.h"
#include "input_formats.h"
#define BUF_SIZE 256

double * input_data[4];
static FILE * fp_infile; /* pointer to the opened input data file */

static int
parse_input_molcas (char * fn_infile) {

  int j,k,l,m; /* control loop variables */
  int n_matches;
  int n_lookup_str = 4; /* number of strings used for searching the input file */
  int a[4] = {200,300,200,300};  /* place-holder integers for allocation of memory to input_data*/

  int matches[4]; /* place to store the indexes of the lines containing the
                   matches */
  FILE * fp_infile;

  /* we are looking for four strings in the output file: one at the beginning
   of each data block, and one at the end. */
  char str_buf[BUF_SIZE];
  char * lookup_str[n_lookup_str];

  char s1[50] = "        MATRIX ELEMENTS OVER SPIN EIGENSTATES FOR:";
  char s2[45] = "  HAMILTONIAN MATRIX FOR THE ORIGINAL STATES:";
  char s3[35] = "         To  From     Osc. strength";
  char s4[37] = "         ----------------------------";

  /* store the string in a more easily accessed format */
  lookup_str[0] = s1;
  lookup_str[1] = s2;
  lookup_str[2] = s3;
  lookup_str[3] = s4;

  /* open the input file */
  if((fp_infile = fopen(fn_infile, "rt")) == NULL) {
    fprintf(stderr,"parse_input: unable to open the input file %s.\n",fn_infile);
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  k = 0; /* index for tmp_string */
  l = 0; /* index for lookup string */
  m = 0; /* index for string matches */
  while(n_matches > 0){
  /* read char by char */
    for (j=0; (c == fgetc(fp_infile)) != EOF; j++) {
      str_buf[k++] = c;
      /* until an entire line has been stored in the temporary str_buf buffer */
      if (c == '\0') {
        k = 0;
        /* check every line for a matching substring */
        if (strstr(str_buf,lookup_str[l++])) {
          matches[m++] = j;
        }
      }
    }
  }

  for (j=0; j<=3; j++) {
    if(( input_data[j] = malloc(a[j]*sizeof(double*))) == NULL ){
      fprintf(stderr, "parse_input_molcas, malloc: failed to allocate memory for \"input_data\"\n");
      printf( "program terminating due to the previous error.\n");
      exit(1);
    }
  }

  printf( "processed molcas2.\n" );
  return 0;
}

int
parse_input (char * fn_infile, /* name of input file */
             int len_fn){

  int j, k, l; /* looping variables */
  FILE * fp_infile;
  char format[BUF_SIZE] = {0};



  /* loop over the input file name and extract the ending */
  for (j=len_fn; j>0; j--) {
    format[j] = fn_infile[j-3];
    if (format[j] = '.') {

      if (strcmp(format,MOLCAS_FORMAT)){
        printf( "found molcas.\n" );
        parse_input_molcas(fn_infile);
        printf( "processed molcas3\n" );
        break; /* for */
      }

      else {
        printf( "format not found in input_formats.h\n" );
      }
    }
  }

  printf( "file processed\n" );
  return 0;
}
