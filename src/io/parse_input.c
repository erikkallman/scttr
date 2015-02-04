#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <regex.h>
#include <ctype.h>
#include "parse_input.h"
#include "input_formats.h"
#define BUF_SIZE 256
#define bin_flip(x) ((x) == 1 ? 0 : 1)

double * input_data[4];
static FILE * fp_infile; /* pointer to the opened input data file */
static char key[] = { ' ', '\n', '\t','-', 0 }; /* key for checking empty strings. */

static int
isempty (char * s,
         int len) {
  int j,k;

  for (j=0; j<len; j++) {
    if (strchr(key, s[j]) == NULL) {
      /* as soon as strchr returns something that isnt a match,
         of the key, the string is not empty */
      return 0;/* match found, it isnt empty */
    }
  }
  /* we got through the string without finding a match to key.
     the string is "empty" */
  return 1;
}

static int
isdigitin (char * s,
           int len) {
  int j;
  for (j=0; j<len; j++) {
    if (isdigit(s[j])){
      return 1;
    }
  }

  return 0;
}

static int
parse_input_molcas (char * fn_infile) {

  int j,k,l,m,i,n; /* control loop variables */
  int mode; /* string matching mode flag */
  int next_match;
  const int n_lookup_str = 4; /* number of strings used for searching the input file */
  const int n_maxmatch = n_lookup_str/2;
  int n_matchfound = 0;
  int match_vals[4] = {0}; /* place to store the indexes of the lines containing the
                              matches */

  /* place to store the indexes of the lines containing the
     matches. each data block will start after,  and end after the offset values.*/
  int match_ln[4];
  FILE * fp_infile;

  int c; /* temporary char for storing input file characters */
  /* we are looking for four strings in the output file: one at the beginning
     of each data block, and one at the end. */
  char * str_buf;

  regex_t c_regex_m1;
  regex_t c_regex_m2;
  regex_t c_regexs[2] = {c_regex_m1, c_regex_m2};

  str_buf = malloc(BUF_SIZE);

  const char s_regex_m1[6] = "[0-9]";
  const char s_regex_m2[6] = "[0-9]";
  const char * s_regexs[2] = {s_regex_m1, s_regex_m1};

  /* compile the regular expression my_regexp to check its compatibility */
  for (j=0; j<2; j++) {
    if (regcomp(&c_regexs[j], s_regexs[j], REG_EXTENDED) != 0) {
      fprintf(stderr, "parse_cfg: failed to compile the regexp pattern %s.\n",s_regexs[j]);
      exit(1);
    }
  }

  const char s1[26] = "        Relative EVac(au)";
  const char s2[36] = " Weights of the five most important";
  const char s3[36] = "         To  From     Osc. strength";
  const char s4[81] = " ###############################################################################";
  const char * lookup_str[4] = {s1,s2,s3,s4};

  /* open the input file */
  if((fp_infile = fopen(fn_infile, "rt")) == NULL) {
    fprintf(stderr,"parse_input: unable to open the input file %s.\n",fn_infile);
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  k = 0; /* index for tmp_string */
  l = 0; /* index for lookup string */
  m = 0; /* index for string matches */
  mode = 0; /* start of in string search mode */

  /* read char by char */
  for (j=0; ((c = fgetc(fp_infile)) != EOF); j++, k++) {
    str_buf[k] = (char) c;

    /* check to see if there are any matches left to catch in the input data */
    if (n_matchfound < n_maxmatch) {

      /* keep extracting characters from the input data until an entire line
         has been stored in the temporary str_buf buffer */
      if (str_buf[k] == '\n') {

        /* check every line for a matching substring */
        if (strstr(str_buf,lookup_str[l]) || (mode == 1)) {

          /* we found the first substring, the data we're looking for is
             inside the coming block of text. switch to mode 1.*/
          if (mode == 0) {
            match_ln[m++] = j;
            l++;
            mode = bin_flip(mode);
          }
          else { /* we're in count the matched lines mode */
            /* if (regexec(&c_regexs[m-1],str_buf,0,NULL,0)) { */
            if (isdigitin(str_buf,k-1) == 1) {
              match_vals[n_matchfound]++;
            }
            else if(isempty(str_buf,k) == 0) { /* mode = 1 and a line match without regexp match means
                                                  that we reached the end of this data block*/
              /* skip empty lines */
              match_ln[m++] = j; /* store the line number of the last delimiting
                                     string of the block */
              l++;
              mode = bin_flip(mode); /* switch back to mode 0 */
              n_matchfound++; /* one data block was found */
            }
          }
        }
        k = 0;
      }
    } else {
      break;
    }
  }

  free(str_buf);
  printf( "\n\n found the following matches\n" );
  for (j=0; j<n_lookup_str/2; j++) {
    printf( "%d:%d\n", j, match_vals[j]);
  }

  printf( "\n\n on the following offsets\n" );
  for (j=0; j<n_lookup_str; j++) {
    printf( "%d:%d\n", j, match_ln[j]);
  }

  for (j=0; j<=2; j++) {
    for (k=0; k<2; k++) {
      if(( input_data[j+k] = malloc(match_vals[j]*sizeof(double*))) == NULL ){
        fprintf(stderr, "parse_input_molcas, malloc: failed to allocate memory for \"input_data\"\n");
        printf( "program terminating due to the previous error.\n");
        exit(1);
      }
    }
  }


  l = 0; /* index for string matches */

  for (j=0; j<n_maxmatch; j++) {

    next_match = match_ln[j];
    fseek(fp_infile, next_match, 0);

    for (k=next_match; k<match_ln[j+1]; k++) {
      c = fgetc(fp_infile);
      str_buf[l++] = (char) c;
      if (str_buf[l] == '\n') {
        for (n=0; n<l; n++) {
          printf( "%c", str_buf[n]);
        }
        printf( "\n" );
        l=0;
      }
    }
  }


  /* the input data structure has been allocated. now read data from the input
     file, using the indexes in match_ln, and store it.  */

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
