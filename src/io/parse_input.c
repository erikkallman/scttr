#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <regex.h>
#include <ctype.h>
#include "parse_input.h"
#include "input_formats.h"
#define BUF_SIZE 256
#define bin_flip(x) ((x) == 1 ? 0 : 1)

int n_states;
int n_trans;
double * input_data[4];
FILE * fp_infile; /* pointer to the opened input data file */

int power(int base,int exp)
{
  int power;
  power=1;
  while(exp-- > 0)
    power *=base;

  return power;
}

double sci_atof(char s[]){
  double val,pow;
  int sign,i,esign,exp;
  int power(int base,int exp);

  for(i=0;isspace(s[i]);i++)
    ;

  sign=(s[i]=='-')?-1:1;

  if(s[i]=='+' || s[i] == '-')
    i++;

  for(val=0.0;isdigit(s[i]);i++)
    val = 10.0 * val + (s[i] - '0');

  if(s[i]=='.')
    i++;

  for(pow=1.0;isdigit(s[i]);i++)
    {
      val = 10.0 * val + (s[i] - '0');
      pow *= 10.0;
    }

  if(s[i]=='e' || s[i] =='E')
    i++;
  if(s[i]=='+' || s[i] =='-')
    {
      esign=s[i];
      i++;
    }

  for(exp=0;isdigit(s[i]);i++)
    exp=10.0 * exp + (s[i] - '0');

  if( esign == '-')
    return sign * (val / pow) / power(10,exp);
  else

    return sign * (val / pow) * power(10,exp);
}

char*
get_numinstr (const char * s,
              int idx,
              int str_len,
              int flag
              ){
  int k = 0;
  int l = 0;
  static int j;
  int n_digits_found = 0;

  static char num_key[] = {'-','.','E',};

  /* mode = 1 corresponds to reading digits, = 0, to reading anything else */
  int mode = 0;
  int last_mode = mode;
  char c;
  char last_c;
  char * numstr; /* a string containing the extracted number */
  char num_buf[BUF_SIZE] = {0};

  /* if we're searching this string more in the future, we can store the index
   of the last digit so that we wont have to loop over the entire string again */
  if (flag == 1) {
    j = 0;
  }

  /* printf( "getnuminstr s=%s", s); */
  /* printf( "str_len=%d", str_len); */

  for (; j<str_len; c = s[j]) {
    /* printf( "read char %c is? %d %d\n", c, isdigit(c), isalpha(c)); */
    /* printf( "last_char? %d", isalpha(last_c)); */
    /* check so that we're not reading a dashed line */
    /* while((last_c == '-') && (c  == '-')){ */
    /*   last_c = c; */
    /*   c = s[j++]; */
    /* } */

    if (((strchr(num_key,c) != NULL) || (isdigit(c) != 0)) &&\
        ((isalpha(last_c) == 0 ) || (last_c == 'E'))) { /* check if we're still reading a number and avoid dashed lines */

      /* check if we're reading the right digit, else ignore the result */
      if (k == idx) {
        num_buf[l++] = c;
      }

      if (mode == 0) {
        /* now proceeding to read a digit */
        /* printf( "flopped at digit %d, of len %d\n",k,l ); */
        last_mode = mode;
        mode = bin_flip(mode);
      }

      if (last_mode == 0) {

      }
    }
    /* if we read a non-number character and we're in read mode, we have
     read the end of a digit. */
    else if(mode == 1) {
      if (++k > idx) { /* increase the counter for read numbers */
        /* printf( "broke at digit %d, of len %d\n",k,l ); */
        break;
      }
      /* printf( "flipped at digit %d, of len %d\n",k,l ); */
      mode = bin_flip(mode); /* go back to reading non numbers */
    }

    /* if in reading number mode and found something non number
       the digit we were reading has now ended. increment k to
       note that we have read a number.
    */
    /* else if((mode == 1)  && (k > idx)){ /\* make sure we passed the last number */
    /*                                      we wanted to read *\/ */
    /*   printf( "broke at digit %d, of len %d\n",k,l ); */
    /*   /\* sleep(2) *\/; */
    /*   break; */
    /* } */
    last_c = c;
    j++;
  }

  numstr = malloc(j);
  /* printf( "\n\n========return\n" ); */
  for (k=0; k<=l; k++) {
    /* store the number and return a pointer to it.
       this gets freed up by the caller. */
    numstr[k] = num_buf[k];
    /* printf( "%c", numstr[k]); */
  }
  /* printf( "\n\n========return\n" ); */

  return numstr;
}

int
get_numsl (char * str,
           int * idxs_out, /* idexes of numbers in string that we want */
           int str_len,
           int n_idxs,
           ...){
  int j,k;
  va_list argv;
  char * numstr;

  double * tmp_num;
  /* store the memory address to the variadic input arguments so that we can
     assign them locally */
  va_start(argv, n_idxs);

  for (j=0; j<n_idxs; j++) { /* loop over indexes */
    /* printf( "\n\n==OOO index %d OOO ==\n\n", j); */

    tmp_num = va_arg(argv, double*); /* grab the next vararg */

    /* find the indexed number in the string and store it.  */
    /* printf( "\n==============================\n"); */
    /* printf( "get_numsl sent string=%s\n", str); */
    /* printf( "of len=%d\n", str_len); */

    numstr = get_numinstr(str,idxs_out[j],str_len,1);
    /* printf( "get_numsl got string %s\n",numstr ); */

    *tmp_num = sci_atof(numstr); /* extract the next memory location */
    /* printf( "resulting in num %le\n",*tmp_num); */
    /* printf( "================================\n\n" ); */
    /* sleep(1); */
    free(numstr);
    /* numstr = NULL; /\* prevent the pointer from becoming a dangling pointer *\/ */
  }

  va_end(argv);

  return 1;
}

int
isempty (char * s,
         int len) {
  int j,k;

  static char key[5] = { ' ', '\n', '\t','-', '\0' }; /* key for checking empty strings. */
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

int
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

double **
parse_input_molcas (char * fn_infile) {

  int j,k,l,m,i,n,k_its; /* control loop variables */
  int mode; /* string matching mode flag */
  int match_start;
  int match_end;
  int tmp_int;
  int idx_from;
  int idx_to;
  int extr_i;
  float extr_f;

  double ** parsed_input;
  double * e_eigval;
  double ** trans_idxs;
  double * t_mom;

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
  const char s2[51] = " Weights of the five most important spin-orbit-free";
  const char s3[18] = "         To  From";
  const char s4[52] = " ##################################################";
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

  /* we now know the number of states used in the molcas calculation, as well
     as the number of possible transitions.*/
  n_states = match_vals[0];
  n_trans = match_vals[1];

  printf( "%d %d\n", n_states,n_trans);

  if((parsed_input = malloc(5*sizeof(double*))) == NULL ){
    fprintf(stderr, "parse_input_molcas, malloc: failed to allocate memory for\
 \"input_data\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  /* allocate space for the "parsed input matrix" that will be filled with data
     in the remaining sections of this function */
  for (j=0; j<5; j++) {
      if((parsed_input[j] = malloc(n_trans*sizeof(double))) == NULL ){
        fprintf(stderr, "parse_input_molcas, malloc: failed to allocate memory \
for pointers in \"input_data\"\n");
        printf( "program terminating due to the previous error.\n");
        exit(1);
    }
  }

  /* storage for the energy eigenvalues */
  if((e_eigval = malloc(n_states*sizeof(double))) == NULL ){
    fprintf(stderr, "parse_input_molcas, malloc: failed to allocate memory for\
 \"e_eigval\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  /* storage for the transition moments */
  if((t_mom = malloc(n_trans*sizeof(double))) == NULL ){
    fprintf(stderr, "parse_input_molcas, malloc: failed to allocate memory for\
 \"e_eigval\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

    /* storage for the transition indexes, column 1 is from a state
     column 2 is to state index */
  if((trans_idxs = malloc(2*sizeof(double*))) == NULL ){
    fprintf(stderr, "parse_input_molcas, malloc: failed to allocate memory for\
 \"e_eigval\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  for (j=0; j<2; j++) {
    if((trans_idxs[j] = malloc(n_trans*sizeof(double))) == NULL ){
        fprintf(stderr, "parse_input_molcas, malloc: failed to allocate memory for\
 \"e_eigval\"\n");
        printf( "program terminating due to the previous error.\n");
        exit(1);
      }
  }

  /* int num_idxs[2] = {0,1}; */
  /* int n_idxs = 2; */

  int num_idxs1[1] = {1};
  int n_idxs1 = 1;

  int num_idxs2[3] = {0,1,8};
  int n_idxs2 = 3;

  for (j=0; j<n_lookup_str; j+=2) {

    l = 0; /* index for string matches */
    m = 0;
    /* if (j == 0) { /\* extract energy eigenvalues and state indexes *\/ */
    /*   printf( "\nextract energy eigenvalues and state indexes\n" ); */

    /* } */

    /* if (j == 2) { /\* extract transition moments and transition indexes *\/ */
    /*   printf( "\nextract transition moments and transition indexes\n" ); */

    /* } */

    match_start = match_ln[j];
    match_end = match_ln[j+1];

    fseek(fp_infile, match_start, 0);

    k_its = match_end-match_start;

    /* printf( "    seeking to %d to iterate for %d\n", match_start, k_its); */


    for (k=0; k<k_its; k++) {
      c = fgetc(fp_infile);
      str_buf[l] = (char)c;
      if ((str_buf[l] == '\n') && (l > 1)) { /* dont send blank lines */

        if (j == 0) { /* extract energy eigenvalues and state indexes */
          /* get_numsl(str_buf,num_idxs,l,n_idxs,&parsed_input[2][m],&parsed_input[3][m]); */
          get_numsl(str_buf,num_idxs1,l,n_idxs1,&e_eigval[m]);
          m++;
        }

        if (j == 2) { /* extract transition moments and transition indexes */
          /* printf( "strlen = %d\n",l); */
          /* sleep(2); */
          /* printf( "\n\nthe string: " ); */
          /* for (n=0; n<l; n++) { */
          /*   printf( "%c", str_buf[n]); */
          /* } */
          /* printf( "\n" ); */
          /* sleep(2); */
          get_numsl(str_buf,num_idxs2,l,n_idxs2,&trans_idxs[0][m],\
                    &trans_idxs[1][m],&t_mom[m]);
          m++;
        }
        l=0; /* reset the buffer write head to start reading a the next line */
      }
      l++;
    }
  }
  /* printf( "sorting data\n" ); */
  /* sort the data */
  for (j=0; j<n_trans; j++) {

    idx_from = trans_idxs[0][j];
    idx_to = trans_idxs[1][j];
    parsed_input[0][j] = idx_from;
    parsed_input[1][j] = idx_to;
    parsed_input[2][j] = e_eigval[idx_from-1];
    parsed_input[3][j] = e_eigval[idx_to-1];
    parsed_input[4][j] = t_mom[j];
  }

  for (k=0; k<n_trans; k++) {
    printf( "%le   %le   %le   %le   %le\n", parsed_input[0][k], parsed_input[1][k], parsed_input[2][k], parsed_input[3][k], parsed_input[4][k]);
  }

  /* the input data structure has been allocated. now read data from the input
     file, using the indexes in match_ln, and store it.  */
  fclose(fp_infile);
  free(str_buf);
  return parsed_input;
}

int
parse_input (char * fn_infile, /* name of input file */
             int len_fn){

  int j, k, l; /* looping variables */
  FILE * fp_infile;
  char format[BUF_SIZE] = {0};

  /* after having been used in a reference call to a parsing function the
     parsed_input array should contain an n_trans * 5 matrix loaded with
     data for each electronic transition that was calculated:
     [idx_f,idx_t,e_f,e_t,osc] where...
     idx_f = the index the electronic state FROM which the transition took place
     idx_t = the index the electronic state where the transition went TO
     e_f = energy corresponding to idx_f
     e_t = energy corresponding to idx_t
     mom = the transition moment for the transition
  */
  double ** parsed_input;

  /* loop over the input file name and extract the ending */
  for (j=len_fn; j>0; j--) {
    format[j] = fn_infile[j-3];
    if (format[j] = '.') {

      if (strcmp(format,MOLCAS_FORMAT)){
        printf( "found molcas.\n" );
        parsed_input = parse_input_molcas(fn_infile);
        printf( "processed molcas\n" );
        break; /* for */
      }

      else {
        printf( "format not found in input_formats.h\n" );
      }
    }
  }

  /* allocate memory space for all data structures used to store the input data. */
  /* for (j=0; j<=2; j++) { */
  /*   for (k=0; k<2; k++) { */
  /*     if((input_data[j+k] = malloc(match_vals[j]*sizeof(double*))) == NULL ){ */
  /*       fprintf(stderr, "parse_input_molcas, malloc: failed to allocate memory for \"input_data\"\n"); */
  /*       printf( "program terminating due to the previous error.\n"); */
  /*       exit(1); */
  /*     } */
  /*   } */
  /* } */

  for (j=0; j<5; j++) {
    free(parsed_input[j]);
  }
  free(parsed_input);
  printf( "file processed\n" );
  return 0;
}
