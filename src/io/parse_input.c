#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <signal.h>
#include <string.h>
#include <fcntl.h>
#include <errno.h>
#include "std_char_ops.h"
#include "std_num_ops.h"
#include "get_numsl.h"
#include "parse_input.h"
#include "input_formats.h"
#include "sci_const.h"
#include "state.h"

int nt;
double ** parsed_input;
double e0;

/* function add_elastic

   * synopsis:

   * algorithm:

   * input:

   * output:

   * side-effects:

   */
int
add_elastic (double * state_er) {
  /* allocate memory for pi_el that is at most the size of ngs*nis*nfs*/
  /* copy the entire pi buffer to pi_el */
  /* make PI point to pi_el instead */
  /* deallocate the memory to PI */

  int j,l;
  int tmp_idx,next_to;
  int nt_el = nt;

  /* all elastic transitions handled so far */
  int * proc_st;

  double ** pi_buf;

  /* PI matrix with enough space to acommodate the elastic transitions */
  double ** pi_el;
  double ** pi_old;

  int n_proc,nb;

  if((pi_buf = malloc(6*sizeof(double*))) == NULL ){
    fprintf(stderr, "parse_input_molcas, malloc: failed to allocate memory for\
 \"input_data\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  for (j=0; j<6; j++) {
    if((pi_buf[j] = malloc(nt*sizeof(double))) == NULL ){
      fprintf(stderr, "parse_input_molcas, malloc: failed to allocate memory \
for pointers in \"input_data\"\n");
      printf( "program terminating due to the previous error.\n");
      exit(1);
    }
  }

  if((pi_el = malloc(6*sizeof(double*))) == NULL ){
    fprintf(stderr, "parse_input_molcas, malloc: failed to allocate memory for\
 \"input_data\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  for (j=0; j<6; j++) {
    if((pi_el[j] = malloc(2*nt*sizeof(double))) == NULL ){
      fprintf(stderr, "parse_input_molcas, malloc: failed to allocate memory \
for pointers in \"input_data\"\n");
      printf( "program terminating due to the previous error.\n");
      exit(1);
    }
  }

  /* free up the memory space occupied by the old matrix */
  if((proc_st = malloc(nt*sizeof(int))) == NULL ){
    fprintf(stderr, "parse_input.c , function add_elastic, malloc: failed to allocate memory for\
 \"proc_st\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  j = nb = n_proc = 0;
  while ((int)parsed_input[0][j] > 0) {
    /* printf( "state %d\n", (int)parsed_input[0][j]); */
    if (!intinint(proc_st, (int)parsed_input[1][j], n_proc) &&
        (((parsed_input[2][j]-e0)*AUTOEV >= state_er[1]) && ((parsed_input[2][j]-e0)*AUTOEV <= state_er[2])) &&
        (((parsed_input[3][j]-e0)*AUTOEV >= state_er[3]) && ((parsed_input[3][j]-e0)*AUTOEV <= state_er[4]))
        ) {

      next_to = (int)parsed_input[1][j];
      /* printf( "found elastic %d\n", next_to); */
      /* sleep(1); */
      /* found a "to" state that has not had its elastic transitions
       added yet. loop over the PI matrix and check if there are transitions
       fromn other states that need to be taken into account. */

      nb = 0;
      pi_buf[0][nb] = parsed_input[1][j];
      pi_buf[1][nb] = parsed_input[0][j];
      pi_buf[2][nb] = parsed_input[3][j];
      pi_buf[3][nb] = parsed_input[2][j];
      pi_buf[4][nb] = parsed_input[4][j];
      pi_buf[5][nb] = parsed_input[5][j];
      nb++;
      /* ^some un-rolling */

      l = get_inext((int)parsed_input[0][j]);
      while ((int)parsed_input[0][l] != -1) {
        if (((int)parsed_input[1][l] == next_to) &&
            ((parsed_input[2][l]-e0)*AUTOEV >= state_er[1]) && ((parsed_input[2][l]-e0)*AUTOEV <= state_er[2])){
          /* printf( "found from state %d\n", (int)parsed_input[1][l]); */
          /* transitions from another state */
          pi_buf[0][nb] = parsed_input[1][l];
          pi_buf[1][nb] = parsed_input[0][l];
          pi_buf[2][nb] = parsed_input[3][l];
          pi_buf[3][nb] = parsed_input[2][l];
          pi_buf[4][nb] = parsed_input[4][l];
          pi_buf[5][nb] = parsed_input[5][l];
          nb++;

          /* jump to the next "from" state, since any given state can only have
           one transition to another specific state */
          l = get_inext(parsed_input[1][j]);
        }
        else {
          /* printf( "didt find %d from state %d\n",next_to, (int)parsed_input[0][l]); */
          l++;
        }
      }

      tmp_idx = get_inext(next_to);
      proc_st[n_proc++] = next_to;

            /* printf( "appending .. " ); */
      /* append the data to the parsed_input matrix */
      for (l=0; l<nb; nt_el++,l++) {

        pi_el[0][nt_el] = pi_buf[0][l];
        pi_el[1][nt_el] = pi_buf[1][l];
        pi_el[2][nt_el] = pi_buf[2][l];
        pi_el[3][nt_el] = pi_buf[3][l];
        pi_el[4][nt_el] = pi_buf[4][l];
        pi_el[5][nt_el] = pi_buf[5][l];
      }
      /* printf( "done\n" ); */
      /* jump to the next state */
      j = get_inext((int)parsed_input[0][j]);
    }
    else {
      j++;
    }
  }
  fprintf(stderr, "\n\n=======Valgrind eject point=======\n\n");
  exit(1);
  /* copy the old pi data to pi_el */
    /* copy all of pi to pi_el */
  for (l=0; l<nt; nt_el++,l++) {
    pi_el[0][nt_el] = parsed_input[0][l];
    pi_el[1][nt_el] = parsed_input[1][l];
    pi_el[2][nt_el] = parsed_input[2][l];
    pi_el[3][nt_el] = parsed_input[3][l];
    pi_el[4][nt_el] = parsed_input[4][l];
    pi_el[5][nt_el] = parsed_input[5][l];
  }

  pi_el[0][nt_el] = -1;

  for (j=0; j<6; j++) {
    free(parsed_input[j]);
  }
  free(parsed_input);

  /* allocate new space for parsed_input, now that we know the total size */
  if((parsed_input = malloc(6*sizeof(double*))) == NULL ){
    fprintf(stderr, "parse_input_molcas, malloc: failed to allocate memory for\
 \"input_data\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  for (j=0; j<6; j++) {
    if((parsed_input[j] = malloc(nt_el*sizeof(double))) == NULL ){
      fprintf(stderr, "parse_input_molcas, malloc: failed to allocate memory \
for pointers in \"input_data\"\n");
      printf( "program terminating due to the previous error.\n");
      exit(1);
    }
  }

  for (l=0; l<6; l++) {
    for (j=0; j<nt_el; j++) {
      parsed_input[l][j] = pi_el[l][j];
    }
  }

  for (j=0; j<6; j++) {
    free(pi_buf[j]);
  }
  free(pi_buf);

  fprintf(stderr, "\n\n=======Valgrind eject point=======\n\n");
  exit(1);
  for (j=0; j<6; j++) {
    free(pi_el[j]);
  }
  free(pi_el);

  printf( "copied.\n" );
  nt = nt_el;
  pi2s();
  fprintf(stderr, "\n\n=======Valgrind eject point=======\n\n");
  exit(1);
  return 1;
}

int
parse_input_bin (char * bin_fpstr
                 ) {
  printf( "\n      processing binary file corresponding to the provided input file ..");
  fflush(stdout);
  int j; /* looping variables */

  int pi_xdim,pi_ydim;
  /* double ** parsed_input;
  /* test writing to binary */
  FILE * fp_bin = fopen(bin_fpstr,"rb");

  fread(&pi_ydim, sizeof(int),1,fp_bin);
  pi_xdim = 6;

  nt = pi_ydim;

  if((parsed_input = malloc(pi_xdim*sizeof(double*))) == NULL ){
    fprintf(stderr, "parse_input_bin, malloc: failed to allocate memory for\
 \"input_data\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  for (j=0; j<pi_xdim; j++) {
    if((parsed_input[j] = malloc((pi_ydim+1)*sizeof(double))) == NULL ){
      fprintf(stderr, "parse_input_bin, malloc: failed to allocate memory \
for pointers in \"input_data\"\n");
      printf( "program terminating due to the previous error.\n");
      exit(1);
    }
  }

  /* fp_bin = fopen(bin_fpstr,"rb"); */
  /* fseek(fp_bin,0,SEEK_SET); */
  /* fread(&nt, sizeof(int),1,fp_bin); */

  for (j=0; j<6; j++) {
    if (fread(parsed_input[j], sizeof(double), nt, fp_bin) != nt) {
      fprintf(stderr, "parse_input.c, function parse_input_bin: unable to read the PI matrix\n");
      printf( "program terminating due to the previous error.\n");
      fclose(fp_bin);
      exit(1);
    }
  }

  e0 = parsed_input[2][0];
  parsed_input[0][pi_ydim] = -1;

  if (fclose(fp_bin) != 0) {
    fprintf(stderr, "parse_input.c, function parse_input_bin: unable to close file:\n%s\n", bin_fpstr);
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }
  printf( " done\n");
  return 1;
}

int
parse_molout (char * fn_relpath,
              char * tmp_fpstr,
              int len_infile
              ) {
  printf( "\n      parsing the molcas .log file .. \n");
  int j,k,l,m,j_last; /* control loop variables */
  int mode; /* string matching mode flag */
  int string_flag = 0;

  FILE * fp_tmpdata;
  FILE * fp_relpath;

  const int n_lookup_str = 6; /* number of strings used for searching the input file */

  int c; /* temporary char for storing input file characters */
  /* we are looking for four strings in the output file: one at the beginning
     of each data block, and one at the end. */
  char * str_buf = malloc(BUF_SIZE*5);

  const char s1[35] = "Eigenvalues of complex Hamiltonian";
  const char s2[40] = "Dipole transition strengths (SO states)";
  const char s3[44] = "Quadrupole transition strengths (SO states)";

  /* create a pointer to the three data block beginners s1,s2,s3 */
  const char * lookup_str[3] = {s1,s2,s3};

  /* open the input file */
  if((fp_relpath = fopen(fn_relpath, "r")) == NULL) {
    fprintf(stderr,"parse_input: unable to open the input file %s.\n",fn_relpath);
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }
  if((fp_tmpdata = fopen(tmp_fpstr, "w+")) == NULL) {
    fprintf(stderr,"parse_input: unable to open the output file %s.\n",tmp_fpstr);
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  k    = 0;                     /* index for tmp_string */
  mode = 0;                     /* start of in string search mode */
  string_flag = 0;

  /* read the Molcas input file */
  for (j=0; ((c = fgetc(fp_relpath)) != EOF); j++, k++) {
    if (j % 100000 == 0) {
      printf( "        %.2f%%\r", (((float)(j*sizeof(char))/(float)sz_inp)*100));
    }

    str_buf[k]                        = (char)c;

    /* keep extracting characters from the input data until an entire line
       has been stored in the temporary str_buf buffer */
    if (str_buf[k] == '\n') {
      if (mode                     == 1) {

        /* check so that the string that was read only contains numbers */
        if ((isanyalpha(str_buf,k) == 0) &&
            (isdashes(str_buf,k)   == 0) &&
            (isempty(str_buf,k)    == 0)){

          if (string_flag == 0) {

            /* we are now reading numbers. accept no strings */
           string_flag =  bin_flip(string_flag);
            /* printf( "%s\n", str_buf); */
          }

          for (m = 0; m<=k; m++) {
            fputc(str_buf[m],fp_tmpdata);
          }
        }
        /* if we find a flag while string_flag ==1 and mode ==1, we jhave
         read beyond the table */
        else if(string_flag == 1){
          /* printf( "read last number\n" ); */
          string_flag = bin_flip(string_flag);

          /* something was read that was not a number flip back
             to searching for new line matches */
          mode    = bin_flip(mode);
          /* printf( "mode flip back!\n" ); */
          /* sleep(1); */
          /* fseek(fp_relpath, -1*(j-j_last),SEEK_CUR); */
        }
      }
      /* check every line for a matching substring */
      /* mode = 1 and a line match means that we reached
         the end of this data block */
      /* if ((l = strscmp(str_buf,lookup_str,6)) != -1) { */
      else if((l = strscmp(str_buf,lookup_str,3)) != -1) {

        /* we found the first substring, the data we're looking for is
           inside the coming table of numbers of text. switch to mode 1.*/
          fprintf(fp_tmpdata,"%s\n",lookup_str[l]);
          mode    = bin_flip(mode);
          /* printf( "mode flip!\n" ); */
          /* printf( "%s\n", lookup_str[l]); */
          /* sleep(1); */
      }
      k = 0;
      j_last = j;
    }
  }
  if ((fclose(fp_tmpdata) != 0) &&\
      (fclose(fp_relpath) != 0)) {
    fprintf(stderr, "parse_input.c, function parse_molout: unable to close files:\n%s\n%s\n", tmp_fpstr,fn_relpath);
    printf("program terminating due to the previous error.\n");
    exit(1);
  }

  printf( "\n      done. \n");
}

/* function

   * synopsis:

   * algorithm:

   * input:

   * output:

   * side-effects:

   */
int
ISINSIDE(double v,
         double r1,
         double r2) {

  if ((v>=r1) && (r2>=v)) {
    return 1;
  } else {
    return 0;
  }
}


#define BUF_SIZE 256

int nt;
int * idxs_map;
double ** parsed_input;
double tmax_d, tmax_q, e0;

int
check_pi (){
  printf( "\n      checking the integrity of the input matrix ..\n");
  int j = 0;

  /* check so that every state in PI can be reached with teh get_i function */
  while (1) {
    printf( "        %.2f%%\r", (((float)j/(float)nt)*100));
    if ((int)parsed_input[0][j] != (int)parsed_input[0][get_i((int)parsed_input[0][j])]) {
      return 0;
    }

    j = get_inext((int)parsed_input[0][j]);
    /* printf( "%d\n", parsed_input[0][j]); */
    if ((int)parsed_input[0][j] == -1) {
      break;
    }
  }
  printf( "\n      done.\n");
  return 1;
}

void
set_tmax(){
  int j = 0;
  tmax_d = tmax_q = -1;

  while (1) {
    if (parsed_input[5][j] == 1) {
      if ( parsed_input[4][j] > tmax_d) {
        tmax_d = parsed_input[4][j];
      }
    }
    else{
      if (parsed_input[4][j] > tmax_q) {
        tmax_q = parsed_input[4][j];
      }
    }
    j++;
    if (parsed_input[0][j] == -1) {
      break;
    }
  }
}

void
state2s(int idx){
  int j = get_i(idx);

  fprintf(stderr, "\n\n================================\n");
  fprintf(stderr, "=======content of state %d=======\n",idx);
  while ((int)parsed_input[0][j] == idx) {
    printf( "\n%le   %le   %le   %le   %le  %le\n",
            parsed_input[0][j],
            parsed_input[1][j],
            parsed_input[2][j],
            parsed_input[3][j],
            parsed_input[4][j],
            parsed_input[5][j]);
    /* printf( "index data %d %d %d\n", idxs_map[last_i-1], idxs_map[last_i-1]+j,(int)parsed_input[0][idxs_map[last_i-1]+j]); */
    j++;
  }
  fprintf(stderr, "\n\n=======content of state %d=======\n",idx);
  fprintf(stderr,     "================================\n\n");

}

void
pi2s () {
  int j = 0;
  fprintf(stderr, "\n\n===========================\n");
  fprintf(stderr, "=======content of PI=======\n\n");
  while ((int)parsed_input[0][j++] != -1) {
    printf( "\n%le   %le   %le   %le   %le  %le\n",
            parsed_input[0][j],
            parsed_input[1][j],
            parsed_input[2][j],
            parsed_input[3][j],
            parsed_input[4][j],
            parsed_input[5][j]);
    j++;
  }
  fprintf(stderr, "\n\n=======content of PI=======\n");
  fprintf(stderr, "===========================\n\n");

}

double
get_efrom (int from) {

  int j = 0;
  while (1) {
    if ((int)parsed_input[0][j] == from) {
      return parsed_input[2][j];
    }
    j++;
    if ((int)parsed_input[0][j] == -1) {
      break;
    }
  }
  fprintf(stderr, "\n\nget_t:ERROR\n");
  printf( "program terminating due to the previous error.\n");
  fprintf(stderr, "\n\n=======Valgrind eject point=======\n\n");
  exit(1);
}

double
get_eto (int to) {

  int j = 0;
  while (1) {
    if ((int)parsed_input[1][j] == to) {
      return parsed_input[3][j];
    }
    j++;
    if (parsed_input[0][j] == -1) {
      break;
    }
  }

  fprintf(stderr, "\n\nget_t:ERROR\n");
  printf( "program terminating due to the previous error.\n");
  fprintf(stderr, "\n\n=======Valgrind eject point=======\n\n");
  exit(1);
}

double
get_t (int from,
       int to) {

  int j = 0;

  while (1) {
    if (((int)parsed_input[0][j] == from) &&
      ((int)parsed_input[1][j] == to)){
      return parsed_input[4][j];
    }
    j++;
    if ((int)parsed_input[0][j] == -1) {
      break;
    }
  }

  fprintf(stderr, "\n\nget_t:ERROR\n");
  printf( "program terminating due to the previous error.\n");
  fprintf(stderr, "\n\n=======Valgrind eject point=======\n\n");
  exit(1);
}

int
get_i (int from
       ) {
  int last_i = (int)parsed_input[0][0];
  int j = 0;
  while (last_i != -1) {

    if ((int)parsed_input[0][j] == from) {
      return j;
    }
    j++;
    last_i = (int)parsed_input[0][j];
  }
  /* fprintf(stderr, "\n\nget_i:ERROR, tried to get state %d\n", from); */
  /* printf( "program terminating due to the previous error.\n"); */
  /* fprintf(stderr, "\n\n=======Valgrind eject point=======\n\n"); */
  /* exit(1); */
  return -1;
}

int
get_inext (int from
           ) {

  int last_i = (int)parsed_input[0][0];
  int j      = 0;
  while (last_i != -1) {

    if ((int)parsed_input[0][j] == from) {
      break;
    }
    j++;
    last_i = (int)parsed_input[0][j];
  }

  if ((int)parsed_input[0][j] != from) {
    fprintf(stderr, "\n\nget_inext:ERROR, cant get state %d, last_i = %d\n", from, last_i);
    printf( "program terminating due to the previous error.\n");
    fprintf(stderr, "\n\n=======Valgrind eject point=======\n\n");
    exit(1);
  }


  while(1){
    if ((int)parsed_input[0][j] != from) {
      return j;
    }
    j++;
  }

  fprintf(stderr, "\n\nget_inext:ERROR\n");
  printf( "program terminating due to the previous error.\n");
  fprintf(stderr, "\n\n=======Valgrind eject point=======\n\n");
  exit(1);
}

int
parse_input_tmp (double * state_er,
                 char * fn_tmpdata,
                 char * bin_fpstr
                 ) {
  printf( "\n      parsing the tmp file .. \n");
  int j,k,l,m,i,n,k_its,j_test; /* control loop variables */
  int mode; /* string matching mode flag */
  int match_start;
  int match_end;
  int idx_from;
  int idx_to;
  int tmp_idx2;
  int n_states;
  int n_trans;
  int last_i;

  int n_proc;

  int trs_type;

  int * num_idxs1;
  int * num_idxs2;
  int * trs_types;
  int * proc_idx;

  double tmp_idx = 0;
  double maxr    = get_maxl(state_er,state_er[0]);

  double *  e_eigval;
  int *     idxs_eigval;
  double ** pi_buf;
  double ** trans_idxs;
  double *  t_mom;

  FILE * fp_tmpdata;
  FILE * fp_bin;

  int c; /* temporary char for storing input file characters */
  /* we are looking for four strings in the output file: one at the beginning
     of each data block, and one at the end. */
  char * str_buf = malloc(BUF_SIZE*2);

  const char s1[38] = "Eigenvalues of complex Hamiltonian";
  const char s2[40] = "Dipole transition strengths (SO states)";
  const char s3[44] = "Quadrupole transition strengths (SO states)";

  /* create a pointer to the three data block beginners s1,s3, */
  const char * lookup_str[3] = {s1,s2,s3};

  /* open the input file */
  if((fp_tmpdata = fopen(fn_tmpdata, "r")) == NULL) {
    fprintf(stderr,"parse_input: unable to open the input file %s.\n",fn_tmpdata);
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  k = 0; /* index for tmp_string */
  l = 0; /* index for lookup string */

  /* storage for the energy eigenvalues */
  if((e_eigval = malloc(5000000*sizeof(double))) == NULL ){
    fprintf(stderr, "parse_input_molcas, malloc: failed to allocate memory for\
 \"e_eigval\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }
  if((trs_types = malloc(5000000*sizeof(double))) == NULL ){
    fprintf(stderr, "parse_input_molcas, malloc: failed to allocate memory for\
 \"e_eigval\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }
  if((idxs_eigval = malloc(5000000*sizeof(int))) == NULL ){
    fprintf(stderr, "parse_input_molcas, malloc: failed to allocate memory for\
 \"idxs_eigval\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }
  for (j=0; j<5000000; j++) {
    idxs_eigval[j] = -1;
  }
  /* storage for the transition moments */
  if((t_mom = malloc(5000000*sizeof(double))) == NULL ){
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
    if((trans_idxs[j] = malloc(5000000*sizeof(double))) == NULL ){
      fprintf(stderr, "parse_input_molcas, malloc: failed to allocate memory for\
 \"e_eigval\"\n");
      printf( "program terminating due to the previous error.\n");
      exit(1);
    }
  }

  if((num_idxs1 = malloc(2*sizeof(int))) == NULL ){
    fprintf(stderr, "parse_input_molcas, malloc: failed to allocate memory for\
 \"e_eigval\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  if((num_idxs2 = malloc(3*sizeof(int))) == NULL ){
    fprintf(stderr, "parse_input_molcas, malloc: failed to allocate memory for\
 \"e_eigval\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  num_idxs1[0] = 0;
  num_idxs1[1] = 1;
  int n_idxs1  = 2;

  num_idxs2[0] = 0;
  num_idxs2[1] = 1;
  num_idxs2[2] = 2;
  int n_idxs2  = 3;

  j = j_test = l = m = 0;

  n_states    = 0;
  n_trans     = 0;
  trs_type    = 1;

  /* now that data structures of the rights size has memory allocated for
     them, start reading data from the temporary file */
  /* rewind(fp_tmpdata); */
  while ((c = fgetc(fp_tmpdata)) != EOF) {

    str_buf[l] = (char)c;

    if ((str_buf[l]  == '\n') && (l > 0)) { /* dont send blank lines */
      str_buf[l+1] = '\0';

      /* if (strcmp(DAT_DELIM, str_buf) <= 0) { */
      if ((j_test = strscmp(str_buf, lookup_str,3)) != -1) {

        j = j_test;
        trs_type = j;

        /* printf( "Found the delimiter j = trs_type = %d\n",j ); */
        /* sleep(1); */
      }
      else{

        if ((j == 0) && (isempty(str_buf,l) != 1)) {

          /* extract energy eigenvalues and state indexes */
          get_numsl(str_buf,num_idxs1,l,n_idxs1,&tmp_idx,&e_eigval[n_states]);

          /* printf( "pre%d\n", j); */
          if ((fabs(e_eigval[n_states]-e_eigval[0])*AUTOEV) \
              < (maxr + 100 )) {
            /* printf( "e_eigval[%d]  = %le, %d\n", n_states+1, e_eigval[n_states],tmp_idx ); */
            idxs_eigval[(int)tmp_idx] = n_states;
            n_states++;

          }
          /* printf( "post%d\n", j); */
        }
        else if ((j > 0) && (isempty(str_buf,l) != 1)) {

          /* extract transition moments and transition indexes */
          get_numsl(str_buf,num_idxs2,l,n_idxs2,&trans_idxs[0][n_trans],\
                    &trans_idxs[1][n_trans],&t_mom[n_trans]);
          trs_types[n_trans] = trs_type;

          if (!((idxs_eigval[(int)(trans_idxs[0][n_trans])] == -1) || (idxs_eigval[(int)(trans_idxs[1][n_trans])] == -1)) &&
              (((fabs(e_eigval[idxs_eigval[(int)(trans_idxs[1][n_trans])]]-e_eigval[0])*AUTOEV) \
                < maxr) &&                \
               ((fabs(e_eigval[idxs_eigval[(int)(trans_idxs[0][n_trans])]]-e_eigval[0])*AUTOEV) \
                < maxr)))
            {
              /* printf( "%le %le %le\n", trans_idxs[0][n_trans], trans_idxs[1][n_trans], t_mom[n_trans]); */
              n_trans++;
          }
        }
      }

      /* reset the buffer write head to start reading a the next line */
      l = 0;
    }
    else{
      l++;
    }
  }


  /* fprintf(stderr, "\n\n=======Valgrind eject point=======\n\n"); */
  /* exit(1); */
  /* allocate space for the "parsed input matrix" that will be filled with data
     in the remaining sections of this function */

  if((proc_idx = malloc(n_states*sizeof(int))) == NULL ){
    fprintf(stderr, "parse_input_molcas, malloc: failed to allocate memory for\
 \"proc_idx\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  if((idxs_map = malloc(n_states*sizeof(int))) == NULL ){
    fprintf(stderr, "parse_input_molcas, malloc: failed to allocate memory for\
 \"idxs_map\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  if((parsed_input = malloc(6*sizeof(double*))) == NULL ){
    fprintf(stderr, "parse_input_molcas, malloc: failed to allocate memory for\
 \"input_data\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  for (j=0; j<6; j++) {
    if((parsed_input[j] = malloc((n_trans+1)*sizeof(double))) == NULL ){
      fprintf(stderr, "parse_input_molcas, malloc: failed to allocate memory \
for pointers in \"input_data\"\n");
      printf( "program terminating due to the previous error.\n");
      exit(1);
    }
  }
  parsed_input[0][n_trans] = -1; /* end of list */

  if((pi_buf = malloc(6*sizeof(double*))) == NULL ){
    fprintf(stderr, "parse_input_molcas, malloc: failed to allocate memory for\
 \"input_data\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  for (j=0; j<6; j++) {
    if((pi_buf[j] = malloc(n_trans*sizeof(double))) == NULL ){
      fprintf(stderr, "parse_input_molcas, malloc: failed to allocate memory \
for pointers in \"input_data\"\n");
      printf( "program terminating due to the previous error.\n");
      exit(1);
    }
  }

  /* finally, store the data in the parsed_input matrix for the parse_input function  */
  n_proc      = 0;
  idxs_map[0] = 1;
  last_i = 1;
  j = k = l = 0;
  e0 = e_eigval[0];

  while (l<n_trans-1) {
    if (l % 5 == 0) {
      printf( "        %.2f%%\r", (((float)l/(float)n_trans)*100));
    }
    if ((j+l) == (n_trans)) {
      /* dont read beyond the last value in PI */
      /* for (k=0; k<j; k++) { */
        /* printf( "%le %le %le %le %le\n", pi_buf[0][k],pi_buf[1][k],pi_buf[2][k],pi_buf[3][k],pi_buf[4][k],pi_buf[5][k]); */
      /* } */
      /* break; */
      pi_buf[0][j] = -1;
      idx_from = -1;
    } else {
      if ((ISINSIDE((e_eigval[idxs_eigval[(int)trans_idxs[1][j+l]]]-e0)*AUTOEV,state_er[3],state_er[4]) &&
           (ISINSIDE((e_eigval[idxs_eigval[(int)trans_idxs[0][j+l]]]-e0)*AUTOEV,state_er[5],state_er[6]))) &&
          (trs_types[j+l] == 1)
          ) {
        idx_from = trans_idxs[1][j+l];
        idx_to = trans_idxs[0][j+l];
        /* printf( "%d %d %d\n", (int)trans_idxs[0][j+l], (int)trans_idxs[1][j+l], (int)trs_types[j+l]); */
        /* fprintf(stderr, "\n\n ======= Valgrind eject point=======\n\n"); */
        /* exit(1); */
      }
      else {
        idx_from = trans_idxs[0][j+l];
        idx_to   = trans_idxs[1][j+l];
      }
      pi_buf[0][j] = idx_from;
      pi_buf[1][j] = idx_to;
      pi_buf[2][j] = e_eigval[idxs_eigval[idx_from]];
      pi_buf[3][j] = e_eigval[idxs_eigval[idx_to]];
      pi_buf[4][j] = t_mom[j+l];
      pi_buf[5][j] = trs_types[j+l];
    }
    if (idx_from != last_i) {

      /* we have read all transitions for a state */
      /* check if the last_i has already been processed */
      if (intinint(proc_idx,last_i,n_proc) == 0) {

        m = l;
        while ((int)pi_buf[0][m-l] != idx_from) {
          parsed_input[0][m] = pi_buf[0][m-l];
          parsed_input[1][m] = pi_buf[1][m-l];
          parsed_input[2][m] = pi_buf[2][m-l];
          parsed_input[3][m] = pi_buf[3][m-l];
          parsed_input[4][m] = pi_buf[4][m-l];
          parsed_input[5][m] = pi_buf[5][m-l];
          m++;
        }
        proc_idx[n_proc++] = last_i;
      }
      else{
        tmp_idx2 = get_inext(last_i);
        fwdsplice(pi_buf,parsed_input,tmp_idx2,l,j,6);
      }

      l += j;
      j = 0;
    } else {
      j++;
    }
    last_i = idx_from;

  }

  parsed_input[0][l] = -1;

  nt = l;

  /* the only way the parse_input_tmp function can get called is if there is if
   no binary file present. we can therefore safely write to the binary file
  without checking if it already exists. */
  if ((fp_bin = fopen(bin_fpstr,"wb")) == NULL) {
    fprintf(stderr, "parse_input.c, function parse_input_bin: unable to open the binary file used to store the PI matrix: %s\n", bin_fpstr);
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  fwrite((const void*)&nt, sizeof(int), 1, fp_bin);
  fflush(fp_bin);
  fclose(fp_bin);

  if ((fp_bin = fopen(bin_fpstr,"ab")) == NULL) {
    fprintf(stderr, "parse_input.c, function parse_input_bin: unable to open the binary file used to store the PI matrix: %s\n", bin_fpstr);
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  for (j=0; j<6; j++) {
    if (fwrite(parsed_input[j], sizeof(double), nt, fp_bin) != nt) {
      fprintf(stderr, "parse_input.c, function parse_input_bin: unable to write the PI matrix\n");
      printf( "program terminating due to the previous error.\n");
      fclose(fp_bin);
      exit(1);
    }
    fflush(fp_bin);
  }

  if ((fclose(fp_bin) != 0) &&\
      (fclose(fp_tmpdata) != 0)) {
    fprintf(stderr, "parse_input.c, function parse_input_tmp: unable to close files:\n%s\n%s\n", bin_fpstr,fn_tmpdata);
    printf("program terminating due to the previous error.\n");
    exit(1);
  }

  free(str_buf);
  free(num_idxs1);
  free(num_idxs2);

  free(e_eigval);
  free(t_mom);

  free(idxs_map);
  free(proc_idx);

  for (j=0; j<6; j++) {
    free(pi_buf[j]);
  }
  free(pi_buf);

  for (j=0; j<2; j++) {
    free(trans_idxs[j]);
  }
  free(trans_idxs);
  printf( "\n      done.\n");
  return EXIT_SUCCESS;

}

int
parse_input (double * state_er,
             char * fn_relpath, /* name of input file */
             char * tmp_fpstr,
             char * format,
             char * bin_fpstr,
             int len_fn){

  int j, k, l; /* looping variables */
  int n_states,n_trans;

  FILE * fp_tmp;
  FILE * fp_bin;
  /* fout = open("test.bin", O_CREAT | O_WRONLY | O_EXCL, S_IRUSR | S_IWUSR); */
  fp_bin = fopen(bin_fpstr,"r");
  if (fp_bin != NULL) {

    /* process the binary file instead */

    fclose(fp_bin);
    parse_input_bin(bin_fpstr);

  } else {

    if (strcmp(format,MOLCAS_FORMAT) <= 0) {
      /* reduce the molcas output to a temp file */
      parse_molout(fn_relpath, tmp_fpstr, len_fn-1);

      /* check so that the tmp file exists, and then read it */
      fp_tmp = fopen(fn_relpath,"r");
      if (fp_tmp == NULL) {
        fprintf(stderr, "parse_input.c, function parse_input: unable to locate file %s for processing.\n",fn_relpath);
        printf( "program terminating due to the previous error.\n");
        fclose(fp_tmp);
        exit(1);
      } else {
        parse_input_tmp(state_er, tmp_fpstr, bin_fpstr);
      }
    }
    else {
      /* the tmp file path was provided in the input */
      parse_input_tmp(state_er, fn_relpath, bin_fpstr);
    }
  }

  if ((state_er[1] == state_er[5]) &&
      (state_er[2] == state_er[6])){
    printf( "adding elastic\n" );
    add_elastic(state_er);
    printf( "added elastic\n" );
  }

  if (check_pi() == 0) {
    fprintf(stderr, "\n\nparse_input.c, function parse_input: input matrix integrity check failure.\n");
    printf( "program terminating due to the previous error.\n");
    exit(EXIT_FAILURE);
  }

  set_tmax();

  return EXIT_SUCCESS;
}
