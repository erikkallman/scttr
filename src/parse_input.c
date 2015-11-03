/* This file is part of Scatter. */

/* Scatter is free software: you can redistribute it and/or modify */
/* it under the terms of the GNU Lesser General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or */
/* (at your option) any later version. */

/* Scatter is distributed in the hope that it will be useful, */
/* but without any warranty; without even the implied warranty of */
/* merchantability or fitness for a particular purpose. See the */
/* GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License */
/* along with Scatter, found in the "license" subdirectory of the root */
/* directory of the Scatter program. If not, see <http://www.gnu.org/licenses/>. */
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <signal.h>
#include <string.h>
#include <fcntl.h>
#include <errno.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "dyn_array.h"
#include "sctr_input.h"
#include "quicksort.h"
#include "std_char_ops.h"
#include "std_num_ops.h"
#include "get_numsl.h"
#include "parse_input.h"
#include "formats.h"
#include "sci_const.h"
#include "transitions.h"

<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD:src/io/parse_input.c
int
get_erange (spec_info s,
            double e)  {
  int j;

  for (j=1; j<s->md->state_er[0]; j +=2 ) {
    if (inrange((e-s->e0)*AUTOEV,s->md->state_er[j],s->md->state_er[j+1])) {
        return j-1;
    }
  }

  fprintf(stderr, "state of energy %le %le is not inside any of the energy ranges\
 provided in the input\n",e, s->e0);
  printf( "program terminating due to the previous error.\n");
  exit(EXIT_FAILURE);
}

/* function add_sym

   * synopsis:

   * algorithm:

   * input:

   * output:

   * side-effects:

   */
int
add_sym (spec_info s) {

    printf( "      adding elastic transitions ..\n");
  /* allocate memory for trs_el that is at most the size of ngs*nis*nfs*/
  /* copy the entire trs buffer to trs_el */
  /* make TRS point to trs_el instead */
  /* deallocate the memory to TRS */

  int j,k,l;
  int n_proc,nb;
  int tmp_idx,next_to;
  int last_i;
  int nt_el = s->n_trans;
  int sz_buf = s->n_tmax;
  /* in the worst case, there is an elastic transition from every intermediate state, to every final state. */
  long int sz_el = ((s->n_trans*s->n_gfs)*2)+1;

  if (nt_el > sz_el) {
    fprintf(stderr, "parse_input.c, function add_sym: input buffer writing outside its memory. nt_el = %d >= sz_el = %ld.\n",nt_el,sz_el);
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  double e0 = s->e0;
  double * state_er = s -> md -> state_er;

  /* which means that at most, we might have to read sz2 states into the buffer */

  /* all sym transitions handled so far */
  int * proc_st;

  double ** trs_buf;

  /* transition matrix with enough space to acommodate the sym transitions */
  double ** trs_el;

  if((trs_buf = malloc(6*sizeof(double*))) == NULL ){
    fprintf(stderr, "parse_input_molcas, malloc: failed to allocate memory for\
 \"trs_buf\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  for (j=0; j<6; j++) {
    if((trs_buf[j] = malloc(sz_buf*sizeof(double))) == NULL ){
      fprintf(stderr, "parse_input_molcas, malloc: failed to allocate memory \
for pointers in \"trs_buf[%d]\"\n",j);
      printf( "program terminating due to the previous error.\n");
      exit(1);
    }
  }

  if((trs_el = malloc(6*sizeof(double*))) == NULL ){
    fprintf(stderr, "parse_input_molcas, malloc: failed to allocate me mory for\
 \"trs_el\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  for (j=0; j<6; j++) {
    if((trs_el[j] = malloc((sz_el+1)*sizeof(double))) == NULL ){
      fprintf(stderr, "parse_input_molcas, malloc: failed to allocate memory \
for pointers in \"trs_el[%d]\"\n",j);
      printf( "program terminating due to the previous error.\n");
      exit(1);
    }
  }

  if((proc_st = malloc(s->n_trans*sizeof(int))) == NULL ){
    fprintf(stderr, "parse_input.c , function add_sym, malloc: failed to allocate memory for\
 \"proc_st\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

<<<<<<< HEAD
  /* copy the old pi data to pi_el */
<<<<<<< HEAD
  for (l=0; l<=nt; /* nt_el++, */l++) {
    if (nt_el > sz_el) {
      fprintf(stderr, "parse_input.c, function add_sym: input buffer writing outside its memory. nt_el = %d >= nt*sz = %ld.\n",nt_el,sz_el);
      printf( "program terminating due to the previous error.\n");
      exit(1);
    }
=======
  for (l=0; l<nt; /* nt_el++, */l++) {

>>>>>>> 17f48ce... first function reformulations on the track towards paralellism

    pi_el[0][l] = parsed_input[0][l];
    pi_el[1][l] = parsed_input[1][l];
    pi_el[2][l] = parsed_input[2][l];
    pi_el[3][l] = parsed_input[3][l];
    pi_el[4][l] = parsed_input[4][l];
    pi_el[5][l] = parsed_input[5][l];
    next_to = (int)parsed_input[0][l];
  }
  pi_el[0][l] = parsed_input[0][l] = -1;
=======
  /* copy the old trans data to trs_el */
  for (l=0; l<s->n_trans; /* nt_el++, */l++) {
    trs_el[0][l] = s->trs[0][l];
    trs_el[1][l] = s->trs[1][l];
    trs_el[2][l] = s->trs[2][l];
    trs_el[3][l] = s->trs[3][l];
    trs_el[4][l] = s->trs[4][l];
    trs_el[5][l] = s->trs[5][l];
    next_to = (int)s->trs[0][l];
  }

<<<<<<< HEAD
  trs_el[0][l] = trs[0][l] = -1;
>>>>>>> 046ed2a... implementing alternative KH formulation
=======
  trs_el[0][l] = s->trs[0][l] = -1;
>>>>>>> 40b7385... branch ready for testing

  j = nb = n_proc = 0;
  while ((int)s->trs[0][j] > 0) {
    printf( "        %.2f%%\r", (((float)j/(float)s->n_trans)*100));
    if ((intinint(proc_st, (int)s->trs[1][j], n_proc) == -1) &&
        (((s->trs[2][j]-e0)*AUTOEV >= state_er[1]) && ((s->trs[2][j]-e0)*AUTOEV <= state_er[2])) &&
        (((s->trs[3][j]-e0)*AUTOEV >= state_er[3]) && ((s->trs[3][j]-e0)*AUTOEV <= state_er[4]))
        ) {

      next_to = (int)s->trs[1][j];

      /* found a "to" state that has not had its sym transitions
       added yet. loop over the trs matrix and check if there are transitions
       fromn other states that need to be taken into account. */

      nb = 0;

      l = j;
      while ((int)s->trs[0][l] != -1) {
        if (((int)s->trs[1][l] == next_to) &&
            (((s->trs[2][l]-e0)*AUTOEV >= state_er[1]) && ((s->trs[2][l]-e0)*AUTOEV <= state_er[2]))
            ){

          /* transitions from another state */
          trs_buf[0][nb] = s->trs[1][l];
          trs_buf[1][nb] = s->trs[0][l];
          trs_buf[2][nb] = s->trs[3][l];
          trs_buf[3][nb] = s->trs[2][l];
          trs_buf[4][nb] = s->trs[4][l];
          trs_buf[5][nb] = s->trs[5][l];

          nb++;

          /* jump to the next "from" state, since any given state can only have
           one transition to another specific state */
          last_i = s->trs[0][l];
          while((int)s->trs[0][l] != -1){
            if ((int)s->trs[0][l] != last_i) {
              break;
            }
            l++;
          }
        }
        else {
          l++;
        }
      }

      /* append the data to the trs matrix */
      fflush(stdout);

      if ((nt_el+nb+1) > sz_el) {
      fprintf(stderr, "parse_input.c, function add_sym: input buffer writing outside its memory. nt_el+nb+1 = %d >= sz_el = %ld.\n",nt_el+nb+1,sz_el);
        printf( "program terminating due to the previous error.\n");
        exit(1);
      }

      /* if the from state cant be found in trs, just store the data in the last available place in trs_el */

      /* otherwise use the fwdsplice function to add it to trs_el */
      if (get_trs(next_to,trs_el) == -1) {
      /* if (get_i(next_to) == -1) { */

        for (k=0; k<nb; nt_el++,k++) {

          trs_el[0][nt_el] = trs_buf[0][k];
          trs_el[1][nt_el] = trs_buf[1][k];
          trs_el[2][nt_el] = trs_buf[2][k];
          trs_el[3][nt_el] = trs_buf[3][k];
          trs_el[4][nt_el] = trs_buf[4][k];
          trs_el[5][nt_el] = trs_buf[5][k];
        }
        trs_el[0][nt_el+1] = -1;

      }
      else {
        tmp_idx = get_trsnext(trs_el,next_to);
        fwdsplice(trs_buf,trs_el,tmp_idx,nt_el,nb,6);
        trs_el[0][nt_el+nb+1] = -1;
        nt_el+=nb;
      }
      proc_st[n_proc++] = next_to;
    }
    j++;
  }

  printf( "          100%%\r" );
  trs_el[0][nt_el] = -1;

  for (j=0; j<6; j++) {
    free(s->trs[j]);
  }
  free(s->trs);
  s->trs = NULL;

  /* allocate new space for trs, now that we know the total size */
  if((s->trs = malloc(6*sizeof(double*))) == NULL ){
    fprintf(stderr, "parse_input_molcas, malloc: failed to allocate memory for\
 \"input_data\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  for (j=0; j<6; j++) {
    if((s->trs[j] = malloc((nt_el+1)*sizeof(double))) == NULL ){
      fprintf(stderr, "parse_input_molcas, malloc: failed to allocate memory \
for pointe rs in \"input_data\"\n");
      printf( "program terminating due to the previous error.\n");
      exit(1);
    }
  }
  for (l=0; l<6; l++) {
    for (j=0; j<=nt_el; j++) {
        s->trs[l][j] = trs_el[l][j];
    }
  }

  s->trs[0][nt_el] = -1;

  for (j=0; j<6; j++) {
    free(trs_buf[j]);
  }
  free(trs_buf);

  for (j=0; j<6; j++) {
    free(trs_el[j]);
  }
  free(trs_el);
  free(proc_st);

  s->n_trans = nt_el;

  printf( "\n      done.\n");
  return 1;
}
=======
>>>>>>> f71c42e... new screening algorithm ready for cleaning and testing:src/parse_input.c

=======
>>>>>>> 8555cb4... splitting the eval_trs function
=======
>>>>>>> 8555cb4... splitting the eval_trs function
int
parse_input_bin (sctr_input s_inp,
                 char * bin_fpstr
                 ) {
  printf( "\n      processing binary file corresponding to the provided input file ..");
  fflush(stdout);
  int j; /* looping variables */

  int pi_xdim,pi_ydim;

  /* test writing to binary */
  FILE * fp_bin = fopen(bin_fpstr,"rb");

  if ( fread(&pi_ydim, sizeof(int),1,fp_bin) != 1) {
    fprintf(stderr, "parse_input.c, parse_input_bin: unable to read pi_ydim from binary output file \n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  pi_xdim = 6;
  s_inp -> n_trans = pi_ydim;

  if((s_inp -> trs = malloc(pi_xdim*sizeof(double*))) == NULL ){
    fprintf(stderr, "parse_input_bin, malloc: failed to allocate memory for\
 \"input_data\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  for (j=0; j<pi_xdim; j++) {
    if((s_inp -> trs[j] = malloc((pi_ydim+1)*sizeof(double))) == NULL ){
      fprintf(stderr, "parse_input_bin, malloc: failed to allocate memory \
for pointers in \"input_data\"\n");
      printf( "program terminating due to the previous error.\n");
      exit(1);
    }
  }

  for (j=0; j<6; j++) {
    if (fread(s_inp -> trs[j], sizeof(double), s_inp -> n_trans, fp_bin) != s_inp -> n_trans) {
      fprintf(stderr, "parse_input.c, function parse_input_bin: unable to read the PI matrix\n");
      printf( "program terminating due to the previous error.\n");
      fclose(fp_bin);
      exit(1);
    }
  }

  s_inp -> e0 = s_inp->trs[2][j];

  for (j=0; j<pi_ydim; j++) {
    if (s_inp -> trs[2][j] < s_inp->e0) {
      s_inp -> e0 = s_inp -> trs[2][j];
    }
  }

  s_inp -> trs[0][pi_ydim] = -1;

  if (fclose(fp_bin) != 0) {
    fprintf(stderr, "parse_input.c, function parse_input_bin: unable to close file:\n%s\n", bin_fpstr);
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }
  printf( " done\n");
  return 1;
}

int
parse_molout (sctr_input s_inp,
              char * fn_relpath,
              char * tmp_fpstr
              ) {
  printf( "\n      parsing the molcas .log file .. \n");
  int j,k,l,m; /* control loop variables */
  int mode; /* string matching mode flag */
  int string_flag = 0;
  int n_match = 0;

  FILE * fp_tmpdata;
  FILE * fp_relpath;

  int c; /* ntemporary char for storing input file characters */
  /* we are looking for four strings in the output file: one at the beginning
     of each data block, and one at the end. */
  char * s1 = NULL;
  char * s2 = NULL;
  char * s3 = NULL;

  char * str_buf = malloc(BUF_SIZE*5);
  if (s_inp->md->so_enrg == 0) {

    /* read spin-orbit data */
    s1 = malloc(35);
    s2 = malloc(40);
    s3 = malloc(44);

    s1 = "Eigenvalues of complex Hamiltonian";
    s2 = "Dipole transition strengths (SO states)";
    s3 = "Quadrupole transition strengths (SO states)";
  }
  else if(s_inp->md->so_enrg == 1){

    /* read spin-orbit free data */
    s1 = malloc(19);
    s2 = malloc(29);
    s3 = malloc(32);

    s1 = "SPIN-FREE ENERGIES";
    s2 = "Dipole transition strengths:";
    s3 = "Quadrupole transition strengths:";
  }

  const char * lookup_str[3] = {s1,s2,s3};

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
  for (j=0; ((c = fgetc(fp_relpath)) != EOF) && (n_match != 3); j++, k++) {
    if (j % 100000 == 0) {
      printf( "        %.2f%%\r", (((float)(j*sizeof(char))/(float)s_inp->md->sz_inp)*100));
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
          }

          for (m = 0; m<=k; m++) {
            fputc(str_buf[m],fp_tmpdata);
          }
        }
        /* if we find a flag while string_flag ==1 and mode ==1, we jhave
           read beyond the table */
        else if(string_flag == 1){
          string_flag = bin_flip(string_flag);

          /* something was read that was not a number flip back
             to searching for new line matches */
          mode    = bin_flip(mode);

          n_match++;

        }
      }
      /* check every line for a matching substring */
      /* mode = 1 and a line match means that we reached
         the end of this data block */
      else if((l = strscmp(str_buf,lookup_str,3)) != -1) {

        /* we found the first substring, the data we're looking for is
           inside the coming table of numbers of text. switch to mode 1.*/
        fprintf(fp_tmpdata,"%s\n",lookup_str[l]);
        mode    = bin_flip(mode);
      }
      k = 0;
    }
  }
  printf( "          100%%\r" );
  if ((fclose(fp_tmpdata) != 0) &&\
      (fclose(fp_relpath) != 0)) {
    fprintf(stderr, "parse_input.c, function parse_molout: unable to close files:\n%s\n%s\n", tmp_fpstr,fn_relpath);
    printf("program terminating due to the previous error.\n");
    exit(1);
  }
  free(str_buf);
  printf( "\n      done. \n");
  return 1;
}

int
parse_input_tmp (sctr_input s_inp,
                 char * fn_tmpdata,
                 char * bin_fpstr
                 ) {
  printf( "\n      parsing the tmp file .. \n");
  int j,k,l,m,j_test; /* control loop variables */
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
  int *     idxs_eigval;
  int * proc_idx;
  double * state_er = s_inp->md->state_er;
  double tmp_idx = 0;
  double maxr    = get_maxl(state_er,state_er[0]);
  double from_state_en, to_state_en, tmp_state_en;

  double *  e_eigval;
  double *  t_mom;

  double ** trs_buf;
  double ** trans_idxs;

  /* variables used later to approximate the amount of input data */
  int sz_tmp;
  struct stat st = {0};

  FILE * fp_tmpdata;

  int c; /* temporary char for storing input file characters */
  /* we are looking for four strings in the output file: one at the beginning
     of each data block, and one at the end. */
  char * str_buf = malloc(BUF_SIZE*5);
  char * s1 = NULL;
  char * s2 = NULL;
  char * s3 = NULL;

  if (s_inp->md->so_enrg == 0) {

    /* read spin-orbit data */

    s1 = "Eigenvalues of complex Hamiltonian";
    s2 = "Dipole transition strengths (SO states)";
    s3 = "Quadrupole transition strengths (SO states)";
  }
  else if(s_inp->md->so_enrg == 1){

    s1 = "SPIN-FREE ENERGIES";
    s2 = "Dipole transition strengths";
    s3 = "Quadrupole transition strengths";
  }

  /* create a pointer to the three data block beginners s1,s3, */
  const char * lookup_str[3] = {s1,s2,s3};

  stat(fn_tmpdata,&st);
  sz_tmp = (int)st.st_size; /* file size in bytes */
  sz_tmp = sz_tmp/34; /* shortest line = 34 chars */

  /* open the input file */
  if((fp_tmpdata = fopen(fn_tmpdata, "r")) == NULL) {
    fprintf(stderr,"parse_input: unable to open the input file %s.\n",fn_tmpdata);
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  k = 0; /* index for tmp_string */
  l = 0; /* index for lookup string */

  /* storage for the energy eigenvalues */
  if((e_eigval = malloc(sz_tmp*sizeof(double))) == NULL ){
    fprintf(stderr, "parse_input_molcas, malloc: failed to allocate memory for\
 \"e_eigval\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }
  if((trs_types = malloc(sz_tmp*sizeof(double))) == NULL ){
    fprintf(stderr, "parse_input_molcas, malloc: failed to allocate memory for\
 \"e_eigval\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }
  if((idxs_eigval = malloc(sz_tmp*sizeof(int))) == NULL ){
    fprintf(stderr, "parse_input_molcas, malloc: failed to allocate memory for\
 \"idxs_eigval\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  for (j=0; j<sz_tmp; j++) {
    idxs_eigval[j] = -1;
  }
  /* storage for the transition moments */
  if((t_mom = malloc(sz_tmp*sizeof(double))) == NULL ){
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
    if((trans_idxs[j] = malloc(sz_tmp*sizeof(double))) == NULL ){
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
  while ((c = fgetc(fp_tmpdata)) != EOF) {

    str_buf[l] = (char)c;

    if ((str_buf[l]  == '\n') && (l > 0)) { /* dont send blank lines */
      str_buf[l+1] = '\0';
      if ((j_test = strscmp(str_buf, lookup_str,3)) != -1) {

        j = j_test;
        trs_type = j;
      }
      else {

        if ((j == 0) && (isempty(str_buf,l) != 1)) {

          /* extract energy eigenvalues and state indexes */
          get_numsl(str_buf,num_idxs1,l,n_idxs1,&tmp_idx,&e_eigval[n_states]);
          if ((fabs(e_eigval[n_states]-e_eigval[0])*AUTOEV) \
              < maxr) {
            idxs_eigval[(int)tmp_idx-1] = n_states+1;
            n_states++;
          }
        }
        else if ((j > 0) && (isempty(str_buf,l) != 1)) {

          if (n_trans == 0) { /* pre-process the eigenvalues */

            /* sort the states in energy */
            quicksort(e_eigval,idxs_eigval,0,n_states-1,n_states);

            s_inp -> e0 = e_eigval[0];
            /* adjust n_states so that it accounts for states not to be read
             due to being outside of the input range */
            m = 0;

            for (k=0; k<n_states; k++) {
              if ((e_eigval[k]-s_inp->e0)*AUTOEV < maxr) {
                m++;
              }
            }
            n_states = m;
          }
          /* extract transition moments and transition indexes */
          get_numsl(str_buf,num_idxs2,l,n_idxs2,&trans_idxs[0][n_trans],\
                    &trans_idxs[1][n_trans],&t_mom[n_trans]);
          trs_types[n_trans] = trs_type;

          from_state_en = get_wi(e_eigval,idxs_eigval,(int)(trans_idxs[1][n_trans]),n_states);

          to_state_en = get_wi(e_eigval,idxs_eigval,(int)(trans_idxs[0][n_trans]),n_states);

          if ((((int)to_state_en != -1) && ((int)from_state_en != -1))
              && (((to_state_en-s_inp->e0)*AUTOEV < maxr) &&   \
              ((from_state_en-s_inp->e0)*AUTOEV < maxr)))
            {
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

  /* allocate space for the "parsed input matrix" that will be filled with data
     in the remaining sections of this function */

  if((proc_idx = malloc(n_states*sizeof(int))) == NULL ){
    fprintf(stderr, "parse_input_molcas, malloc: failed to allocate memory for\
 \"proc_idx\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  if((s_inp->trs = malloc(6*sizeof(double*))) == NULL ){
    fprintf(stderr, "parse_input_molcas, malloc: failed to allocate memory for\
 \"input_data\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  for (j=0; j<6; j++) {
    if((s_inp->trs[j] = malloc((n_trans+1)*sizeof(double))) == NULL ){
      fprintf(stderr, "parse_input_molcas, malloc: failed to allocate memory \
for pointers in \"input_data\"\n");
      printf( "program terminating due to the previous error.\n");
      exit(1);
    }
  }
  s_inp->trs[0][n_trans] = -1; /* end of list */

  if((trs_buf = malloc(6*sizeof(double*))) == NULL ){
    fprintf(stderr, "parse_input_molcas, malloc: failed to allocate memory for\
 \"input_data\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  for (j=0; j<6; j++) {
    if((trs_buf[j] = malloc((n_trans+1)*sizeof(double))) == NULL ){
      fprintf(stderr, "parse_input_molcas, malloc: failed to allocate memory \
for pointers in \"input_data\"\n");
      printf( "program terminating due to the previous error.\n");
      exit(1);
    }
  }

  /* finally, store the data in the trs matrix for the parse_input function  */
  n_proc      = 0;
  last_i = 1;
  j = k = l = 0;

  /* l = index of line in the data arrays */
  /* j = index of data read into trs_buf, reset upon reading the transitions
     from a new state */

  while (l<n_trans-1) {
    if (l % 1 == 0) {
      printf( "        %.2f%%\r", (((float)l/(float)n_trans)*100));
    }
    if ((j+l) == (n_trans)) {
      /* dont read beyond the last value in PI */
      trs_buf[0][j] = -1;
      idx_from = -1;

    }
    else {
      from_state_en = get_wi(e_eigval,idxs_eigval,(int)(trans_idxs[0][j+l]),n_states);
      to_state_en = get_wi(e_eigval,idxs_eigval,(int)(trans_idxs[1][j+l]),n_states);

      if (inrange((to_state_en-s_inp->e0)*AUTOEV,state_er[3],state_er[4]) &&
           (inrange((from_state_en-s_inp->e0)*AUTOEV,state_er[5],state_er[6])) &&
          (((int)state_er[1] != (int)state_er[5]) && ((int)state_er[2] != (int)state_er[6]))
          ) {

        idx_from = trans_idxs[1][j+l];
        idx_to = trans_idxs[0][j+l];

        tmp_state_en = from_state_en;
        from_state_en = to_state_en;
        to_state_en = tmp_state_en;
      }
      else {

        idx_from = trans_idxs[0][j+l];
        idx_to   = trans_idxs[1][j+l];
      }

      trs_buf[0][j] = idx_from;
      trs_buf[1][j] = idx_to;
      trs_buf[2][j] = from_state_en;
      trs_buf[3][j] = to_state_en;
      trs_buf[4][j] = t_mom[j+l];
      trs_buf[5][j] = trs_types[j+l];

    }
    if (idx_from != last_i) {

      tmp_idx2 = get_inext(s_inp,last_i);
      /* we have read all transitions for a state */
      /* check if the last_i has already been processed */
      if (intinint(proc_idx,last_i,n_proc) == -1) {
        m = l;
        while ((int)trs_buf[0][m-l] != idx_from) {

          s_inp->trs[0][m] = trs_buf[0][m-l];
          s_inp->trs[1][m] = trs_buf[1][m-l];
          s_inp->trs[2][m] = trs_buf[2][m-l];
          s_inp->trs[3][m] = trs_buf[3][m-l];
          s_inp->trs[4][m] = trs_buf[4][m-l];
          s_inp->trs[5][m] = trs_buf[5][m-l];

          m++;
        }

        s_inp->trs[0][m] = -1;
        proc_idx[n_proc++] = last_i;
      }
      else{
        tmp_idx2 = get_inext(s_inp,last_i);
        fwdsplice(trs_buf,s_inp->trs,tmp_idx2,l,j,6);
        fflush(stdout);
      }
      l += j;
      j = 0;

    } else {
      j++;
    }
    last_i = idx_from;
  }

  s_inp->trs[0][l] = -1;
  s_inp -> n_trans = l;
  printf( "          100%%\r");

  if (fclose(fp_tmpdata) != 0) {
    fprintf(stderr, "parse_input.c, function parse_input_tmp: unable to close file:\n%s\n", fn_tmpdata);
    printf("program terminating due to the previous error.\n");
    exit(1);
  }

  free(e_eigval);
  free(trs_types);
  free(idxs_eigval);
  free(t_mom);

  for (j=0; j<2; j++) {
    free(trans_idxs[j]);
  }
  free(trans_idxs);

  free(num_idxs1);
  free(num_idxs2);
  free(proc_idx);

  for (j = 0; j<6; j++) {
    free(trs_buf[j]);
  }
  free(trs_buf);


  free(str_buf);

  printf( "\n      done.\n");
  return EXIT_SUCCESS;

}

int
parse_input (sctr_input s_inp
             ){
  int j;
  int rc;                       /* return code */

  FILE * fp_tmp;
  FILE * fp_bin;

  metadata md = s_inp -> md;

  char * bin_fpstr = concs(3,md->outpath,md->inp_fn,bin_sfx);
  char * tmp_fpstr = concs(3,md->outpath,md->inp_fn,tmp_sfx);
  char * format = md->inp_sfx;
  char * inpath = md -> inpath;

  fp_bin = fopen(bin_fpstr,"r");
  if (fp_bin != NULL) {

    /* process the binary file instead */
    fclose(fp_bin);
    parse_input_bin(s_inp,bin_fpstr);
  }
  else {

    if (strcmp(format,MOLCAS_FORMAT) <= 0) {

      /* reduce the molcas output to a temp file */
      parse_molout(s_inp,inpath, tmp_fpstr );

      fp_tmp = fopen(inpath,"r");
      if (fp_tmp == NULL) {
        fprintf(stderr, "parse_input.c, function parse_input: unable to locate\
 file %s for processing.\n",inpath);
        printf( "program terminating due to the previous error.\n");
        fclose(fp_tmp);
        exit(1);
      } else {
        parse_input_tmp(s_inp, tmp_fpstr, bin_fpstr);
      }
    }
    else {
      /* the tmp file path was provided in the input */
      parse_input_tmp(s_inp, inpath, bin_fpstr);
    }
    if ((md->state_er[1] == md->state_er[5]) &&
        (md->state_er[2] == md->state_er[6])){

      count_states(s_inp);
      add_sym(s_inp);
    }

    /* the only way the parse_input_tmp function can get called is if there is if
       no binary file present. we can therefore safely write to the binary file
       without checking if it already exists. */
    if ((fp_bin = fopen(bin_fpstr,"wb")) == NULL) {
      fprintf(stderr, "parse_input.c, function parse_input_bin: unable to open\
the binary file used to store the PI matrix: %s\n", bin_fpstr);
      printf( "program terminating due to the previous error.\n");
      exit(1);
    }

    fwrite((const void*)&(s_inp->n_trans), sizeof(int), 1, fp_bin);
    fflush(fp_bin);
    fclose(fp_bin);

    if ((fp_bin = fopen(bin_fpstr,"ab")) == NULL) {
      fprintf(stderr, "parse_input.c, function parse_input_bin: unable to open\
the binary file used to store the PI matrix: %s\n", bin_fpstr);
      printf( "program terminating due to the previous error.\n");
      exit(1);
    }

    for (j=0; j<6; j++) {
      if (fwrite(s_inp->trs[j], sizeof(double), s_inp->n_trans, fp_bin) != s_inp->n_trans) {
        fprintf(stderr, "parse_input.c, function parse_input_bin: unable to \
write the PI matrix\n");
        printf( "program terminating due to the previous error.\n");
        fclose(fp_bin);
        exit(1);
      }
      fflush(fp_bin);
    }

    if (fclose(fp_bin) != 0){
      fprintf(stderr, "parse_input.c, function parse_input_tmp: unable to close\
files:\n%s\n", bin_fpstr);
      printf("program terminating due to the previous error.\n");
      exit(1);
    }
  }

  if ((rc = eval_trs(s_inp)) != 0) {

    fprintf(stderr, "\n\nparse_input.c, function eval_trs: input matrix \
integrity check failure.\n");

    if (rc == -1) {
      fprintf(stderr, "\nthe get_i function was unable to obtain the element index for some states.\n");
    }
    else if(rc >= 0){
      fprintf(stderr, "\nstate %d occured multiple times in the input matrix. \
Symmetric transitions were erronously added, most likely a failure in the fwdsplice function.\n",rc);
    }
    printf( "program terminating due to the previous error.\n");
    exit(EXIT_FAILURE);
  }

  free(bin_fpstr);
  free(tmp_fpstr);

  return EXIT_SUCCESS;
}
