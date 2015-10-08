#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <signal.h>
#include <string.h>
#include <fcntl.h>
#include <errno.h>
#include <unistd.h>
#include "structs.h"
#include "spec_info.h"
#include "quicksort.h"
#include "std_char_ops.h"
#include "std_num_ops.h"
#include "get_numsl.h"
#include "parse_input.h"
#include "formats.h"
#include "sci_const.h"

/* function find_range

   * synopsis:

   * algorithm:

   * input:

   * output:

   * side-effects:

   */
int
get_erange (spec_info s,
            double e)  {
  int j;

  for (j=1; j<s->md->state_er[0]; j +=2 ) {
    if (inrange((e-s->e0)*AUTOEV,s->md->state_er[j],s->md->state_er[j+1])) {
        return j-1;
    }
  }

  fprintf(stderr, "state of energy %le is not inside any of the energy ranges\
 provided in the input\n",e);
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
  double ** trs = s -> trs;
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

  /* free up the memory space occupied by the old matrix */
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

    trs_el[0][l] = trs[0][l];
    trs_el[1][l] = trs[1][l];
    trs_el[2][l] = trs[2][l];
    trs_el[3][l] = trs[3][l];
    trs_el[4][l] = trs[4][l];
    trs_el[5][l] = trs[5][l];
    next_to = (int)trs[0][l];
  }

  trs_el[0][l] = trs[0][l] = -1;
>>>>>>> 046ed2a... implementing alternative KH formulation

  j = nb = n_proc = 0;
  while ((int)trs[0][j] > 0) {
    printf( "        %.2f%%\r", (((float)j/(float)s->n_trans)*100));
    if ((intinint(proc_st, (int)trs[1][j], n_proc) == -1) &&
        (((trs[2][j]-e0)*AUTOEV >= state_er[1]) && ((trs[2][j]-e0)*AUTOEV <= state_er[2])) &&
        (((trs[3][j]-e0)*AUTOEV >= state_er[3]) && ((trs[3][j]-e0)*AUTOEV <= state_er[4]))
        ) {

      next_to = (int)trs[1][j];

      /* found a "to" state that has not had its sym transitions
       added yet. loop over the trs matrix and check if there are transitions
       fromn other states that need to be taken into account. */

      nb = 0;

      l = j;
      while ((int)trs[0][l] != -1) {
        if (((int)trs[1][l] == next_to) &&
            (((trs[2][l]-e0)*AUTOEV >= state_er[1]) && ((trs[2][l]-e0)*AUTOEV <= state_er[2]))
            ){

          /* transitions from another state */
          trs_buf[0][nb] = trs[1][l];
          trs_buf[1][nb] = trs[0][l];
          trs_buf[2][nb] = trs[3][l];
          trs_buf[3][nb] = trs[2][l];
          trs_buf[4][nb] = trs[4][l];
          trs_buf[5][nb] = trs[5][l];

          nb++;

          /* jump to the next "from" state, since any given state can only have
           one transition to another specific state */
          last_i = trs[0][l];
          while((int)trs[0][l] != -1){
            if ((int)trs[0][l] != last_i) {
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
    free(trs[j]);
  }
  free(trs);
  trs = NULL;

  /* allocate new space for trs, now that we know the total size */
  if((trs = malloc(6*sizeof(double*))) == NULL ){
    fprintf(stderr, "parse_input_molcas, malloc: failed to allocate memory for\
 \"input_data\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  for (j=0; j<6; j++) {
    if((trs[j] = malloc((nt_el+1)*sizeof(double))) == NULL ){
      fprintf(stderr, "parse_input_molcas, malloc: failed to allocate memory \
for pointe rs in \"input_data\"\n");
      printf( "program terminating due to the previous error.\n");
      exit(1);
    }
  }
  for (l=0; l<6; l++) {
    for (j=0; j<=nt_el; j++) {
        trs[l][j] = trs_el[l][j];
    }
  }

  trs[0][nt_el] = -1;

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

int
parse_input_bin (spec_info s,
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
  s -> n_trans = pi_ydim;

  if((s -> trs = malloc(pi_xdim*sizeof(double*))) == NULL ){
    fprintf(stderr, "parse_input_bin, malloc: failed to allocate memory for\
 \"input_data\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  for (j=0; j<pi_xdim; j++) {
    if((s -> trs[j] = malloc((pi_ydim+1)*sizeof(double))) == NULL ){
      fprintf(stderr, "parse_input_bin, malloc: failed to allocate memory \
for pointers in \"input_data\"\n");
      printf( "program terminating due to the previous error.\n");
      exit(1);
    }
  }

  for (j=0; j<6; j++) {
    if (fread(s -> trs[j], sizeof(double), s -> n_trans, fp_bin) != s -> n_trans) {
      fprintf(stderr, "parse_input.c, function parse_input_bin: unable to read the PI matrix\n");
      printf( "program terminating due to the previous error.\n");
      fclose(fp_bin);
      exit(1);
    }
  }

  s -> e0 = -1;

  for (j=0; j<pi_ydim; j++) {
    if (s -> trs[2][j] > s->e0) {
      s -> e0 = s -> trs[2][j];
    }
  }

  s -> trs[0][pi_ydim] = -1;

  if (fclose(fp_bin) != 0) {
    fprintf(stderr, "parse_input.c, function parse_input_bin: unable to close file:\n%s\n", bin_fpstr);
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }
  printf( " done\n");
  return 1;
}

int
parse_molout (spec_info s,
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
  if (s->md->etype == 0) {

    /* read spin-orbit data */
    s1 = malloc(35);
    s2 = malloc(40);
    s3 = malloc(44);

    s1 = "Eigenvalues of complex Hamiltonian";
    s2 = "Dipole transition strengths (SO states)";
    s3 = "Quadrupole transition strengths (SO states)";
  }
  else if(s->md->etype == 1){

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
      printf( "        %.2f%%\r", (((float)(j*sizeof(char))/(float)s->md->sz_inp)*100));
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

  printf( "\n      done. \n");
  return 1;
}

int
check_trs (spec_info s){

  int j;
  int last_i,curr_i;
  int prog_step;

  /* keep check on what states that have been processed to make sure that no
     "to" state appears in the pi matrix more than once  */
  int * proc_st;
  int n_proc = 0;
  s -> n_states = 0;

  if((proc_st = malloc(s -> n_trans*sizeof(double))) == NULL ){
    fprintf(stderr, "parse_input.c, malloc: failed to allocate memory for\
 \"proc_st\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  if((s -> idx_map = malloc((s -> n_trans+1)*sizeof(int))) == NULL ){
    fprintf(stderr, "parse_input.c, malloc: failed to allocate memory for\
 \"s -> idx_map\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  for (j=0; j<s -> n_trans; j++) {
    s -> idx_map[j] = -1;
  }

  s -> idx_map[s->n_trans] = -2; /* mark the end of the list */
  s -> idx_map[0] = 0;

  prog_step = s->n_trans/10;

  printf( "\n      checking the integrity of the input matrix ..\n");

  /* check so that every state in PI can be reached with the get_i function */
  last_i = -2;
  curr_i = 0;
  j = 0;

  s -> tmax_d = s -> tmax_q = -1;

  while((int)s->trs[0][j] != -1) {

    curr_i = (int)s->trs[0][j];

    /* make sure the transition is not taking place between states
     in the same energy range interval */
    if (get_erange(s,s->trs[2][j]) != get_erange(s,s->trs[3][j])) {
      /* printf( "other %d %d %d %d\n", (int)s->trs[0][j],(int)s->trs[1][j],get_erange(s,s->trs[2][j]),get_erange(s,s->trs[3][j])); */
      /* sleep(1); */
      if ((int)s->trs[5][j] == 1) {
        if ( s->trs[4][j] > s->tmax_d) {
          s->tmax_d = s->trs[4][j];
        }
      }
      else{
        if (s->trs[4][j] > s->tmax_q) {
          s->tmax_q = s->trs[4][j];
        }
      }
    } /* else { */
    /*   printf( "same %d %d %d %d\n", (int)s->trs[0][j],(int)s->trs[1][j],get_erange(s,s->trs[2][j]),get_erange(s,s->trs[3][j])); */
    /*   sleep(1); */
    /* } */


    if (j % prog_step == 0) {
      printf( "        %.2f%%\r", (((float)j/(float)s->n_trans)*100));
    }

    if (curr_i != last_i) {
      if (curr_i != (int)s->trs[0][get_i(s,curr_i)]) {
        return -1;
      }
      else if(intinint(proc_st, last_i, n_proc) != -1){
        return last_i;
      }

      s -> idx_map[curr_i-1] = j;
      proc_st[n_proc++] = last_i;
    }

    last_i = (int)s->trs[0][j];
    j++;
  }

  s -> n_states = n_proc;
  printf( "          100%%\r");
  printf( "\n      done.\n");

  free(proc_st);

  return 0;
}

int
get_i (spec_info s,
       int from
       ) {
  int last_i = (int)s->trs[0][0];
  int j = 0;
  while (last_i != -1) {

    if ((int)s->trs[0][j] == from) {
      return j;
    }
    j++;
    last_i = (int)s->trs[0][j];
  }

  return -1;
}

int
get_trs (int from,
        double ** trs
       ) {
  int last_i = (int)trs[0][0];
  int j = 0;
  while (last_i != -1) {

    if ((int)trs[0][j] == from) {
      return j;
    }
    j++;
    last_i = (int)trs[0][j];
  }

  return -1;
}


int
get_il (spec_info s,
        int from
       ) {

  if (from>s->n_states) {
    return -1;
  } else {
    return s->idx_map[from-1];
  }
}

int
get_inext (spec_info s,
           int from
           ) {

  int j      = 0;

  while ((int)s->trs[0][j] != -1) {
    if ((int)s->trs[0][j] == from) {
      break;
    }
    j++;
  }

  if ((int)s->trs[0][j] != from) {

    return (int)s->trs[0][j];
  }

  while((int)s->trs[0][j] == from){
    j++;
  }

  return j;

}

int
get_ilnext (spec_info s,
            int from
           ) {
  int j = 0;

  if (from > s->n_states) {
    return -1;
  }

  if ((j = get_il(s,from)) == -1) {
    return -1;
  }

  if ((int)s->trs[0][j] != from) {

    return (int)s->trs[0][j];
  }

  while((int)s->trs[0][j] != -1){
    if ((int)s->trs[0][j] != from) {
      return j;
    }
    j++;
  }

  return -1;
}

int
get_trsnext (double ** trs,
           int from
           ) {

  int j      = 0;

  while ((int)trs[0][j] != -1) {
    if ((int)trs[0][j] == from) {
      break;
    }
    j++;
  }

  if ((int)trs[0][j] != from) {

    return (int)trs[0][j];
  }

  while((int)trs[0][j] == from){
    j++;
  }

  return j;
}

int
parse_input_tmp (spec_info s,
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
  double * state_er = s->md->state_er;
  double tmp_idx = 0;
  double maxr    = get_maxl(state_er,state_er[0]);
  double from_state_en, to_state_en, tmp_state_en;

  double *  e_eigval;
  double *  t_mom;

  double ** trs_buf;
  double ** trans_idxs;

  FILE * fp_tmpdata;

  int c; /* temporary char for storing input file characters */
  /* we are looking for four strings in the output file: one at the beginning
     of each data block, and one at the end. */
  char * str_buf = malloc(BUF_SIZE*5);
  char * s1 = NULL;
  char * s2 = NULL;
  char * s3 = NULL;

  if (s->md->etype == 0) {

    /* read spin-orbit data */
    s1 = malloc(35);
    s2 = malloc(40);
    s3 = malloc(44);

    s1 = "Eigenvalues of complex Hamiltonian";
    s2 = "Dipole transition strengths (SO states)";
    s3 = "Quadrupole transition strengths (SO states)";
  }
  else if(s->md->etype == 1){
    s1 = malloc(19);
    s2 = malloc(28);
    s3 = malloc(32);

    s1 = "SPIN-FREE ENERGIES";
    s2 = "Dipole transition strengths";
    s3 = "Quadrupole transition strengths";
  }

  /* const char s1[38] = "Eigenvalues of complex Hamiltonian"; */
  /* const char s2[40] = "Dipole transition strengths (SO states)"; */
  /* const char s3[44] = "Quadrupole transition strengths (SO states)"; */

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
            s -> e0 = e_eigval[0];
            /* adjust n_states so that it accounts for states not to be read
             due to being outside of the input range */
            m = 0;

            for (k=0; k<n_states; k++) {
              if ((e_eigval[k]-s->e0)*AUTOEV < maxr) {
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
              && (((to_state_en-s->e0)*AUTOEV < maxr) &&   \
              ((from_state_en-s->e0)*AUTOEV < maxr)))
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

  if((s->trs = malloc(6*sizeof(double*))) == NULL ){
    fprintf(stderr, "parse_input_molcas, malloc: failed to allocate memory for\
 \"input_data\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  for (j=0; j<6; j++) {
    if((s->trs[j] = malloc((n_trans+1)*sizeof(double))) == NULL ){
      fprintf(stderr, "parse_input_molcas, malloc: failed to allocate memory \
for pointers in \"input_data\"\n");
      printf( "program terminating due to the previous error.\n");
      exit(1);
    }
  }
  s->trs[0][n_trans] = -1; /* end of list */

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

      if (inrange((to_state_en-s->e0)*AUTOEV,state_er[3],state_er[4]) &&
           (inrange((from_state_en-s->e0)*AUTOEV,state_er[5],state_er[6])) &&
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

      tmp_idx2 = get_inext(s,last_i);
      /* we have read all transitions for a state */
      /* check if the last_i has already been processed */
      if (intinint(proc_idx,last_i,n_proc) == -1) {
        m = l;
        while ((int)trs_buf[0][m-l] != idx_from) {

          s->trs[0][m] = trs_buf[0][m-l];
          s->trs[1][m] = trs_buf[1][m-l];
          s->trs[2][m] = trs_buf[2][m-l];
          s->trs[3][m] = trs_buf[3][m-l];
          s->trs[4][m] = trs_buf[4][m-l];
          s->trs[5][m] = trs_buf[5][m-l];

          m++;
        }

        s->trs[0][m] = -1;
        proc_idx[n_proc++] = last_i;
      }
      else{
        tmp_idx2 = get_inext(s,last_i);
        fwdsplice(trs_buf,s->trs,tmp_idx2,l,j,6);
        fflush(stdout);
      }
      l += j;
      j = 0;

    } else {
      j++;
    }
    last_i = idx_from;
  }

  s->trs[0][l] = -1;
  s -> n_trans = l;
  printf( "          100%%\r");

  if (fclose(fp_tmpdata) != 0) {
    fprintf(stderr, "parse_input.c, function parse_input_tmp: unable to close file:\n%s\n", fn_tmpdata);
    printf("program terminating due to the previous error.\n");
    exit(1);
  }

  free(num_idxs1);
  free(num_idxs2);
  free(trs_types);
  free(proc_idx);

  free(e_eigval);
  free(idxs_eigval);

  for (j = 0; j<6; j++) {
    free(trs_buf[j]);
  }
  free(trs_buf);

  for (j=0; j<2; j++) {
    free(trans_idxs[j]);
  }
  free(trans_idxs);
  free(t_mom);

  free(str_buf);

  printf( "\n      done.\n");
  return EXIT_SUCCESS;

}

void
count_states (spec_info s) {

  int j = 0; /* looping variables */
  int last_i = -2;
  s -> n_gfs = 0;
  s -> n_is = 0;
  s -> n_tmax = 0;

  int s_idx;
  int t_max;

  int n_proc = 0;
  int * t;
  int * proc_is;

  double * state_er = s -> md -> state_er;

  if((proc_is = malloc(s -> n_trans*sizeof(int))) == NULL ){
    fprintf(stderr, "parse_input.c, count_states: malloc: failed to allocate memory for\
 \"proc_is\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(EXIT_FAILURE);
  }

  if((t = malloc(s -> n_trans * sizeof(int))) == NULL ){
    fprintf(stderr, "parse_input.c, count_states: malloc: failed to allocate memory for\
p \"t\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(EXIT_FAILURE);
  }

  for (j=0; j<s->n_trans; j++) {
    t[j] = 0;
  }

  j = 0;

  /* gount all ground and final states */
  while((int)(s->trs[0][j]) != -1){
    if(((s->trs[2][j]-s->e0)*AUTOEV >= state_er[1]) && ((s->trs[2][j]-s->e0)*AUTOEV <= state_er[2])){
      if (last_i != (int)s->trs[0][j]) {
        s -> n_gfs++;
      }
      if(((s->trs[3][j]-s->e0)*AUTOEV >= state_er[3]) && ((s->trs[3][j]-s->e0)*AUTOEV <= state_er[4])){
        if ((s_idx = intinint(proc_is,(int)s->trs[1][j],n_proc)) == -1) {
          proc_is[n_proc] = (int)s->trs[1][j];
          t[n_proc] += 1;
          n_proc++;
        }
        else{
          t[s_idx] += 1;
        }
        s -> n_is++;
      }
    }

    last_i = (int)s->trs[0][j++];
  }

  t_max = 0;
  for (j=0; j<n_proc; j++) {
    if (t[j] > t_max) {
      t_max = t[j];
    }
  }

  s -> n_tmax = t_max;
  free(proc_is);
  free(t);
}

int
parse_input (spec_info s
             ){
  int j;
  int rc;                       /* return code */

  FILE * fp_tmp;
  FILE * fp_bin;

  metadata md = s -> md;

  char * bin_fpstr = concs(3,md->outpath,md->inp_fn,bin_sfx);
  char * tmp_fpstr = concs(3,md->outpath,md->inp_fn,tmp_sfx);
  char * format = md->inp_sfx;
  char * inpath = md -> inpath;

  fp_bin = fopen(bin_fpstr,"r");
  if (fp_bin != NULL) {

    /* process the binary file instead */
    fclose(fp_bin);
    parse_input_bin(s,bin_fpstr);
  }
  else {

    if (strcmp(format,MOLCAS_FORMAT) <= 0) {

      /* reduce the molcas output to a temp file */
      parse_molout(s,inpath, tmp_fpstr );

      fp_tmp = fopen(inpath,"r");
      if (fp_tmp == NULL) {
        fprintf(stderr, "parse_input.c, function parse_input: unable to locate\
 file %s for processing.\n",inpath);
        printf( "program terminating due to the previous error.\n");
        fclose(fp_tmp);
        exit(1);
      } else {
        parse_input_tmp(s, tmp_fpstr, bin_fpstr);
      }
    }
    else {
      /* the tmp file path was provided in the input */
      parse_input_tmp(s, inpath, bin_fpstr);
    }
    if ((md->state_er[1] == md->state_er[5]) &&
        (md->state_er[2] == md->state_er[6])){

      count_states(s);
      add_sym(s);
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

    fwrite((const void*)&(s->n_trans), sizeof(int), 1, fp_bin);
    fflush(fp_bin);
    fclose(fp_bin);

    if ((fp_bin = fopen(bin_fpstr,"ab")) == NULL) {
      fprintf(stderr, "parse_input.c, function parse_input_bin: unable to open\
the binary file used to store the PI matrix: %s\n", bin_fpstr);
      printf( "program terminating due to the previous error.\n");
      exit(1);
    }

    for (j=0; j<6; j++) {
      if (fwrite(s->trs[j], sizeof(double), s->n_trans, fp_bin) != s->n_trans) {
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

  if ((rc = check_trs(s)) != 0) {

    fprintf(stderr, "\n\nparse_input.c, function check_trs: input matrix \
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
