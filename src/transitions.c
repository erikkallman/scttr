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
#include <math.h>
#include "quicksort.h"
#include "transitions.h"
#include "std_num_ops.h"
#include "sci_const.h"
#include "sctr_input.h"

int
get_erange (sctr_input s_inp,
            double e){

  int j;

  for (j=1; j<s_inp->md->state_er[0]; j +=2 ) {
    if (inrange((e-s_inp->e0)*AUTOEV,s_inp->md->state_er[j],s_inp->md->state_er[j+1])) {
        return j-1;
    }
  }

  fprintf(stderr, "state of energy %le %le is not inside any of the energy ranges\
 provided in the input\n",e, s_inp->e0);
  printf( "program terminating due to the previous error.\n");
  exit(EXIT_FAILURE);
}

int
get_i (sctr_input s_inp,
       int from){

  int last_i = (int)s_inp->trs[0][0];
  int j = 0;
  while (last_i != -1) {

    if ((int)s_inp->trs[0][j] == from) {
      return j;
    }
    j++;
    last_i = (int)s_inp->trs[0][j];
  }

  return -1;
}

int
get_trs (int from,
        double ** trs){

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
get_il (sctr_input s_inp,
        int from){

  if (from>s_inp->n_states) {
    return -1;
  } else {
    return s_inp->idx_map[from-1];
  }
}

int
get_inext (sctr_input s_inp,
           int from){

  int j      = 0;

  while ((int)s_inp->trs[0][j] != -1) {
    if ((int)s_inp->trs[0][j] == from) {
      break;
    }
    j++;
  }

  if ((int)s_inp->trs[0][j] != from) {

    return (int)s_inp->trs[0][j];
  }

  while((int)s_inp->trs[0][j] == from){
    j++;
  }

  return j;

}

int
get_ilnext (sctr_input s_inp,
            int from) {
  int j = 0;

  if (from > s_inp->n_states) {
    return -1;
  }

  if ((j = get_il(s_inp,from)) == -1) {
    return -1;
  }

  if ((int)s_inp->trs[0][j] != from) {

    return (int)s_inp->trs[0][j];
  }

  while((int)s_inp->trs[0][j] != -1){
    if ((int)s_inp->trs[0][j] != from) {
      return j;
    }
    j++;
  }

  return -1;
}

int
get_trsnext (double ** trs,
           int from){

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
eval_trs (sctr_input s_inp){

  int j,k,l,m;
  int n_fs,n_is;

  int last_i,curr_i;

  int n_sfs; /* number of screened final states */

  int is_num = 0;
  int n_is_tmp = 0;
  int is_pos = 0;
  int is_idx = 0;
  int nt = 0; /* number of transitions fround in the provided energy range */
  int wh = 0; /* write head for the gs2is array */

  float int_thrsh = s_inp->md->state_t[1];
  float bw_thrsh = s_inp->md->state_t[2];

  spectrum spec = init_spec(s_inp, 10, 10);
  spectrum spec1 = init_spec(s_inp, 10, 10);

  double e_gs,e_is,e_fs;
  double tmom_gi,tmom_if;
  double bw;
  double int_tot,int_scr;
  double ** tmp_evals = malloc(sizeof(double*)*2);

  /* array to accumulate the distribution of transition intensities */
  double ** int_dist  = malloc(sizeof(double*)*2);

  /* keep check on what states that have been processed to make sure that no
     "to" state appears in the pi matrix more than once  */
  int * proc_st;
  int n_proc = 0;
  int is_proc;
  s_inp -> n_states = 0;

  if((proc_st = malloc(s_inp -> n_trans*sizeof(double))) == NULL ){
    fprintf(stderr, "parse_input.c, malloc: failed to allocate memory for\
 \"proc_st\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  if((s_inp -> idx_map = malloc((s_inp -> n_trans+1)*sizeof(int))) == NULL ){
    fprintf(stderr, "parse_input.c, malloc: failed to allocate memory for\
 \"s_inp -> idx_map\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  for (j=0; j<s_inp -> n_trans; j++) {
    s_inp -> idx_map[j] = -1;
  }

  s_inp -> idx_map[s_inp->n_trans] = -2; /* mark the end of the list */
  s_inp -> idx_map[0] = 0;

  printf( "\n      checking the integrity of the input matrix ..");
  fflush(stdout);
  /* check so that every state in PI can be reached with the get_i function */
  last_i = -2;
  curr_i = 0;
  j = 0;

  s_inp -> tmax_d = s_inp -> tmax_q = -1;

  while((int)s_inp->trs[0][j] != -1) {

    curr_i = (int)s_inp->trs[0][j];

    /* make sure the transition is not taking place between states
       in the same energy range interval */
    if (get_erange(s_inp,s_inp->trs[2][j]) != get_erange(s_inp,s_inp->trs[3][j])) {
      if ((int)s_inp->trs[5][j] == 1) {
        if ( s_inp->trs[4][j] > s_inp->tmax_d) {
          s_inp->tmax_d = s_inp->trs[4][j];
        }
      }
      else{
        if (s_inp->trs[4][j] > s_inp->tmax_q) {
          s_inp->tmax_q = s_inp->trs[4][j];
        }
      }
    }

    if (curr_i != last_i) {
      if (curr_i != (int)s_inp->trs[0][get_i(s_inp,curr_i)]) {
        return -1;
      }
      s_inp -> idx_map[curr_i-1] = j;
    }

    e_fs = s_inp->trs[3][j];
    e_is = s_inp->trs[2][j];

    /* is this a transition from an intermediate to final state? */
    if (inrange((e_fs - s_inp->e0)*AUTOEV,s_inp->md->state_er[5],s_inp->md->state_er[6]) &&
        inrange((e_is - s_inp->e0)*AUTOEV,s_inp->md->state_er[3],s_inp->md->state_er[4])) {

      da_append(spec->is2fs,j);
      is_num = s_inp->trs[0][j];

      /* have the gs-is transitions for this particular is already
         been processed? */
      is_proc = intinint(proc_st,is_num,n_proc);

      if (is_proc == -1) {

        /* loop over the list of transitions looking for corresponding
           g.state -> i-state transitions */
        k = 0;
        n_is_tmp = spec->gs2is->n_el;
        while((int)s_inp->trs[0][k] != -1) {
          e_gs = s_inp->trs[2][k];
          bw = get_rbdist(s_inp->e0,e_gs);
          if ((((int)s_inp->trs[1][k]) == is_num) &&
              inrange((e_gs - s_inp->e0)*AUTOEV,s_inp->md->state_er[1],\
                      s_inp->md->state_er[2]) &&
              (bw > bw_thrsh)
              ) {

            da_append(spec->is_idxs,spec->gs2is->n_el);
            da_append(spec->gs2is,k);
            nt++;
          }
          k++;
        }
        if (n_is_tmp == spec->gs2is->n_el) {
          /* no gs2is transitions found for that is2fs transition */
          spec->is2fs->n_el--;
        }
        else{
          /* gs2is transitions were added */

          da_append(spec->is_idxs,-1);
          da_append(spec->ii_start,wh);

          wh = spec->is_idxs->n_el;
          proc_st[n_proc++] = is_num;
        }
      }
      else {

        /* find the intermediate state in the gs2trs matrix */
        for (l=0; l<spec->is_idxs->n_el; l++) {
          is_idx = spec->is_idxs->a[l];
          if (is_idx != -1) {
            if ((int)s_inp->trs[1][spec->gs2is->a[is_idx]] == is_num) {
              break;
            }
          }
        }
        da_append(spec->ii_start,l);
        for (l=l; spec->is_idxs->a[l] != -1; l++) {
          nt++;
        }
      }
    }
    last_i = (int)s_inp->trs[0][j];
    j++;
  }
  spec -> n_st = nt;
  n_sfs = spec->is2fs->n_el;

  for (j=0; j<2; j++) {
    if((int_dist[j] = malloc(nt*sizeof(double))) == NULL ){
      fprintf(stderr, "transitions.c:function eval_trs, malloc: failed \
to allocate memory for \"int_dist[%d]\"\n",j);
      printf( "program terminating due to the previous error.\n");
      exit(1);
    }
  }

  l = 0;
  int_tot = 0;
  for (j=0; j<n_sfs; j++) {

    is_pos = spec->ii_start->a[j]; /* position at is_idx where gs2is
                                     idx will start */
    is_idx = spec->is_idxs->a[is_pos];
    tmom_if = s_inp->trs[4][spec->is2fs->a[j]];

    for (k=is_pos; ((is_idx = spec->is_idxs->a[k]) != -1); k++){

      tmom_gi = s_inp->trs[4][spec->gs2is->a[is_idx]];
      bw = get_rbdist(s_inp->e0,s_inp->trs[2][spec->gs2is->a[is_idx]]);
      int_dist[0][l] = l;

      int_dist[1][l] = tmom_if*tmom_gi*bw;
      int_tot += int_dist[1][l];
      l++;
    }
  }

  /* sort ascendingly according to energy */
  quicksort_d(int_dist[1],int_dist[0],0,nt-1,nt);

  /* screen until retaining x% of the intensity of all transitions,
   set all "sceened" states to the n_trans+1 */
  int_scr = 0;

  for (j=0; j<nt; j++) {
    if ((int_scr/int_tot) <= int_thrsh) {
      int_scr += int_dist[1][j];
    }
    else {
      j--;
      int_dist[0][j] = -1;
      break;
    }
  }

  quicksort_d(int_dist[0],int_dist[1],j,nt-1,nt);

  tmp_evals[0] = malloc(spec->is2fs->n_el*sizeof(double));
  tmp_evals[1] = malloc(spec->is2fs->n_el*sizeof(double));

  /* finally, generate the new screen by only including those states that were not screened out above */
  n_fs = n_is = 0;
  l = 0;
  m = j+1;
  wh = 0;

  for (j=0; j<n_sfs; j++) {

    is_pos = spec->ii_start->a[j]; /* position at is_idx where gs2is
                                     idx will start */
    is_idx = spec->is_idxs->a[is_pos];
    is_num = (int)s_inp->trs[0][spec->is2fs->a[j]];

    n_is_tmp = spec1->is_idxs->n_el;

    for (k=is_pos; ((is_idx = spec->is_idxs->a[k]) != -1); k++){
      if (l == (int)int_dist[0][m]) {

        /* is there an index in spec1 that already maps to
           this gs2is transition?*/
        if ((is_proc = intinint(spec1->gs2is->a,spec->gs2is->\
                                a[is_idx],spec1->gs2is->n_el)) != -1) {
          da_append(spec1->is_idxs,is_proc);
        }
        else {
          da_append(spec1->is_idxs,spec1->gs2is->n_el);
          da_append(spec1->gs2is,spec->gs2is->a[is_idx]);
        }
        spec1->n_st++;
        m++; /* start looking for the next state in int_dist */
        if (m == nt) {
          break;
        }
      }
      l++;
    }
    /* check if any gs2is transitions were added */
    if (n_is_tmp != spec1->is_idxs->n_el) {
      da_append(spec1->is_idxs,-1);
      da_append(spec1->is2fs,spec->is2fs->a[j]);
      da_append(spec1->ii_start,wh);

      wh = spec1->is_idxs->n_el;
      tmp_evals[0][n_is++] = (s_inp->trs[3][spec->gs2is->a[spec->is_idxs->a[is_pos]]] - s_inp->e0)*AUTOEV;
      tmp_evals[1][n_fs++] = (s_inp->trs[3][spec->is2fs->a[j]] - s_inp->e0)*AUTOEV;
    }
    if (m == nt) {
      break;
    }
  }

  spec1->emin_x = floor((get_minl(tmp_evals[0],n_is)-1))/AUTOEV;
  spec1->emax_x = ceil((get_maxl(tmp_evals[0],n_is) + 1))/AUTOEV;
  spec1->emin_y = floor((get_minl(tmp_evals[1],n_fs)-1))/AUTOEV;
  spec1->emax_y = ceil((get_maxl(tmp_evals[1],n_fs) + 1))/AUTOEV;

  s_inp -> n_states = n_proc;

  s_inp -> scr = spec1;
  free_spec(scr);

  if ((spec1->is2fs->n_el == 0) || (spec1->gs2is->n_el == 0)) {
    fprintf(stderr, "calc_spec.c, : no intermediate or final states states were\
 found in the energy range you provided. (spec1->is2fs->n_el = %d,\
 spec->gs2is->n_el = %d)\n",spec->is2fs->n_el,spec->gs2is->n_el);
    printf( "program terminating due to the previous error.\n");
    exit(EXIT_FAILURE);
  }

  printf( "          100%%\r");
  printf( "\n      done.\n");

  free(tmp_evals[0]);
  free(tmp_evals[1]);
  free(tmp_evals);
  free(proc_st);

  return 0;
}

void
count_states (sctr_input s_inp){

  int j = 0; /* looping variables */
  int last_i = -2;
  s_inp -> n_gfs = 0;
  s_inp -> n_is = 0;
  s_inp -> n_tmax = 0;

  int s_idx;
  int t_max;

  int n_proc = 0;
  int * t;
  int * proc_is;

  double * state_er = s_inp -> md -> state_er;

  if((proc_is = malloc(s_inp -> n_trans*sizeof(int))) == NULL ){
    fprintf(stderr, "parse_input.c, count_states: malloc: failed to allocate memory for\
 \"proc_is\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(EXIT_FAILURE);
  }

  if((t = malloc(s_inp -> n_trans * sizeof(int))) == NULL ){
    fprintf(stderr, "parse_input.c, count_states: malloc: failed to allocate memory for\
p \"t\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(EXIT_FAILURE);
  }

  for (j=0; j<s_inp->n_trans; j++) {
    t[j] = 0;
  }

  j = 0;

  /* gount all ground and final states */
  while((int)(s_inp->trs[0][j]) != -1){
    if(((s_inp->trs[2][j]-s_inp->e0)*AUTOEV >= state_er[1]) && ((s_inp->trs[2][j]-s_inp->e0)*AUTOEV <= state_er[2])){
      if (last_i != (int)s_inp->trs[0][j]) {
        s_inp -> n_gfs++;
      }
      if(((s_inp->trs[3][j]-s_inp->e0)*AUTOEV >= state_er[3]) && ((s_inp->trs[3][j]-s_inp->e0)*AUTOEV <= state_er[4])){
        if ((s_idx = intinint(proc_is,(int)s_inp->trs[1][j],n_proc)) == -1) {
          proc_is[n_proc] = (int)s_inp->trs[1][j];
          t[n_proc] += 1;
          n_proc++;
        }
        else{
          t[s_idx] += 1;
        }
        s_inp -> n_is++;
      }
    }

    last_i = (int)s_inp->trs[0][j++];
  }

  t_max = 0;
  for (j=0; j<n_proc; j++) {
    if (t[j] > t_max) {
      t_max = t[j];
    }
  }

  s_inp -> n_tmax = t_max;
  free(proc_is);
  free(t);
}

int
add_sym (sctr_input s_inp) {

    printf( "      adding elastic transitions ..\n");
  /* allocate memory for trs_el that is at most the size of ngs*nis*nfs*/
  /* copy the entire trs buffer to trs_el */
  /* make TRS point to trs_el instead */
  /* deallocate the memory to TRS */

  int j,k,l;
  int n_proc,nb;
  int tmp_idx,next_to;
  int last_i;
  int nt_el = s_inp->n_trans;
  int sz_buf = s_inp->n_tmax;
  /* in the worst case, there is an elastic transition from every intermediate state, to every final state. */
  long int sz_el = ((s_inp->n_trans*s_inp->n_gfs)*2)+1;

  if (nt_el > sz_el) {
    fprintf(stderr, "parse_input.c, function add_sym: input buffer writing outside its memory. nt_el = %d >= sz_el = %ld.\n",nt_el,sz_el);
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  double e0 = s_inp->e0;
  double * state_er = s_inp -> md -> state_er;

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

  if((proc_st = malloc(s_inp->n_trans*sizeof(int))) == NULL ){
    fprintf(stderr, "parse_input.c , function add_sym, malloc: failed to allocate memory for\
 \"proc_st\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  /* copy the old trans data to trs_el */
  for (l=0; l<s_inp->n_trans; /* nt_el++, */l++) {
    trs_el[0][l] = s_inp->trs[0][l];
    trs_el[1][l] = s_inp->trs[1][l];
    trs_el[2][l] = s_inp->trs[2][l];
    trs_el[3][l] = s_inp->trs[3][l];
    trs_el[4][l] = s_inp->trs[4][l];
    trs_el[5][l] = s_inp->trs[5][l];
    next_to = (int)s_inp->trs[0][l];
  }

  trs_el[0][l] = s_inp->trs[0][l] = -1;

  j = nb = n_proc = 0;
  while ((int)s_inp->trs[0][j] > 0) {
    printf( "        %.2f%%\r", (((float)j/(float)s_inp->n_trans)*100));
    if ((intinint(proc_st, (int)s_inp->trs[1][j], n_proc) == -1) &&
        (((s_inp->trs[2][j]-e0)*AUTOEV >= state_er[1]) && ((s_inp->trs[2][j]-e0)*AUTOEV <= state_er[2])) &&
        (((s_inp->trs[3][j]-e0)*AUTOEV >= state_er[3]) && ((s_inp->trs[3][j]-e0)*AUTOEV <= state_er[4]))
        ) {

      next_to = (int)s_inp->trs[1][j];

      /* found a "to" state that has not had its sym transitions
       added yet. loop over the trs matrix and check if there are transitions
       fromn other states that need to be taken into account. */

      nb = 0;

      l = j;
      while ((int)s_inp->trs[0][l] != -1) {
        if (((int)s_inp->trs[1][l] == next_to) &&
            (((s_inp->trs[2][l]-e0)*AUTOEV >= state_er[1]) && ((s_inp->trs[2][l]-e0)*AUTOEV <= state_er[2]))
            ){

          /* transitions from another state */
          trs_buf[0][nb] = s_inp->trs[1][l];
          trs_buf[1][nb] = s_inp->trs[0][l];
          trs_buf[2][nb] = s_inp->trs[3][l];
          trs_buf[3][nb] = s_inp->trs[2][l];
          trs_buf[4][nb] = s_inp->trs[4][l];
          trs_buf[5][nb] = s_inp->trs[5][l];

          nb++;

          /* jump to the next "from" state, since any given state can only have
           one transition to another specific state */
          last_i = s_inp->trs[0][l];
          while((int)s_inp->trs[0][l] != -1){
            if ((int)s_inp->trs[0][l] != last_i) {
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

  trs_el[0][nt_el] = -1;

  for (j=0; j<6; j++) {
    free(s_inp->trs[j]);
  }
  free(s_inp->trs);
  s_inp->trs = NULL;

  /* allocate new space for trs, now that we know the total size */
  if((s_inp->trs = malloc(6*sizeof(double*))) == NULL ){
    fprintf(stderr, "parse_input_molcas, malloc: failed to allocate memory for\
 \"input_data\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  for (j=0; j<6; j++) {
    if((s_inp->trs[j] = malloc((nt_el+1)*sizeof(double))) == NULL ){
      fprintf(stderr, "parse_input_molcas, malloc: failed to allocate memory \
for pointe rs in \"input_data\"\n");
      printf( "program terminating due to the previous error.\n");
      exit(1);
    }
  }
  for (l=0; l<6; l++) {
    for (j=0; j<=nt_el; j++) {
        s_inp->trs[l][j] = trs_el[l][j];
    }
  }

  s_inp->trs[0][nt_el] = -1;

  for (j=0; j<6; j++) {
    free(trs_buf[j]);
  }
  free(trs_buf);

  for (j=0; j<6; j++) {
    free(trs_el[j]);
  }
  free(trs_el);
  free(proc_st);

  s_inp->n_trans = nt_el;

  printf( "\n      done.\n");
  return 1;
}
