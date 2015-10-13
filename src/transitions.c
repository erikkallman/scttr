#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "quicksort.h"
#include "transitions.h"
#include "std_num_ops.h"
#include "sci_const.h"
#include "spectrum_info.h"

int
get_erange (spec_info s,
            double e){

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

int
get_i (spec_info s,
       int from){

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
get_il (spec_info s,
        int from){

  if (from>s->n_states) {
    return -1;
  } else {
    return s->idx_map[from-1];
  }
}

int
get_inext (spec_info s,
           int from){

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
            int from) {
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
eval_trs (spec_info s){

  int j,k,l;
  int n_fs,n_is;
  n_fs = n_is = 0;

  int last_i,curr_i;
  int prog_step;
  int is_num;
  int n_is_tmp = 0;
  int wh = 0; /* write head for the gs2is array */

  screen scr = init_screen(s);
  scr -> is2fs = da_init(10,10);
  scr -> is_idxs = da_init(10,10);
  scr -> gs2is = da_init(10,10);

  double e_gs,e_is,e_fs;
  double tmp_tmax1,tmp_tmax2;
  double tmom_gi,tmom_if;
  double bw;
  double ** tmp_evals = malloc(sizeof(double*)*2);

  for (j=0; j<2; j++) {
    if((tmp_evals[j] = malloc(s->n_trans*sizeof(double))) == NULL ){
      fprintf(stderr, "calc_spec.c:function calc_spec, malloc: failed \
to allocate memory for \"tmp_evals[%d]\"\n",j);
      printf( "program terminating due to the previous error.\n");
      exit(1);
    }
  }

  /* keep check on what states that have been processed to make sure that no
     "to" state appears in the pi matrix more than once  */
  int * proc_st;
  int n_proc = 0;
  int is_proc;
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
    }

    if (curr_i != last_i) {
      if (curr_i != (int)s->trs[0][get_i(s,curr_i)]) {
        return -1;
      }

      s -> idx_map[curr_i-1] = j;
    }

    e_fs = s->trs[3][j];
    e_is = s->trs[2][j];

    /* is this a transition from an intermediate to final state? */
    if (inrange((e_fs - s->e0)*AUTOEV,s->md->state_er[5],s->md->state_er[6]) &&
        inrange((e_is - s->e0)*AUTOEV,s->md->state_er[3],s->md->state_er[4])) {
      tmom_if = s->trs[4][j];

      if (s->trs[5][j] == 1) {
        tmp_tmax1 = s->tmax_d;
      } else {
        tmp_tmax1 = s->tmax_q;
      }

      /* scr the i.state -> f.state transition */
      if ((tmom_if/tmp_tmax1) > s->md->state_t[3]){

        da_append(scr->is2fs,j);
        is_num = s->trs[0][j];

        /* have the gs-is transitions for this particular is already
           been processed? */
        is_proc = intinint(proc_st,is_num,n_proc);
        tmp_evals[1][n_fs++] = (e_fs-s->e0)*AUTOEV;
        tmp_evals[0][n_is++] = (e_is-s->e0)*AUTOEV;

        if (is_proc == -1) {

          /* loop over the list of transitions looking for corresponding
             g.state -> i-state transitions */
          k = 0;
          n_is_tmp = scr->gs2is->n_el;
          while((int)s->trs[0][k] != -1) {
            e_gs = s->trs[2][k];

            if ((((int)s->trs[1][k]) == is_num) &&
                inrange((e_gs - s->e0)*AUTOEV,s->md->state_er[1],\
                        s->md->state_er[2])) {

              tmom_gi = s->trs[4][k];

              /* scr intermediate state transition */
              if (s->trs[5][k] == 1) {
                tmp_tmax2 = s->tmax_d;
              } else {
                tmp_tmax2 = s->tmax_q;
              }

              bw = get_rbdist(s->e0,e_gs);
              if (((tmom_gi/tmp_tmax2) > s->md->state_t[2]) &&
                  (bw >= s->md->state_t[1])){
                da_append(scr->gs2is,k);
              }
            }
            k++;
          }

          da_append(scr->is_idxs,wh);
          wh = scr->gs2is->n_el;
          proc_st[n_proc++] = is_num;

          if (n_is_tmp == scr->gs2is->n_el) {

            /* no gs2is transitions found for that is2fs transition */
            scr->is2fs->n_el--;
            scr->is_idxs->n_el--;
            n_is--;
            n_fs--;
            n_proc--;
          }
        }
        else{

          /* find the intermediate state in the gs2trs matrix */
          for (l=0; l<scr->gs2is->n_el; l++) {
            if ( (int)s->trs[1][scr->gs2is->a[l]] == is_num) {
              break;
            }
          }
          da_append(scr->is_idxs,l);
        }
      }
    }

    last_i = (int)s->trs[0][j];
    j++;

    if (j % prog_step == 0) {
      printf( "        %.2f%%\r", (((float)j/(float)s->n_trans)*100));
    }
  }
  da_append(scr->gs2is,-1);

  scr->emin_x = floor((get_minl(tmp_evals[0],n_is)-1))/AUTOEV;
  scr->emax_x = ceil((get_maxl(tmp_evals[0],n_is) + 1))/AUTOEV;
  scr->emin_y = floor((get_minl(tmp_evals[1],n_fs)-1))/AUTOEV;
  scr->emax_y = ceil((get_maxl(tmp_evals[1],n_fs) + 1))/AUTOEV;

  s -> n_states = n_proc;
  s -> scr = scr;

  if ((scr->is2fs->n_el == 0) || (scr->gs2is->n_el == 0)) {
    fprintf(stderr, "calc_spec.c, : no intermediate or final states states were\
 found in the energy range you provided. (scr->is2fs->n_el = %d,\
 scr->gs2is->n_el = %d)\n",scr->is2fs->n_el,scr->gs2is->n_el);
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

int
eval_trs_mult (spec_info s){

  int j,k,l,m;
  int n_fs,n_is;

  int last_i,curr_i;
  int prog_step;
  int n_sfs; /* number of screened final states */
  int n_sis;
  int is_num = 0;
  int n_is_tmp = 0;
  int is_pos = 0;
  int is_idx = 0;
  int nt = 0; /* number of transitions fround in the provided energy range */
  int nt2=0;
  int wh = 0; /* write head for the gs2is array */

  screen scr = init_screen(s);
  scr -> is2fs = da_init(10,10);
  scr -> is_idxs = da_init(10,10);
  scr -> gs2is = da_init(10,10);
  scr -> ii_start = da_init(10,10);

  screen scr1 = init_screen(s);
  scr1 -> is2fs = da_init(10,10);
  scr1 -> is_idxs = da_init(10,10);
  scr1 -> gs2is = da_init(10,10);
  scr1 -> ii_start = da_init(10,10);

  double e_gs,e_is,e_fs;
  double tmp_tmax1,tmp_tmax2;
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
    }

    if (curr_i != last_i) {
      if (curr_i != (int)s->trs[0][get_i(s,curr_i)]) {
        return -1;
      }
      s -> idx_map[curr_i-1] = j;
    }

    e_fs = s->trs[3][j];
    e_is = s->trs[2][j];

    /* is this a transition from an intermediate to final state? */
    if (inrange((e_fs - s->e0)*AUTOEV,s->md->state_er[5],s->md->state_er[6]) &&
        inrange((e_is - s->e0)*AUTOEV,s->md->state_er[3],s->md->state_er[4])) {

      da_append(scr->is2fs,j);
      is_num = s->trs[0][j];

      /* have the gs-is transitions for this particular is already
         been processed? */
      is_proc = intinint(proc_st,is_num,n_proc);

      if (is_proc == -1) {

        /* loop over the list of transitions looking for corresponding
           g.state -> i-state transitions */
        k = 0;
        n_is_tmp = scr->gs2is->n_el;
        while((int)s->trs[0][k] != -1) {
          e_gs = s->trs[2][k];

          if ((((int)s->trs[1][k]) == is_num) &&
              inrange((e_gs - s->e0)*AUTOEV,s->md->state_er[1],\
                      s->md->state_er[2])) {

            da_append(scr->is_idxs,scr->gs2is->n_el);
            da_append(scr->gs2is,k);
            /* nt++; */
          }
          k++;
        }
        if (n_is_tmp == scr->gs2is->n_el) {
          /* no gs2is transitions found for that is2fs transition */
          scr->is2fs->n_el--;
        }
        else{
          /* gs2is transitions were added */

          da_append(scr->is_idxs,-1);
          da_append(scr->ii_start,wh);

          wh = scr->is_idxs->n_el;
          proc_st[n_proc++] = is_num;
        }
      }
      else {

        /* find the intermediate state in the gs2trs matrix */
        for (l=0; l<scr->is_idxs->n_el; l++) {
          is_idx = scr->is_idxs->a[l];
          if (is_idx != -1) {
            if ((int)s->trs[1][scr->gs2is->a[is_idx]] == is_num) {
              break;
            }
          }
        }
        da_append(scr->ii_start,l);
      }
    }

    last_i = (int)s->trs[0][j];
    j++;

    if (j % prog_step == 0) {
      printf( "        %.2f%%\r", (((float)j/(float)s->n_trans)*100));
    }
  }

  n_sfs = scr->is2fs->n_el;

  /* for (j=0; j<scr->ii_start->n_el; j++) { */
  /*   printf( "%d\n", scr->ii_start->a[j]); */
  /* } */
  /* for (j=0; j<scr->is_idxs->n_el; j++) { */
  /*   printf( "%d\n", scr->is_idxs->a[j]); */
  /* } */


  /* printf( "%d, %d\n", n_sfs, n_sis); */

  /* fprintf(stderr, "\n\n=======Valgrind eject point=======\n\n"); */
  /* exit(1); */
  nt2 = 0;
  for (j=0; j<n_sfs; j++) {

    is_pos = scr->ii_start->a[j]; /* position at is_idx where gs2is
                                     idx will start */
    is_idx = scr->is_idxs->a[is_pos];
    is_num = (int)s->trs[0][scr->is2fs->a[j]];

    /* printf( "\n new vals: %d %d %d %d %d == %d \n",is_pos,is_idx,is_idx-1,is_idx+1,is_num,(int)s->trs[1][scr->gs2is->a[is_idx]]); */

    for (k=is_pos; ((is_idx = scr->is_idxs->a[k]) != -1); k++) {

      /* printf( "g:%d i:%d i:%d f:%d\n",(int)s->trs[0][scr->gs2is->a[is_idx]],(int)s->trs[1][scr->gs2is->a[is_idx]],(int)s->trs[0][scr->is2fs->a[j]],(int)s->trs[1][scr->is2fs->a[j]] ); */
      nt2++;
    }
  }

  for (j=0; j<2; j++) {
    if((int_dist[j] = malloc(nt2*sizeof(double))) == NULL ){
      fprintf(stderr, "transitions.c:function eval_trs, malloc: failed \
to allocate memory for \"int_dist[%d]\"\n",j);
      printf( "program terminating due to the previous error.\n");
      exit(1);
    }
  }

  l = 0;
  int_tot = 0;
  for (j=0; j<n_sfs; j++) {

    is_pos = scr->ii_start->a[j]; /* position at is_idx where gs2is
                                     idx will start */
    is_idx = scr->is_idxs->a[is_pos];
    is_num = (int)s->trs[0][scr->is2fs->a[j]];

    tmom_if = s->trs[4][scr->is2fs->a[j]];

    for (k=is_pos; ((is_idx = scr->is_idxs->a[k]) != -1); k++){

      tmom_gi = s->trs[4][scr->gs2is->a[is_idx]];
      bw = get_rbdist(s->e0,s->trs[2][scr->gs2is->a[is_idx]]);
      int_dist[0][l] = l;

      int_dist[1][l] = tmom_if*tmom_gi*bw;
      int_tot += int_dist[1][l];
      l++;
    }
  }

  /* fflush(stdout); */
  /* printf( "int_dist pre:\n" ); */
  /* for (j=0; j<nt2; j++) { */
  /*   printf( "%le %le\n", int_dist[0][j],int_dist[1][j]); */
  /* } */
  /* sort ascendingly according to energy */
  quicksort_d(int_dist[1],int_dist[0],0,nt2-1,nt2);

  /* printf( "int_dist post:\n" ); */
  /* for (j=0; j<nt2; j++) { */
  /*   printf( "%le %le\n", int_dist[0][j],int_dist[1][j]); */
  /* } */
  /* fflush(stdout); */
  /* fprintf(stderr, "\n\n=======Valgrind eject point=======\n\n"); */
  /* exit(1); */


  /* screen until retaining x% of the intensity of all transitions,
   set all "sceened" states to the n_trans+1 */
  int_scr = int_dist[1][nt2-1];

  for (j=nt2-2; j>=0; j--) {
    if ((int_scr/int_tot) <= 0.99) {
      /* printf( "\nBreaking at state %d/%d  %le %le => %le div = %le\n",j,nt2, int_dist[0][j],int_dist[1][j],int_scr,int_tot); */
      int_scr += int_dist[1][j];
    }
    else {
      int_dist[0][j] = nt2;
    }
  }

  /* printf( "\ntotal:%le partial:%le quota:%le\n", int_tot, int_scr, int_scr/int_tot); */
  /* fflush(stdout); */
  /* fprintf(stderr, "\n\n=======Valgrind eject point=======\n\n"); */
  /* exit(1); */
  /* printf( "PRE\n" ); */
  /* for (j=0; j<nt2; j++) { */
  /*   printf( "%le %le\n", int_dist[0][j],int_dist[1][j]); */
  /* } */
  quicksort_d(int_dist[0],int_dist[1],0,nt2-1,nt2);
  /* printf( "POST\n" ); */
  /* for (j=0; j<nt2; j++) { */
  /*   printf( "%le %le\n", int_dist[0][j],int_dist[1][j]); */
  /* } */
  /* fflush(stdout); */
  /* fprintf(stderr, "\n\n=======Valgrind eject point=======\n\n"); */
  /* exit(1); */

  tmp_evals[0] = malloc(scr->is2fs->n_el*sizeof(double));
  tmp_evals[1] = malloc(scr->is2fs->n_el*sizeof(double));
  /* finally, generate the new screen by only including those states that were not screened out above */
  n_fs = n_is = 0;
  l = m = 0;
  wh = 0;
  for (j=0; j<n_sfs; j++) {

    is_pos = scr->ii_start->a[j]; /* position at is_idx where gs2is
                                     idx will start */
    is_idx = scr->is_idxs->a[is_pos];
    is_num = (int)s->trs[0][scr->is2fs->a[j]];

    n_is_tmp = scr1->is_idxs->n_el;

    for (k=is_pos; ((is_idx = scr->is_idxs->a[k]) != -1); k++){
      /* printf( "checking : %d == %d\n", l,(int)int_dist[0][m]); */
      if (l == (int)int_dist[0][m]) {
        /* printf( "g:%d i:%d i:%d f:%d\n",(int)s->trs[0][scr->gs2is->a[is_idx]],(int)s->trs[1][scr->gs2is->a[is_idx]],(int)s->trs[0][scr->is2fs->a[j]],(int)s->trs[1][scr->is2fs->a[j]] ); */
        /* is there an index in scr1 that already maps to
           this gs2is transition?*/
        if ((is_proc = intinint(scr1->gs2is->a,scr->gs2is->\
                                a[is_idx],scr1->gs2is->n_el)) != -1) {
          da_append(scr1->is_idxs,is_proc);
        }
        else {


          da_append(scr1->is_idxs,scr1->gs2is->n_el);
          da_append(scr1->gs2is,scr->gs2is->a[is_idx]);
        }
        m++; /* start looking for the next state in int_dist */
      }
      l++;
    }
    /* check if any gs2is transitions were added */
    if (n_is_tmp != scr1->is_idxs->n_el) {
      da_append(scr1->is_idxs,-1);
      da_append(scr1->is2fs,scr->is2fs->a[j]);
      da_append(scr1->ii_start,wh);
       fflush(stdout);
      wh = scr1->is_idxs->n_el;
      tmp_evals[0][n_is++] = (s->trs[3][scr->gs2is->a[scr->gs2is->a[scr->is_idxs->a[is_pos]]]] - s->e0)*AUTOEV;
      tmp_evals[1][n_fs++] = (s->trs[3][scr->is2fs->a[j]] - s->e0)*AUTOEV;
    }
  }

  /* printf( "screening finished\n" ); */
  /* fflush(stdout); */
  /* n_sfs = scr1->is2fs->n_el; */
  /* for (j=0; j<n_sfs; j++) { */

  /*   is_pos = scr1->ii_start->a[j]; /\* position at is_idx where gs2is */
  /*                                    idx will start *\/ */
  /*   is_idx = scr1->is_idxs->a[is_pos]; */
  /*   is_num = (int)s->trs[0][scr1->is2fs->a[j]]; */
  /*   printf( "idx:%d\n",is_idx); */
  /*   for (k=is_pos; ((is_idx = scr1->is_idxs->a[k]) != -1); k++) { */
  /*     printf( "idx:%d g:%d i:%d i:%d f:%d\n",is_idx,(int)s->trs[0][scr1->gs2is->a[is_idx]],(int)s->trs[1][scr1->gs2is->a[is_idx]],(int)s->trs[0][scr1->is2fs->a[j]],(int)s->trs[1][scr1->is2fs->a[j]] ); */
  /*   } */
  /* } */
  /* fflush(stdout); */
  /* fprintf(stderr, "\n\n=======Valgrind eject point=======\n\n"); */
  /* exit(1); */
  /* n_sis = scr1->gs2is->n_el; */
  /* for (j=0; j<n_sis; j++) { */
  /*     printf( "g:%d i:%d\n",(int)s->trs[0][scr1->gs2is->a[j]],(int)s->trs[1][scr1->gs2is->a[j]] ); */
  /* } */

  /* for (j=0; j<scr1->is2fs->n_el; j++) { */
  /*     printf( "i:%d f:%d\n",(int)s->trs[0][scr1->is2fs->a[j]],(int)s->trs[1][scr1->is2fs->a[j]]); */
  /* } */
  /* fflush(stdout); */
  /* for (j=0; j<scr1->is_idxs->n_el; j++) { */

  /*   if (scr1->is_idxs->a[j] != -1) { */
  /*     printf( "is_idx:%d => %d\n",scr1->is_idxs->a[j],(int)s->trs[0][scr1->gs2is->a[scr1->is_idxs->a[j]]]); */
  /*       fflush(stdout); */
  /*   } */
  /*   else { */
  /*     printf( "is_idx:%d \n",scr1->is_idxs->a[j]); */
  /*       fflush(stdout); */
  /*   } */
  /* } */

  /* fflush(stdout); */
  /* fprintf(stderr, "\n\n=======Valgrind eject point=======\n\n"); */
  /* exit(1); */
  scr1->emin_x = floor((get_minl(tmp_evals[0],n_is)-1))/AUTOEV;
  scr1->emax_x = ceil((get_maxl(tmp_evals[0],n_is) + 1))/AUTOEV;
  scr1->emin_y = floor((get_minl(tmp_evals[1],n_fs)-1))/AUTOEV;
  scr1->emax_y = ceil((get_maxl(tmp_evals[1],n_fs) + 1))/AUTOEV;

  s -> n_states = n_proc;
  s -> scr = scr1;
  free_screen(scr);

  if ((scr1->is2fs->n_el == 0) || (scr1->gs2is->n_el == 0)) {
    fprintf(stderr, "calc_spec.c, : no intermediate or final states states were\
 found in the energy range you provided. (scr1->is2fs->n_el = %d,\
 scr->gs2is->n_el = %d)\n",scr->is2fs->n_el,scr->gs2is->n_el);
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
count_states (spec_info s){

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

  trs_el[0][l] = s->trs[0][l] = -1;

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
