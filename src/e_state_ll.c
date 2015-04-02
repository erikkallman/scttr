#include <stdlib.h>
#include <stdio.h>
#include "e_state_ll.h"
#include "rmap_structs.h"
#include "info_ll.h"
#include "std_f.h"
#include "sci_const.h"
#include "std_num_ops.h"

int
is_state_inlist (info_node inode,
                 int idx
                 ) {

  e_state next_state = inode -> root_e_state;
  e_state curr_state;

  while(next_state != NULL){
    curr_state = next_state;
    if (curr_state -> state_idx == idx) {
      return 1;
    }

    next_state = curr_state -> next;
  }
  return 0;
}

void
swapd_estate(e_state e1,
             e_state e2){

  e2 -> next = e1 -> next;
  e2 -> last = e1 -> last;
  e2 -> info = e1 -> info;
  e2 -> list_idx = e1 -> list_idx;

  dstruct_estate(e1);
}

void
dstruct_estate(e_state e){

  free(e -> idxs_to);
  free(e -> t_moms);
  free(e -> e_vals);
  free(e);
  e = NULL;
}

double
get_trans (e_state es, /* root of the electronic state llist */
           int idx_to /* index of the state transitioning to */
           ) {
  /*                if (a1 == NULL) { */
  /*   fprintf(stderr, "Found NULL\n"); */
  /*   printf( "program terminating due to the previous error.\n"); */
  /*   exit(1); */
  /* } */
  int j; /* looping variables */
  int trs_max = es -> n_tfrom;

  int * ti = es -> idxs_to;
  /* printf( "looping\n" ); */
  for (j=0; j<trs_max; j++) {
    /* printf( "in state %d: to = %d, from = %d, %d\n", es -> state_idx,idx_to, ti[j], trs_max); */
    if (ti[j] == idx_to) {
      return (es -> t_moms)[j];
    }
  }
  /* we looped over the entire list of electronic states but didnt find the
   one that was requested */
  fprintf(stderr, "e_state_ll.c, function get_trans: unable to locate state of \
index %d in the list of transitions for state %d.\n", idx_to, es -> state_idx);
  e_state2s(es, 1);
  printf( "program terminating due to the previous error.\n");
  exit(1);
}

double
get_eval (e_state es, /* root of the electronic state llist */
           int idx_to /* index of the state transitioning to */
           ) {
  int j; /* looping variables */
  int trs_max = es -> n_tfrom;

  int * ti = es -> idxs_to;
  /* printf( "looping\n" ); */
  for (j=0; j<trs_max; j++) {
    /* printf( "to = %d, from = %d, %d\n",idx_to, ti[j], trs_max); */
    if (ti[j] == idx_to) {
      return (es -> e_vals)[j];
    }
  }
  /* we looped over the entire list of electronic states but didnt find the
   one that was requested */
  fprintf(stderr, "e_state_ll.c, function get_trans: unable to locate state of \
index %d in the list of transitions for state %d.\n", idx_to, es -> state_idx);
  printf( "program terminating due to the previous error.\n");
  exit(1);
}

e_state
get_state_si (info_node inode, /* the info node at root of the state ll */
           int s_idx /* index of the state to get */
           ){

  int n_s = inode -> n_states;
  int n_t = inode -> n_trans;

  e_state curr_st = (inode -> root_e_state);
  e_state next_st = curr_st;

  while((curr_st = next_st) != NULL){

    if ((curr_st -> state_idx) == s_idx) {
      return curr_st;
    }
    next_st = curr_st -> next;
  }

  /* unable to locate state of in the list of electronic states */
  return NULL;
}

e_state
get_state_li (info_node inode, /* the info node at root of the state ll */
              int l_idx /* list index of the state to get */
              ){

  int n_s = inode -> n_states;
  int n_t = inode -> n_trans;

  e_state curr_st = (inode -> root_e_state);
  e_state next_st = curr_st;

  while((curr_st = next_st) != NULL){

    if ((curr_st -> list_idx) == l_idx) {
      return curr_st;
    }
    next_st = curr_st -> next;
  }

  /* unable to locate state of in the list of electronic states */
  return NULL;
}

double
get_ediff (info_node inode, /* root of the electronic state llist */
           int idx_es1,
           int idx_es2
           ) {
  int j;
  e_state es = get_state_si(inode, idx_es1);
  int trs_max = es -> n_tfrom;
  /* printf( "\neval = %le\n", (es -> e_val)); */
  int * ti = es -> idxs_to;
  /* printf( "looping\n" ); */
  for (j=0; j<=trs_max; j++) {
    /* printf( "%le\n", (es -> e_vals)[j]); */
    if (ti[j] == idx_es2) {
      return ((es -> e_vals)[j] - (es -> e_val));
    }
  }
  fprintf(stderr, "e_state_ll.c, get_ediff: unable to locate state of \
index %d in the list of transitions from state %d.\n", idx_es2, idx_es1);
  printf( "program terminating due to the previous error.\n");
  exit(1);
}

void
e_statelist2s(info_node inode,
              int flag){

  int j;

  char * fn_in = inode -> str_id;
  int n_s = inode -> n_states;
  int n_t = inode -> n_trans;
  double bw_sum = inode -> bw_sum;
  double gs_eval =  (inode -> root_e_state) -> e_val;
  e_state curr_st = (inode -> root_e_state);
  e_state next_st = curr_st;
  printf( "\n  -printing the content of the electronic state list:");

  for (j=0; j<n_s; j++) {

    curr_st = next_st;
    if (flag == 1) {
      e_state2s(curr_st,1);
    } else {

      printf( "\n   state[%d/%d], type %d\n", curr_st->list_idx + 1, n_s, curr_st->type);
      printf( "     state_idx = %d\n", curr_st->state_idx);

      if ((curr_st->type) != 2) {
        printf( "     boltzmann weight = %le\n", curr_st->bw);
      }

      printf( "     n_tfrom = %d\n", curr_st->n_tfrom);
      printf( "     e_val(au) = %le\n", curr_st->e_val);
      printf( "     max_tmom = %le\n", curr_st->max_tmom);
      printf( "     delta e(ev) = %le\n", (gs_eval - (curr_st->e_val))*(double)AUTOEV);

    }
    next_st = curr_st -> next;
  }
}

void
  e_state2s(e_state es,
            int flag){
  int j;

  printf( "\n   state[%d], type %d\n", es->list_idx + 1, es->type);
  printf( "     state_idx = %d\n", es->state_idx);

  if ((es->type) != 2) {
    printf( "     boltzmann weight = %le\n", es->bw);
  }

  printf( "     n_tfrom = %d\n", es->n_tfrom);
  printf( "     e_val(au) = %le\n", es->e_val);
  printf( "     max_tmom = %le\n", es->max_tmom);
  if (flag == 1) {
    for (j=0; j<es->n_tfrom; j++) {
      printf( "     -> %d e = %le, tmom = %le\n", es->idxs_to[j], es->e_vals[j], es->t_moms[j]);
    }
  }
}

void
reset_info_maxvals (info_node inode) {

  double maxt_fs, maxt_is;
  double tmp_max;

  double * trans;

  e_state curr_st = (inode -> root_e_state);
  e_state next_st = curr_st;
  inode -> mt_is = 0;
  inode -> mt_fs = 0;

  maxt_fs = -0.1;
  maxt_is = -0.1;

  while((curr_st = next_st) != NULL){
    tmp_max = get_maxl(curr_st -> t_moms, curr_st -> n_tfrom);
    curr_st -> max_tmom = tmp_max;
    if ((curr_st -> state_idx) != 2) {
      if (tmp_max >= maxt_fs) {
        maxt_fs = tmp_max;
      }
    }
    else {
      if (tmp_max >= maxt_is) {
        maxt_is = tmp_max;
      }
    }
    next_st = curr_st -> next;
  }

  inode -> mt_is = maxt_is;
  inode -> mt_fs = maxt_fs;

}
