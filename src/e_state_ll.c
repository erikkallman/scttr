#include <stdlib.h>
#include <stdio.h>
#include "e_state_ll.h"
#include "rmap_structs.h"
#include "info_ll.h"
#include "sci_const.h"

double
get_trans (e_state es, /* root of the electronic state llist */
           int idx_to /* index of the state transitioning to */
           ) {
  int j; /* looping variables */
  int trs_max = es -> n_tfrom;

  int * ti = es -> idxs_to;
  /* printf( "looping\n" ); */
  for (j=0; j<trs_max; j++) {
    /* printf( "to = %d, from = %d, %d\n",idx_to, ti[j], trs_max); */
    if (ti[j] == idx_to) {
      return (es -> t_moms)[j];
    }
  }
  /* we looped over the entire list of electronic states but didnt find the
   one that was requested */
  fprintf(stderr, "e_state_ll.c, function get_trans: unable to locate state of \
index %d in the list of transitions for state %d.\n", idx_to, es -> state_idx);
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
get_state (info_node inode, /* the info node at root of the state ll */
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

double
get_ediff (info_node inode, /* root of the electronic state llist */
           int idx_es1,
           int idx_es2
           ) {
  int j;
  e_state es = get_state(inode, idx_es1);
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
e_state2s(info_node inode){

  int j;

  char * fn_in = inode -> str_id;
  int n_s = inode -> n_states;
  int n_t = inode -> n_trans;

  double gs_eval =  (inode -> root_e_state) -> e_val;
  e_state curr_st = (inode -> root_e_state);
  e_state next_st = curr_st;
  printf( "\n  -printing the content of the electronic state list:");

  for (j=0; j<n_s; j++) {
    curr_st = next_st;
    printf( "\n   state[%d/%d], type %d\n", curr_st->list_idx + 1, n_s, curr_st->type);
    printf( "     state_idx = %d\n", curr_st->state_idx);

    if ((curr_st->type) != 2) {
      printf( "     boltzmann weight = %le\n", curr_st->bw);
    }

    printf( "     n_tfrom = %d\n", curr_st->n_tfrom);
    printf( "     e_val(au) = %le\n", curr_st->e_val);
    printf( "     max_tmom = %le\n", curr_st->max_tmom);
    printf( "     delta e(ev) = %le\n", (gs_eval - (curr_st->e_val))*(double)AUTOEV);
    next_st = curr_st -> next;
  }
}
