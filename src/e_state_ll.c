#include <stdlib.h>
#include <stdio.h>
#include "e_state_ll.h"
#include "rmap_structs.h"
#include "info_ll.h"


double
get_trans (e_state es, /* root of the electronic state llist */
           int idx_to /* index of the state transitioning to */
           ) {
  int j; /* looping variables */
  int trs_max = es -> n_tfrom;

  int * ti = es -> idxs_to;
  /* printf( "looping\n" ); */
  for (j=0; j<=trs_max; j++) {
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

e_state
get_state (info_node inode, /* the info node at root of the state ll */
           int s_idx /* index of the state to get */
           ){
  e_state curr_st = (inode -> root_e_state);
  e_state next_st = curr_st;

  while((curr_st = next_st) != NULL){

    if ((curr_st -> state_idx) == s_idx) {
      return curr_st;
    }
    next_st = curr_st -> next;
  }
  fprintf(stderr, "e_state_ll.c, function get_state: unable to locate state of \
index %d in the list of electronic states.\n", s_idx);
  printf( "program terminating due to the previous error.\n");
  exit(1);
}

double
get_ediff (info_node inode, /* root of the electronic state llist */
           int idx_es1,
           int idx_es2
           ) {
  int j;
  e_state es = get_state(inode, idx_es1);
  int trs_max = es -> n_tfrom;

  int * ti = es -> idxs_to;
  /* printf( "looping\n" ); */
  for (j=0; j<=trs_max; j++) {
    if (ti[j] == idx_es2) {
      return ((es -> e_val) - (es -> e_vals)[j]);
    }
  }
  fprintf(stderr, "e_state_ll.c, get_editt: unable to locate state of \
index %d in the list of electronic states.\n", idx_es1);
  printf( "program terminating due to the previous error.\n");
  exit(1);
}
