#ifndef E_STATE_LL_H
#define E_STATE_LL_H
#include "rmap_structs.h"
#include "info_ll.h"

/* function get_trans

   * synopsis:

   * algorithm:

   * input:

   * output:

   * side-effects:

   */
double
get_trans (e_state es_root, /* root of the electronic state llist */
           int idx_to /* index of the state transitioning to */
           );

/* function get_state

   * synopsis:

   * algorithm:

   * input:

   * output:

   * side-effects:

   */
e_state
get_state (info_node inode,
           int s_idx /* index of the state to get */
           );

double
get_ediff (info_node inode, /* root of the electronic state llist */
           int idx_es1,
           int idx_es2
           );

void
e_state2s(info_node inode);
#endif /* E_STATE_LL_H */
