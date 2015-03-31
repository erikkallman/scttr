#ifndef E_STATE_LL_H
#define E_STATE_LL_H
#include "rmap_structs.h"
#include "info_ll.h"

/* function is_state_inlist

   * synopsis:

   * algorithm:

   * input:

   * output:

   * side-effects:

   */
int
is_state_inlist (info_node inode,
                 int idx
                 );

/* function swapd_estate
   destructive swap. swaps state e1 with e2 and destructs e2
   * synopsis:

   * algorithm:

   * input:

   * output:

   * side-effects:

   */
void
swapd_estate(e_state e1,
             e_state e2);

void
dstruct_estate(e_state e);

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

double
get_eval (e_state es_root, /* root of the electronic state llist */
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
e_statelist2s(info_node inode,
              int flag);

void
e_state2s(e_state es,
          int flag);

/* function reset_info_maxvals

   * synopsis:

   * algorithm:

   * input:

   * output:

   * side-effects:

   */
void
reset_info_maxvals (info_node inode);

#endif /* E_STATE_LL_H */
