#ifndef E_STATE_LL_H
#define E_STATE_LL_H
#include "structs.h"

/* function set_ttypes

   * synopsis:

   * algorithm:

   * input:

   * output:

   * side-effects:

   */
void
set_ttypes (double * state_er,
            info_node inode);

/* function sort_states

   * synopsis:

   * algorithm:

   * input:

   * output:

   * side-effects:

   */
void
sort_states (double * state_er,
             double * a,
             int ** groups,
             int n_vals);

/* function set_estate

   * synopsis:
   uses the parsed_input defined in parse_input to declare the variables
   in the linked list of states.
   * algorithm:

   * input:

   * output:

   * side-effects:

   */
int
set_estate (estate st,
            int s_idx,
            int * idxs_buf,
            double * evals_buf,
            double * moms_buf,
            int n_trs_from,
            double e_rel,
            double e
            );

int
set_estate_list (double * state_er,
                 double ** parsed_input,
                 int n_states,
                 int n_trans,
                 char * id);

/* function init_estate_list

   * synopsis:

   * algorithm:

   * input:

   * output:

   * side-effects:

   */
estate
init_estate_list (char * str_id,
                  int n_states,
                  int n_trans
                  );

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
swapd_estate(estate e1,
             estate e2);

void
dstruct_estate(estate e);

/* function get_trans

   * synopsis:

   * algorithm:

   * input:

   * output:

   * side-effects:

   */
double
get_trans (estate es_root, /* root of the electronic state llist */
           int idx_to /* index of the state transitioning to */
           );

double
get_eval (estate es_root, /* root of the electronic state llist */
           int idx_to /* index of the state transitioning to */
           );

/* function get_state_si

   * synopsis:
   Get a state from the list based on its state index, obtained from the input to the program.
   * algorithm:

   * input:

   * output:

   * side-effects:

   */
estate
get_state_si (info_node inode,
              int s_idx /* state index of the state to get */
              );


/* function get_state_si

   * synopsis:
   Get a state from the list based on its list index.
   * algorithm:

   * input:

   * output:

   * side-effects:

   */
estate
get_state_li (info_node inode,
              int l_idx /* list index of the state to get */
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
e_state2s(estate es,
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
