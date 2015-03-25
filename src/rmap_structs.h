#ifndef RMAP_STRUCTS_H
#define RMAP_STRUCTS_H

struct e_state_s;
typedef struct e_state_s * e_state;

struct info_node_s;
typedef struct info_node_s * info_node;

struct e_state_s{

  int list_idx;
  int state_idx;

  /* 0=not yet assigned, 1=ground state, 2=intermediate state, 3=final state,
   23= either 2 or 3 */
  int type;

  double bw; /* boltzmann weight */
  double max_tmom;
  double e_val; /* energy of the state on the node (the "from" energy) */

  /* number of transitions from this state */
  int n_tfrom;

  /* indices of states in transitions occuring from this state*/
  int * idxs_to;

   /* transition moments values for each transition found in idxs_to */
  double * t_moms;

  /* energy eigenvalues for the "to" state */
  double * e_vals;

  /* this list is doubly linked. */
  info_node info;
  e_state next;
  e_state last;
};

typedef struct e_state_s * e_state;

/* each input file loaded into the program has an associated info node
 containing general data about that specific input file*/
struct info_node_s{

  int idx; /* default index value */

  int n_states; /* number of electronic states */
  int n_trans; /* number of transitions */

  int n_gs;
  int n_is;
  int n_fs;

  /* sum of the boltzmann weights of all states in the system*/
  double mt_is; /* maximum transition moment for the intermediate states */
  double mt_fs; /* maximum transition moment for all final states */
  double bw_sum; /* sum of all boltz  */

  char * str_id; /* the input file name identifying this info node */
  e_state root_e_state;
  info_node next;
  info_node last;
};

typedef struct info_node_s * info_node;

#endif /* RMAP_STRUCTS_H */

  int idx;
  int n_states;
  int n_trans;

  int n_os;
  int n_is;

  double mt_is;
  double mt_fs;
  double bw_sum;

  char * str_id; /* the input file name identifying this info node */
  e_state root_e_state;
