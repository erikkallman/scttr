#ifndef ESTATE_S_H
#define ESTATE_S_H

/* the structs used by the program are not logically separable and are all
   mutually dependent, which is why they share the same scope */

struct e_state_s;
typedef struct e_state_s * estate;

struct info_node_s;
typedef struct info_node_s * info_node;

struct spec_s;
typedef struct spec_s * spec;

static int n_inodes;

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

   /* transition moments types for each transition found in idxs_to */
  double * ttypes;

  /* energy eigenvalues for the "to" state */
  double * e_vals;

  /* this list is doubly linked. */
  info_node info;
  estate next;
  estate last;
};

/* each input file loaded into the program has an associated info node
 containing general data about that specific input file*/
struct info_node_s{

  int idx; /* default index value */

  int n_states; /* number of electronic states */
  int n_trans; /* number of transitions */
  int n_spec; /* number of spectra */

  int n_gfs;
  int n_is;

  int * ev_idxs; /* indices of what states reside in what element of e_vals */

  /* sum of the boltzmann weights of all states in the system*/
  double mt_is; /* maximum transition moment for the intermediate states */
  double mt_fs; /* maximum transition moment for all final states */
  double bw_sum; /* sum of all boltz  */

  char * str_id; /* the input file name identifying this info node */

  estate root_e_state;
  estate * e_states;

  spec root_spec;

  info_node next;
  info_node last;

};

struct spec_s{

  int idx;
  int layer;
  int n_layers;
  int height;
  int length;

  double ** sdat; /* the actual spec data */

  /* next spec for the information node */
  spec next;
  spec last;

  /* next layer of this spec */
  spec next_l;
  spec last_l;
  spec root_l;
};

#endif /* ESTATE_S_H */
