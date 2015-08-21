#ifndef ESTATE_S_H
#define ESTATE_S_H

/* the structs used by the program are not logically separable and are all
   mutually dependent, which is why they share the same scope */
struct info_node_s;
typedef struct info_node_s * info_node;

struct spec_s;
typedef struct spec_s * spec;

static int n_inodes;

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
  int * mom_types; /* types of momenta in the plot 1/2/3 di/quad/octo-pole */
  /* sum of the boltzmann weights of all states in the system*/
  double ** mt; /* pointers to the maximum transition moment for the
                     intermediate and final/ground states */
  double bw_sum; /* sum of all boltz  */

  char * str_id; /* the input file name identifying this info node */

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
