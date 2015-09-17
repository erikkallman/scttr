#ifndef ESTATE_S_H
#define ESTATE_S_H

/* the structs used by the program are not logically separable and are all
   mutually dependent, which is why they share the same scope */
struct spec_info_s;
typedef struct spec_info_s * spec_info;

struct spec_s;
typedef struct spec_s * spec;

struct metadata_s; /* metadata container */
typedef struct metadata_s * metadata;

/* each input file loaded into the program has an associated spectrum information node
 containing all data about that specific input file*/
struct spec_info_s{

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

  metadata s_md; /* spectrum information node metadata */

  spec root_spec;

  spec_info next;
  spec_info last;

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


struct metadata_s{
  /* the metadata struct contains all data that is read from the command line
     input, like file paths, energy ranges, and so on. */

  /* path to .. */
  char * outpath; /* smap data file */
  char * inpath; /* plot output file */
  char * inp_fn; /* input filename */
  char * inp_sfx; /* input file suffix */

  double * state_er; /* expected range of energie eigenvalues in the input */
  double * state_t; /* transition intensity screening thresholds */
  double * res; /* spectral resolution (eV) of the produced spectrum */
  double * fwhm; /* full-width half-maximum value for broadening the peaks */
};
#endif /* ESTATE_S_H */
