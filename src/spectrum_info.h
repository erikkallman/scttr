#ifndef SPEC_INFO_H
#define SPEC_INFO_H
#include "dyn_array.h"

/* suffixes for the output and input files */
extern const char * dat_sfx;
extern const char * dat_sfx;
extern const char * plot_sfx;
extern const char * log_sfx;
extern const char * bin_sfx;
extern const char * tmp_sfx;
extern const char * log_sfx;

/* All the following struct are documented below  */
struct spec_info_s;
typedef struct spec_info_s * spec_info;

struct metadata_s;
typedef struct metadata_s * metadata;

struct screen_s;
typedef struct screen_s * screen;

/* Global variable pointing to the first spectrum information node
 in the linked list of nodes. */
spec_info root_sinfo;
int n_sinfo; /* total number of information nodes */

/* The screen_s struct:
   This struct contains data defined in the eval_trs function.
   It is used in the calc_spec function to increase the
   execution speed of the rixs_map calculation loop. After all transitions have
   been read from the input, the transitions are screened such that only
   the transitions that are strong enough, compared to the treshold provided
   by the user, are included in the calculation of the RIXS spectrum */
struct screen_s {

  /* the 0th indexed screen contains all states, unscreened, but sorted
   into the energy ranges provided by the user */
  int idx;

  /* maximum and minimum energies in the range of... */
  double emin_x,emax_x; /* .. final states */
  double emin_y,emax_y; /* .. intermediate states. */
  /* the default values of the above variables are defined from user input
   (see the init_screen function)*/

  /* As a result of a call to the eval_trs function, two dynamic arrays (da)
   are defined, that each respectively contains indices of a screened transition
  inside the trs array from ..  */
  da is2fs; /* .. intermediate to final state */
  da gs2is; /* .. ground to intermediate state */

  /* For a given intermediate (I) to final state transition in is2fs, transitions
   from a given ground state to that intermediate state I are located in the
  gs2is matrix, starting from the index is_idxs[j], where j is the index of the
  corresponding final state transition in is2fs. is_idxs is a mapping between
  transitions in is2fs and gs2is which is used to avoid having to search
  through gs2is for the right intermediate state index. */
  da is_idxs;
  da ii_start; /* intermediate state transition index start */
};

/* metadata_s:
  The metadata_s struct contains all data that is read from the command line
  input, like file paths, energy ranges, and so on, which the user provides. */
struct metadata_s{

  int sz_inp; /* size (in bytes) of the input file */
  int so_enrg; /* if == 1, the program reads spin-orbit energies,
                  if not, reads spin-free*/

  /* path to .. */
  char * outpath; /* .. output directory */
  char * inpath; /* .. input directory */

  char * inp_fn; /* input filename */
  char * inp_sfx; /* input file suffix */

  double * state_er; /* user-provided ranges of energy eigenvalues in the input */
  double * state_t; /* transition intensity screening thresholds */
  double * res; /* spectral resolution (eV) of the produced spectrum */
  double * fwhm; /* full-width half-maximum value for broadening the peaks */
};

/* the spec_info_s struct:
   Each input file loaded into the program has an associated spectrum
   information node containing all data about, and extracted from, that
   specific input file. That node is stored in a spec_info_s struct. */
struct spec_info_s{

  int idx; /* index of node in the list value */

  int n_states; /* number of electronic states */
  int n_trans; /* number of transitions between those states */
  int n_spec; /* number of spectra contained in this node */
  int n_tmax; /* maximum number of intermediate state transitions from
               any given electronic state */

  int n_gfs; /* number of ground and final electronic states */
  int n_is; /* number of intermediate electronic states */

  int * idx_map; /* a mapping between a state of a given index, to the position
                  of its transitions stored in the trs variable (see below) */

  /* sum of the boltzmann weights of all states in the system*/
  double bw_sum; /* sum of all boltz  */
  double tmax_q,tmax_d; /* maximum transition moment of all quadrupole (q) and
                           dipole (tmax_d) transitions */
  double e0; /* lowest energy eigenstate */
  double sfac; /* scaling factor for the spectrum */
  /* char * str_id; */ /* the input file name identifying this info node */

  /* the trs variable contains an array of six rows and n_trans columns, each
     row stores the following data, for a transition from state x to state y:
     [0]: index of state x
     [1]: index of state y
     [2]: energy eigenvalue of state x
     [3]: energy eigenvalue of state y
     [4]: transition moment for the transition from x to y
     [5]: type of transition, 1=dipole, 2=quadrupole
  */
  double ** trs;

  metadata md; /* spectrum information node metadata */
  screen scr; /* transition screening information */

  /* spec root_spec; */
  spec_info next;
  spec_info last;

};


metadata
init_md ();

int
free_md (metadata md);


spec_info
get_sinfo (char * idx);


screen
init_screen (spec_info s);

int
free_screen (screen scr);

spec_info
init_sinfo (metadata md);

int
free_sinfo (spec_info s);
#endif /* SPEC_INFO_H */
