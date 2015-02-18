#ifndef PARSE_INPUT_H
#define PARSE_INPUT_H

extern int n_states; /* number of states used in the calculation */
extern int n_trans; /* number of electronic transitions involved */
extern double * input_data[4]; /* defined in parse_input.c */
extern int ** state_indices; /* post-screening state indices */

struct e_state{

  int list_idx;
  int state_idx;
  int type; /* 1=ground state, 2=intermediate state, 3=final state */

  /* indices of states in transitions occuring from this state*/
  int * idxs_to;

  double bw; /* boltzmann weight */
  double e_val; /* energy */

   /* transition moments values for each transition found in idxs_to */
  double * t_moms;

  /* this list is doubly linked. */
  struct e_state * next;
  struct e_state * last;
};

/* each input file loaded into the program has an associated info node
 containing general data about that specific input file*/
struct info_node{
  int n_gs;
  int n_is;
  int n_fs;
  struct info_node * next;
};

/* function init_state_node

   * synopsis:
   by default, a call to the init_root_node function also creates a general info node for one data tree.

   * algorithm:

   * input:

   * output:

   * side-effects:

   */
struct e_state *
init_root_node ();

/* function init_state_node

   * synopsis:

   * algorithm:

   * input:

   * output:

   * side-effects:

   */
struct e_state *
init_state_node ();

/* function init_state_ll

   * synopsis:

   * algorithm:

   * input:

   * output:

   * side-effects:

   */
struct e_state *
init_state_ll (double ** parsed_input
              );


/* function get_bdist

   * synopsis:

   * algorithm:

   * input:

   * output:

   * side-effects:

   */
double
get_bdist (double e_val,
           double temp);

/* function sort_states

 * synopsis:
 The sort_states function is used to reduce the number of states used
 for generating the RIXS map. Through analyzing the boltzmann distribution
 of initial states, and screening of intensities for transitions,
 states that are unlikely to contribute to the RIXS process are disregarded.
 This effectively lowers the time needed to generate a given RIXS map. The
 thresholdsused to exclude states from the calculation can be used to find
 an optimal compromise between number of states included and the time needed
 to run the program.

 * algorithm:
 reads all initial states (k) contained in the array parsed_input[2]. from the
 these energy values, the boltzmann distribution is calculated. All k states
 with a weight value below the thrsh_b value gets excluded from the RIXS map
 calculation.

 * input:

 * output:

 * side-effects:

 */
int
sort_states (double ** parsed_input);

/* function getinput_molcas:

   * algorithm:

   * input:

   * output:

   * side-effects:

   */

double **
parse_input_molcas (char * fn_infile
                    );


/* function parse_input:
 * synopsis:
 parse_input is the output data file interface to rmap. For a given output type,
 be it from Molcas, RACAH or other electron energy calculation programs,
 an input parsing function needs to be defined. This function has to
 return a multi-dimensional matrix "parsed_input", containing data as follows
 ([row index][column index], ":" read as "all indexes"):

 parsed_input[0][:] : idexes of initial states for a given transition
 parsed_input[1][:] : idexes of final states for a given transition
 parsed_input[2][:] : energy eigenvalues, relative to the lowest energy ground
 state of the system, for the indexes in parsed_input[0]
 parsed_input[3][:] : -||- , for the indexes in parsed_input[1]
 parsed_input[4][:] : transition moment values for each transition extracted
 from the input

 * algorithm:

 * input:

 * output:

 * side-effects:

 */

int
parse_input (char * fn_infile,
                   int fn_len);

#endif /* PARSE_INPUT_H */
