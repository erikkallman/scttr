#ifndef PARSE_INPUT_H
#define PARSE_INPUT_H

extern int n_states; /* number of states used in the calculation */
extern int n_trans; /* number of electronic transitions involved */
extern double * input_data[4]; /* defined in parse_input.c */
extern int ** state_indices; /* post-screening state indices */

/* each input file loaded into the program has an associated info node
 containing general data about that specific input file*/
struct info_node_s{

  int idx; /* default index value */

  int n_gs;
  int n_is;
  int n_fs;

  char * str_id; /* the input file name identifying this info node */
  struct e_state * root_e_state;
  struct info_node_s * next;
  struct info_node_s * last;
};

typedef struct info_node_s * info_node;

struct e_state{

  int list_idx;
  int state_idx;

  /* 0=not yet assigned, 1=ground state, 2=intermediate state, 3=final state */
  int type;

  double bw; /* boltzmann weight */
  double e_val; /* energy of the state on the node (the "from" energy) */

  /* indices of states in transitions occuring from this state*/
  int * idxs_to;

   /* transition moments values for each transition found in idxs_to */
  double * t_moms;

  /* energy eigenvalues for the "to" state */
  double * e_vals;

  /* this list is doubly linked. */
  info_node info;
  struct e_state * next;
  struct e_state * last;
};


/* function init_info_node

   * synopsis:

   * algorithm:

   * input:

   * output:

   * side-effects:

   */
info_node
init_info_node (char * s);

/* function set_state_node

   * synopsis:
   uses the parsed_input defined in parse_input to declare the variables
   in the linked list of states.
   * algorithm:

   * input:

   * output:

   * side-effects:

   */
int
set_state_node (struct e_state * st,
                int s_idx,
                int * idxs_buf,
                double * evals_buf,
                double * moms_buf,
                int n_trs_from,
                double e_rel,
                double e
                );

/* function set_state_ll

   * synopsis:
   uses the parsed_input defined in parse_input to declare the variables
   in the linked list of states.
   * algorithm:

   * input:

   * output:

   * side-effects:

   */
int
set_state_ll (double ** parsed_input,
              char * id);


/* function init_state_ll

   * synopsis:

   * algorithm:

   * input:

   * output:

   * side-effects:

   */
struct e_state *
init_state_ll (char * str_id);


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
