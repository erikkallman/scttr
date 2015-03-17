#ifndef PARSE_INPUT_H
#define PARSE_INPUT_H
#include "dynarray.h"
#include "rmap_structs.h"
/* function get_inode

   * synopsis:

   * algorithm:

   * input:

   * output:

   * side-effects:

   */
info_node
get_inode (char * fn_infile);

/* function screen_states
   The screen_states function is used to reduce the number of states used
   for generating the RIXS map. Through analyzing the boltzmann distribution
   of initial states, and screening of intensities for transitions,
   states that are unlikely to contribute to the RIXS process are disregarded.
   This effectively lowers the time needed to generate a given RIXS map. The
   thresholds used to exclude states from the calculation can be used to find
   an optimal compromise between number of states included and the time needed
   to run the program.

   * algorithm:
   1. For each initial state, if its relative boltzmann weight is above threshold,
   include its state index in the array of screened initial states.

   2. For each transition from that ground state, if its transition moment
   relative to the maximum transition moment is above the intermediate state
   threshold, include that intermediate state in the screened intermediate
   states specific to that ground state.

   3. For each intermediate state that passes stage 2 above, screen its final
   state transitions if that has not already been done in subsequent screening
   iterations.

   * input:
   - n_args: the number of threshold values included in the function call.
   the function uses three screening stages (one for the ground, intermediate
   and final states respectively) and one threshold value for each stage.
   for each excluded threshold, a default value is used instead.
   -(...): a number of double values for the thresholds in the order specified above.
   * output:

   * side-effects:

   */
mdda_s *
screen_states (char * fn_infile,
               int n_args,
             ...);

/* function reduce_input

   * synopsis:

   * algorithm:

   * input:

   * output:

   * side-effects:

   */
double **
reduce_input (mdda_s * sidxs /* state indices */
              );


int
init_data_branch(double ** pi, /* parsed input */
                 int ns, /* n states */
                 int nt, /* n transitions */
                 char * fs);

/* function init_info_node

   * synopsis:

   * algorithm:

   * input:

   * output:

   * side-effects:

   */
info_node
init_info_node (char * s,
                int ns,
                int nt);

/* function init_state_ll

   * synopsis:

   * algorithm:

   * input:

   * output:

   * side-effects:

   */
e_state
init_state_ll (char * str_id,
               int n_states,
               int n_trans
               );

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
set_state_node (e_state st,
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
              int n_states,
              int n_trans,
              char * id /* info node id string */
              );


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
/* int */
/* sort_states (double ** parsed_input); */

/* function getinput_molcas:

   * algorithm:

   * input:

   * output:

   * side-effects:

   */

int
parse_input_molcas (char * fn_infile
                    );


/* function parse_input:
 * synopsis:
 parse_input is the output data file interface to rmap. For a given output type,
 be it from Molcas, RACAH or other electron energy calculation programs,
 an input parsing function needs to be defined. this function will need to
 define; 1. the number of states (n_states); 2.transitions (n_trans) in found in
 the input file as well as; 3.a matrix (here called "parsed_input"), of 5 rows
 and n_trans columns, where:

 parsed_input[0][:] : idexes of initial states for a given transition
 parsed_input[1][:] : idexes of final states for a given transition
 parsed_input[2][:] : energy eigenvalues, relative to the lowest energy ground
 state of the system, for the indexes in parsed_input[0]
 parsed_input[3][:] : -||- , for the indexes in parsed_input[1]
 parsed_input[4][:] : transition moment values for each transition extracted
 from the input

 this matrix and the state and transition numbers are used, at the end of the
 function as arguments for the function constructing the linked list of data.
 although this might seem like an overtly complex way of extracting data from
 the input, a developer that wants to add input parsing functionality for a
 different output format only needs to write code for defining the three
 variables described above. note: the first series of initial state values in
 the matrix is expected to contain data for the lowest energy ground state.

 * algorithm:

 * input:

 * output:

 * side-effects:

 */

int
parse_input (char * fn_infile,
                   int fn_len);

#endif /* PARSE_INPUT_H */
