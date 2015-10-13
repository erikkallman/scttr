#ifndef PARSE_INPUT_H
#define PARSE_INPUT_H
#include "spectrum_info.h"
#include "dyn_array.h"

/* function init_screen

   * synopsis:

   * algorithm:

   * input:

   * output:

   * side-effects:

   */
screen
init_screen (spec_info s);

int
get_erange (spec_info s,
            double e);
int
add_sym (spec_info s);

int
parse_input_bin (spec_info s,
                 char * bin_fpstr);
int
parse_molout (spec_info s,
              char * fn_relpath,
              char * tmp_fpstr );

/* function check_trs

   * synopsis:

   * algorithm:

   * input:

   * output:

   * side-effects:

   */
int
check_trs (spec_info s);

int
get_i (spec_info s,
       int from);
int
get_trs (int from,
        double ** trs);

int
get_il (spec_info s,
        int from);

int
get_inext (spec_info s,
           int from);
int
get_ilnext (spec_info s,
            int from);

int
get_trsnext (double ** trs,
           int from);
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
/* mdda_s * */
/* screen_states (char * fn_relpath, */
/*                double * state_t, */
/*                double * state_er */
/*                ); */


/* function parse_input:
 * synopsis:
 parse_input is the output data file interface to smap. For a given output type,
 be it from Molcas, RACAH or other electron energy calculation programs,
 an input parsing function needs to be defined. this function will need to
 define; 1. the number of states (n_states); 2.transitions (n_trans) in found in
 the input file as well as; 3.a matrix (here called "trs"), of 5 rows
 and n_trans columns, where:

 trs[0][:] : idexes of initial states for a given transition
 trs[1][:] : idexes of final states for a given transition
 trs[2][:] : energy eigenvalues, relative to the lowest energy ground
 state of the system, for the indexes in trs[0]
 trs[3][:] : -||- , for the indexes in trs[1]
 trs[4][:] : transition moment values for each transition extracted
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
parse_input_tmp (spec_info s,
                 char * fn_tmpdata,
                 char * bin_fpstr);

/* function count_states

   * synopsis:

   * algorithm:

   * input:

   * output:

   * side-effects:

   */
void
count_states (spec_info s);

int
parse_input (spec_info s);

#endif /* PARSE_INPUT_H */
