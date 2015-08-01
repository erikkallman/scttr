#ifndef PARSE_INPUT_H
#define PARSE_INPUT_H
#include "dynarray.h"
#include "structs.h"

int
add_sym (double * state_er);

int
parse_input_bin (char * bin_fpstr
                 );
int
parse_molout (char * fn_relpath,
              char * fn_infile,
              int len_infile
              );

int ISINSIDE(double v,
             double r1,
             double r2);

double
get_efrom (int from);

double
get_eto (int to);

/* function check_pi

   * synopsis:

   * algorithm:

   * input:

   * output:

   * side-effects:

   */
int
check_pi ();

void
set_tmax();

void
state2s(int idx);

void
astate2s(double ** pi,
         int idx);
/* function pi2s

   * synopsis:

   * algorithm:

   * input:

   * output:

   * side-effects:

   */
void
pi2s ();

void
pi2f (char * fn);

/* function get_t

   * synopsis:

   * algorithm:

   * input:

   * output:

   * side-effects:

   */
double
get_t (int from,
       int to);
int
get_i (int from
       );

int
get_il (int from
        );

int
get_inext (int from
       );

int
get_ilnext (int from
       );

int
get_pinext (double ** pi,
           int from
            );
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

/* function get_es

   * synopsis:

   * algorithm:

   * input:

   * output:

   * side-effects:

   */
double*
get_es (int idx);

/* function parse_input:
 * synopsis:
 parse_input is the output data file interface to smap. For a given output type,
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
parse_input_tmp (double * state_er,
                 char * fn_tmpdata,
                 char * bin_fpstr
                 );

int
parse_input (double * state_er,
             char * fn_relpath, /* name of input file */
             char * tmp_fpstr,
             char * format,
             char * bin_fpstr,
             int len_fn);

#endif /* PARSE_INPUT_H */
