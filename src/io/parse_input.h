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
  double e_val;

   /* transition moments values for each transition found in idxs_to */
  double * t_moms;

  /* this list is doubly linked. */
  struct e_state * next;
  struct e_state * last;
};

/* function intinint

   * synopsis:

   * algorithm:

   * input:

   * output:

   * side-effects:

   */
int
intinint (int * a,
          int num,
          int n_el);

/* function getintinint
   for a given array @a1, returns all elements of @a2 not present in a1.
   * synopsis:

   * algorithm:

   * input:

   * output:

   * side-effects:
   - allocates memory space for res, which needs to be freed up by the callee

   */
int *
getintinint (int * a1,
             int * a2,
             int n);

/* function get_numinstr:
return the nth number in a string, as a string, that can be cast into any
given type using, for instance, atoi or atof. a "number" is any series of
digits possibly containing, but not ending with, a dot. any other character
in the string is considered a separator between two numbers.
* algorithm:

* input:

* output:

* side-effects:

*/
char*
get_numinstr (char * s,
              int idx,
              int str_len,
              int flag /* variable to store the last location in the string */
              );

/* function get_nums:

get_numsl is a function used to, for any given string containing unsorted
numbers, extract the numbers and store

* algorithm:
-takes a string @str containing any number of numbers, reads the first
values in idxs_in (n) and idxs_out (m), extracts the number in @str whos
index corresponds to n and stores it in @out at the index specified in
@idxs_out.

* input:
n_idxs corresponds to the number of values that are to be extracted, and thus
also the number of variadic input arguments to the function (...).
* output:

* side-effects:

*/
int
get_numsl (char * str,
           int * idxs_out,
           int str_len,
           int n_idxs,
           ...);

/* function isempty:

   * algorithm:

   * input:

   * output:

   * side-effects:

   */
int
isempty (char * s,
         int len);


/* function isdigitin:

   * algorithm:

   * input:

   * output:

   * side-effects:

   */
int
isdigitin (char * s,
           int len
           );

/* function get_state_ll

   * synopsis:

   * algorithm:

   * input:

   * output:

   * side-effects:

   */
struct e_state *
get_state_ll (double ** parsed_input,
              int flag);


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
