#ifndef PARSE_INPUT_H
#define PARSE_INPUT_H

extern int n_states; /* number of states used in the calculation */
extern int n_trans; /* number of electronic transitions involved */
extern double * input_data[4]; /* defined in parse_input.c /*


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
get_numinstr (const char * s,
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
           int len);


/* function parse_input_molcas:

   * algorithm:

   * input:

   * output:

   * side-effects:

   */

double **
parse_input_molcas (char * fn_infile
                    );


/* function parse_input:

   * algorithm:
   read the job file defined by the argument string following the -j flag
   * input:
   @job_fn : a BUF_SIZE sized caracter array containing the name of the job file
   @fn_len : the length of the job file name

   * output: 0 if successful, 1 otherwise

   * side-effects: over-writes @job_fn and stores the string of the input file
   found in the job file, starting at the 1st element of the @job_fn array.
   element 0 is a numerical flag used to identify which simulation program
   to call with the obtained input file string identifier.

*/

int
parse_input (char * fn_infile,
                   int fn_len);

#endif /* PARSE_INPUT_H */
