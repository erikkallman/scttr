#ifndef PARSE_INPUT_H
#define PARSE_INPUT_H

/* function tau_parse_jobfile:

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
