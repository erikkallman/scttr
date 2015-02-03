#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "parse_input.h"
#define BUF_SIZE 256

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
                   int fn_len) {

  /* int j, k, l; /\* looping variables *\/ */
  /* FILE * sim_fp; /\* pointer to the opened simulation data file *\/ */
  /* FILE * job_fp; /\* pointer to the job file *\/ */

  /* char simf_lines[BUF_SIZE][BUF_SIZE]; */
  /* char jobf_lines[BUF_SIZE][BUF_SIZE]; */
  /* char tmp_str1[BUF_SIZE]; */
  /* char tmp_str2[BUF_SIZE]; */

  /* /\* read the simulation data file  *\/ */
  /* if((sim_fp = fopen("sim.dat", "rt")) == NULL) { */
  /*   printf("tau_parse_jobfile: can not open simulation data file. check formating\n"); */
  /* } */

  /* /\* read the job file  *\/ */
  /* if((job_fp = fopen(job_fn, "rt")) == NULL) { */
  /*   printf("tau_parse_jobfile: can not open job file. check formating\n"); */
  /* } */

  /* /\* scan the file line by line until we reach the end of the file. */
  /*    store the content. *\/ */
  /* for(j=0; fgets(simf_lines[j],BUF_SIZE,sim_fp) != NULL; j++){} */
  /* for(j=0; fgets(jobf_lines[j],BUF_SIZE,job_fp) != NULL; j++){} */

  /* fclose(sim_fp); */
  /* fclose(job_fp); */

  /* strcpy(tmp_str1,jobf_lines[0]); */
  /* /\* /\\* reset the content of job_fn *\\/ *\/ */
  /* /\* for (j=0; j<BUF_SIZE; j++) { *\/ */
  /* /\*   job_fn[j] = '\0'; *\/ */
  /* /\* } *\/ */

  /* /\* printf( "str1\n" ); *\/ */
  /* /\* for (k=0; tmp_str1[k] != '\0'; k++) { *\/ */
  /* /\*   printf( "%c", tmp_str1[k]); *\/ */
  /* /\* } *\/ */

  /* /\* loop over the simf_lines and identify the simulation method string */
  /*    identifier that was found in the job file, which is located on the first */
  /*    line of the job file *\/ */

  /* for(j=0;;j++) { */
  /*   strcpy(tmp_str2,simf_lines[j]+2); */

  /*   /\* printf( "str2\n" ); *\/ */
  /*   /\* for (k=0; tmp_str2[k] != '\0'; k++) { *\/ */
  /*   /\*   printf( "%c", tmp_str2[k]); *\/ */
  /*   /\* } *\/ */
  /*   /\* sleep(1); *\/ */

  /*   if (strcmp(tmp_str1, tmp_str2) == 0) { */

  /*     /\* store the numerical identifier of the simulation method in the first */
  /*        element of the @job_fn array .. *\/ */
  /*     job_fn[0] = simf_lines[j][0]; */

  /*     /\* .. and the simulation method input file name string in the remaining part of the */
  /*        array *\/ */
  /*     strcpy(&job_fn[1],jobf_lines[1]); */

  /*     /\* add the end of line identifier since it is lost from the string copy */
  /*      operation above *\/ */
  /*     job_fn[str_size(&job_fn[1])] = '\0'; */

  /*     break; */
  /*   } */
  /*   if (tmp_str2[0] == '\0') { */
  /*     strcpy(&job_fn[0],tmp_str1); */
  /*     return 1; */
  /*   } */
  /* } */
  return 0;
}
