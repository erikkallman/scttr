#ifndef CHAR_OPS_H
#define CHAR_OPS_H

/* function sci_atof

   * synopsis:

   * algorithm:

   * input:

   * output:

   * side-effects:

   */
double
sci_atof(char * s);

/* function isdashes:

   * synopsis:

   * algorithm:

   * input:

   * output:

   * side-effects:

   */

int
isdashes (char * s,
          int len
          );

/* function isempty:

   * algorithm:

   * input:

   * output:

   * side-effects:

   */
int
isempty (char * s,
         int len
         );

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

#endif /* CHAR_OPS_H */
