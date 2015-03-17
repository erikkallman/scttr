#ifndef CHAR_OPS_H
#define CHAR_OPS_H

int
satopow(char * s,
        int len
        );

double
satof(char * s,
      int len);

/* function sci_atof

   * synopsis:
   this function converts an array of characters to a double-sized number
   , even if its in scientific notation.

   * algorithm:

   * input:
   @s - the character array that contains the string to be converted to a number

   * output:
   double - the double-sized number resulting from applying the *algorithm
   to the input string @s

   * side-effects:
   -none

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
