/* This file is part of std char ops. */

/* std char ops is free software: you can redistribute it and/or modify */
/* it under the terms of the GNU Lesser General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or */
/* (at your option) any later version. */

/* std char ops is distributed in the hope that it will be useful, */
/* but without any warranty; without even the implied warranty of */
/* merchantability or fitness for a particular purpose. See the */
/* GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License */
/* along with std char ops, found in the "license" subdirectory of the root */
/* directory of any program using the std char ops library.*/
/*   If not, see <http://www.gnu.org/licenses/>. */
#ifndef CHAR_OPS_H
#define CHAR_OPS_H

/* function concs

   * synopsis:

   * algorithm:

   * input:

   * output:

   * side-effects:

   */
char *
concs (int n_args, ...);

int
strscmp (const char *str1, const char **str2, int n_str);

/* function isanychar

   * synopsis:

   * algorithm:

   * input:

   * output:

   * side-effects:

   */
int
isanyalpha (char *s, int len);

/* function charinstr

   * synopsis:

   * algorithm:

   * input:

   * output:
   returns the char in str equal to c if found, NULL otherwise.
   * side-effects:

   */
int
charinstr (char *str, char c);

int
satopow(char *s, int len);

double
satof(char *s, int len);

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
sci_atof(char *s);

/* function isdashes:

   * synopsis:

   * algorithm:

   * input:

   * output:

   * side-effects:

   */

int
isdashes (char *s, int len);

/* function isempty:

   * algorithm:

   * input:

   * output:

   * side-effects:

   */
int
isempty (char *s, int len);

/* function isdigitin:

   * algorithm:

   * input:

   * output:

   * side-effects:

   */
int
isdigitin (char *s, int len);

#endif /* CHAR_OPS_H */
