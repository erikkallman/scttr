/* This file is part of get nums. */

/* Get nums is free software: you can redistribute it and/or modify */
/* it under the terms of the GNU Lesser General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or */
/* (at your option) any later version. */

/* get nums is distributed in the hope that it will be useful, */
/* but without any warranty; without even the implied warranty of */
/* merchantability or fitness for a particular purpose. See the */
/* GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License */
/* along with get nums, found in the "license" subdirectory of the root */
/* directory of any program using the get nums library.*/
/*   If not, see <http://www.gnu.org/licenses/>. */
#ifndef GET_NUMSL_H
#define GET_NUMSL_H
#define BUF_SIZE 256
#define BIN_FLIP(x) ((x) == 1 ? 0 : 1)

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
int
get_numinstr (char *s, char *buf, int idx, int str_len);

/* function get_nums:

get_nums is a function used to, for any given string containing unsorted
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
get_nums (char *str, int *idxs_out, int str_len, int n_idxs,...);

#endif /* GET_NUMSL_H */
