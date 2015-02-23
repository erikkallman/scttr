#ifndef GET_NUMSL_H
#define GET_NUMSL_H
#define BUF_SIZE 256
#define bin_flip(x) ((x) == 1 ? 0 : 1)

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
get_numinstr (char * s,
              char * buf,
              int idx,
              int str_len
              );

/* function get_numsl:

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
           int * idxs_out, /* idexes of numbers in string that we want */
           int str_len,
           int n_idxs,
           ...
           );

#endif /* GET_NUMSL_H */
