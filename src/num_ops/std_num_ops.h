#ifndef STD_NUM_OPS_H
#define STD_NUM_OPS_H

/* function

   * synopsis:
   Check if value @v is inside of the closed range bound by @r1 and @r2
   */
int inrange(double v,
            double r1,
            double r2);

/* function get_wi

   * synopsis:
   get long type value from array, using an input array of integer indices and a requested integer

   * side-effects: none

   */
double
get_wi (double * vals,
        int * idxs_map,
        int wanted_idx,
        int n_idxs);

/* function get_wi
 * synopsis:
 splice the values from the array @from, into the array @into by shifting the
 latter forwards, n_vals elements. assumes that from is of the same
 dimensionality as into

*/

int
fwdsplice (double ** from,
           double ** into,
           int start, /* on what element in @into to start splicing */
           int end, /* current last allocated element in into*/
           int s,
           int n_dims
           );

/* function cupto

   * synopsis:
   copy the content of @a1 to @a2 upto the element @upto
   * algorithm:

   * input:

   * output:

   * side-effects:

   */
int
cupto (double * a1,
       double * a2,
       int upto);

/* function get_maxl

   * synopsis:

   * algorithm:

   * input:

   * output:

   * side-effects:

   */
double
get_maxl (double * v,
          int n);

double
get_min (double * v,
          int n);
/* function arit_meanl

   * synopsis:

   * algorithm:

   * input:

   * output:

   * side-effects:

   */
double
arit_meanl (double * a,
            int n_vals);

double
pyth_distl (double q,
            double p);

int power(int base,
          int exp
          );

double
powerl(double base,
       int exp);
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
          int n_el
          );

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
             int n
             );

double
get_bdist (double e_val);

double
get_rbdist (double e_rel,
            double e_val);

/* function lorz

   * synopsis:

   * algorithm:

   * input:

   * output:

   * side-effects:

   */
double
lorz (double x_diff,
      double fwhm);
#endif /* STD_NUM_OPS_H */
