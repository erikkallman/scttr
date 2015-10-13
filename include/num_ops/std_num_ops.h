#ifndef STD_NUM_OPS_H
#define STD_NUM_OPS_H
#include <complex.h>

double
fabsc (double complex c1);

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
get_minl (double * v,
          int n);

int
power(int base,
      int exp);

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
