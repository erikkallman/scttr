#ifndef STD_NUM_OPS_H
#define STD_NUM_OPS_H

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
#endif /* STD_NUM_OPS_H */
