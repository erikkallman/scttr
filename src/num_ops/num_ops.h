#ifndef NUM_OPS_H
#define NUM_OPS_H

int power(int base,
          int exp
          );

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
#endif /* NUM_OPS_H */
