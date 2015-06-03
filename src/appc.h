#ifndef STD_F_H
#define STD_F_H

int *
appc_dd (int * a1,
         int * a2,
         int sz1,
         int sz2);

int *
appc_d (int * a1,
        double * a2,
        int sz1,
        int sz2);

/* function appc
   append array a2 to a1 through combining them into a new array a3
   * synopsis:

   * algorithm:

   * input:

   * output:
   a3: pointer to the memory address of the newly allocated array
   * side-effects:

   */
double *
appc (double * a1,
      double * a2,
      int sz1,
      int sz2);



#endif /* STD_F_H */
