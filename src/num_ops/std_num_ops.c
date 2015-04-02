#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "std_num_ops.h"
#include "sci_const.h"

int
cupto (double * a1,
       double * a2,
       int upto){
  int j = 0;
  for (j=0; j<upto; j++) {
    a2[j] = a1[j];
  }
  return 0;
}

double
pyth_distl (double q,
            double p){
  return sqrt((q-p)*(q-p));
}

double
get_maxl (double * v,
          int n){
  int j;
  double m = v[0];

  for (j=0; j<n; j++)
    if (v[j] > m)
      m = v[j];

  return m;
}

double
arit_meanl (double * a,
            int n_vals){
  int j; /* looping variables */
  double sum = 0;

  for (j=0; j<n_vals; j++)
    sum += a[j];

  return sum/n_vals;
}

int
power(int base,
      int exp){
  int power;
  power=1;
  while(exp-- > 0)
    power *=base;

  return power;
}

double
powerl(double base,
       int exp){
  double power;
  power = 1;
  while(exp-- > 0)
    power *=base;

  return power;
}

int
intinint (int * a,
          int num,
          int n_el){
  int j; /* looping variables */

  for (j=0; j<n_el; j++){
    if (a[j] == num) {
      return 1;
    }
  }

  return 0;
}

int *
getintinint (int * a1,
             int * a2,
             int n){
  int j;
  int k=0;
  int * res;

  if((res = malloc(n*sizeof(int))) == NULL ){
    fprintf(stderr, "parse_input.c:function getintinint, malloc: failed to \
allocate memory for \"state_indices\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  for (j=0; j<n; j++) {
    if (intinint(a2, a1[j], n)) {
      res[k++] = a1[j];
    }
  }

  return res;
}

double
get_bdist (double e_val){
  /* return exp((-e_val/(double)AUTOEV)/((1.380648813e-23)*(double)(TEXP*TTOEV))); */
  return exp((-e_val*(double)AUTOEV)/((8.6173324e-05)*(double)(300)));
}

double
get_rbdist (double e_rel,
           double e_val){
  return exp((-(e_val-e_rel)*(double)AUTOEV)/((8.6173324e-05)*TEXP));
}
