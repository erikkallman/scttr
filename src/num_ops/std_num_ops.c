#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "std_num_ops.h"

double
pyth_distl (double q,
            double p){
  return sqrt((q-p)*(q-p));
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

int
intinint (int * a,
          int num,
          int n_el){
  int j; /* looping variables */
  for (j=0; j<n_el; j++)
    if (a[j] == num) {
      return 1;
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
get_bdist (double e_val,
           double temp){
  return exp(e_val/((1.380648813e-23)*300));
}
