#include <stdio.h>
#include <stdlib.h>
#include "std_num_ops.h"

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
  int * res;

  if((res = malloc(n*sizeof(int))) == NULL ){
    fprintf(stderr, "parse_input.c:function getintinint, malloc: failed to \
allocate memory for \"state_indices\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  for (j=0; j<n; j++) {
    if (intinint(a2, a1[j], n)) {
      res[n] = a1[j];
    }
  }

  return res;
}
