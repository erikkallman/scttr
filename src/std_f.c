#include <stdio.h>
#include <stdlib.h>
#include "std_f.h"

int *
appc_d (int * a1,
        double * a2,
        int sz1,
        int sz2) {

  int j,k; /* looping variables */
  int sz3 = sz1+sz2;
  int * a3;
  /* if (a1 == NULL) { */
  /*   fprintf(stderr, "Found NULL\n"); */
  /*   printf( "program terminating due to the previous error.\n"); */
  /*   exit(1); */
  /* } */
  if((a3 = malloc((sz1+sz2)*sizeof(int))) == NULL ){
    fprintf(stderr, "std_f.c:function comb_array, malloc: failed \
to allocate memory for \"tmp_dat\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }
  /* printf( "%d %d %d\n", sz1, sz2, sz3); */
  for (j=0; j<sz1; j++) {
    a3[j] = a1[j];
    /* printf( "a1[%d] = %d\n", j, a1[j]); */
  }

  k=j;
  free(a1);

  for (j=0; j<sz2; j++) {
    a3[j+k] = (int)a2[j];
    /* printf( "a2[%d] = %d\n", j, (int)a2[j]); */
    /* printf( "jk=%d\n", j+k); */
  }
  /* printf( "\n" ); */

  return a3;
}


double *
appc (double * a1,
      double * a2,
      int sz1,
      int sz2) {

  int j,k; /* looping variables */
  int sz3 = sz1+sz2;
  double * a3;

  if((a3 = malloc((sz1+sz2)*sizeof(double))) == NULL ){
    fprintf(stderr, "std_f.c:function comb_array, malloc: failed \
to allocate memory for \"tmp_dat\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  for (j=0; j<sz1; j++) {
    a3[j] = a1[j];
  }

  k=j;
  free(a1);

  for (j=0; j<sz2; j++) {
    a3[j+k] = a2[j];
  }

  return a3;
}
