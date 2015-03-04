#include <stdio.h>
#include <stdlib.h>
#include "smap.h"

double **
calc_smap (char * method,
           double ** trans_data
           ) {



  double ** test = malloc(sizeof(double*));
  test[0] = malloc(sizeof(double));

  free(trans_data[0]);
  free(trans_data);
  printf( "calc_smap got this method: %s \n", method);
  return test;
}
