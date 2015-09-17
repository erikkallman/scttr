#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "std_num_ops.h"
#include "sci_const.h"

double
get_wi (double * vals,
        int * idxs_map,
        int wanted_idx,
        int n_idxs) {

  int j;
  for (j=0; j<n_idxs; j++) {
    if (idxs_map[j] == wanted_idx-1) {
      return vals[j];
    }
  }
  return -1;
}

int
fwdsplice (double ** from,
           double ** into,
           int start, /* on what element in @into to start splicing */
           int end, /* current last allocated element in into*/
           int s,
           int n_dims
           ) {

  int j,k; /* looping variables */

  if (end < start) {
    fprintf(stderr, "\n\nERROR: std_num_ops.c, function fwdsplice: there is not enough room in the @into array for the specified data in @from to fit. end = %d < start = %d, shift = %d.\n\n",end,start,s);
    printf( "program terminating due to the previous error.\n");
    exit(EXIT_FAILURE);
  }
  else if(s>(end-start)){
    end += s;
  }
  else if((start < 0) || (end < 0) || (s < 0)){
    fprintf(stderr, "\n\nERROR: std_num_ops.c, function fwdsplice: negatigve input arguments: start = %d, end = %d, s = %d .\n\n",start,end,s);
    printf( "program terminating due to the previous error.\n");
    exit(EXIT_FAILURE);
  }

  /* make space for the new values by shifting all values forward */
  /*  we know how many values are already in the into array */
  /* start++; */
  for (j=end; j>=start; j--) {
    for (k=0; k<n_dims; k++) {

      if (memcpy(&into[k][j+s],&into[k][j],sizeof(into[0][0])) != &into[k][j+s]) {
        fprintf(stderr, "parse_input.c, function sharr_fwd: the double %le stored at memory location %p cannot be copied to location %p.\n",into[j][k],&into[j][k+s],&into[j][k]);
        printf( "program terminating due to the previous error.\n");
        exit(1);
      }

      /* into[k][j] = from[k][j-start]; */
    }
  }

  for (j=0; j<s; j++) {
    for (k=0; k<n_dims; k++) {
      into[k][j+start] = from[k][j];
    }
  }


  return 1;
}

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
get_min (double * v,
          int n){
  int j;
  double m = v[0];

  for (j=0; j<n; j++)
    if (v[j] < m)
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
      return j;
    }
  }

  return -1;
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
    if (intinint(a2, a1[j], n) != -1) {
      res[k++] = a1[j];
    }
  }

  return res;
}

double
get_bdist (double e_val){
  /* return exp(((e_val+1000)/(double)AUTOEV)/((1.380648813e-23)*(double)(TEXP*TTOEV))); */
  return exp(((e_val)*(double)AUTOEV)/((8.6173324e-05)*(double)TEXP));
}

double
get_rbdist (double e_rel,
           double e_val){
  return exp((-(e_val-e_rel)*(double)AUTOEV)/((8.6173324e-05)*TEXP));
}

double
lorz (double x_diff,
      double fwhm){
  return (1/PI)*((0.5*fwhm)/((x_diff*x_diff) + (0.25*fwhm*fwhm)));
}
