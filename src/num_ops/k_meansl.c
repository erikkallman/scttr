#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "std_num_ops.h"
#include "k_meansl.h"

int
k_meansl (double * a,
          int ** groups,
          int n_vals) {

  int j,k,l,m,n; /* looping variables */

  /* group indices used to write to the correct element of the groups array */
  int g_idxs[3] = {0,0,0};

  /* initialize the reference points */
  double sum,change,last_change,step; /* difference between two data points*/
  double ref_p[3] = {a[0], arit_meanl(a, n_vals), a[n_vals-1]};

  /* distnance from each reference point to each data point */
  double ** dist;

  if((dist = malloc(3*sizeof(double*))) == NULL ){
      fprintf(stderr, "k_meansl, malloc: failed to allocate memory for\
 \"dist\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  for (j=0; j<3; j++) {
    if((dist[j] = malloc(n_vals*sizeof(double))) == NULL ){
      fprintf(stderr, "k_meansl, malloc: failed to allocate memory for\
 \"dist\"\n");
      printf( "program terminating due to the previous error.\n");
      exit(1);
    }
  }

  change = 1;
  last_change = 0;
  while(fabs(change-last_change) > 0.2){
  /* while(1){ */
    /* printf( "\n\nchange=%le\nlast_change=%le\ndiff=%le\n\n",change,last_change,fabs(change-last_change)); */
    last_change = change;
    /* sleep(1); */
    /* for each reference point, calculate the distance from it to all
   data points */
    for (j=0; j<3; j++) {
      for (k=0; k<n_vals; k++) {
        dist[j][k] = pyth_distl(ref_p[j],a[k]);
      }
    }

    /* compare the distances between the reference points and add data points
     to the corresponding group */
    for (j=0; j<3; j++) {
      for (m=0; m<n_vals; m++) { /* each value we want to group */
        for (k=0; k<3; k++) {
          /* only compare to differencenses not equal to the current */
          if (k != j) {
            if (dist[j][m] < dist[k][m]) {
              n = 1;
            } else {
              n = 0;
            }
          }
        }
        if (n == 1) {
          /* we looped over all differences and found no other reference
           point k to be closer to the m:th data point than j */
          groups[j][g_idxs[j]+1] = m+1;
          g_idxs[j] += 1;
        }
      }
    }
    /* calculate new mean values */
    for (j=0; j<3; j++) {

      sum = 0;
      for (k=0; k < g_idxs[j]; k++) {
        sum += a[groups[j][k+1]-1];
      }

      step = sum/(g_idxs[j]+1);
      if (step > ref_p[j]) { /* remember that all energies are negative */
        ref_p[j] -= step/1000;
      } else {
        ref_p[j] += step/1000;
      }
      /* ref_p[j] = sum/(g_idxs[j]+1); */
    }

    change = arit_meanl(ref_p,3);
    /* printf( "\n\nchange2=%le\nlast_change2=%le\ndiff2=%le\n\n",change,last_change,fabs(change-last_change)); */
    for (j=0; j<3; j++) {
      printf( "\ng_idxs[%d] = %d\n",j,g_idxs[j]);
      printf( "ref_p[%d] = %le\n",j,ref_p[j]);
      for (k=0; k<g_idxs[j]; k++) {
      printf( "groups[%d][%d] = %d\n",j,k,groups[j][k+1]);
        printf( "e_vals[%d] = %le\n",groups[j][k+1],a[groups[j][k+1]-1]);
      }

      /* store the number of indices in each group in the 0th element of the group
         matrix so that the callee can extract this value */
      groups[j][0] = g_idxs[j];
      g_idxs[j] = 0;
      if (ref_p[j] > 0) {
        exit(EXIT_FAILURE);
      }
    }
  }

  for (j=0; j<3; j++) {
    free(dist[j]);
  }
  free(dist);
  fprintf(stderr, "\n\n=======Valgrind eject point=======\n\n");
  exit(1);
  return 0;
}
