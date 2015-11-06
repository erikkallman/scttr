#include <stdio.h>
#include "iquicks.h"

void
iquicks(double *x, int *sorted_idx, int first, int last, int n_el)
{

  int p; /* pivot element */
  int j, i;
  double temp1;
  int temp2;

  if(first < last) {
    p = first;
    i = first;
    j = last;

    while(i < j) {
      while(x[i] <= x[p] && i < last)
        i++;
      while(x[j] > x[p])
        j--;
      if(i < j) {
        temp1 = x[i];
        x[i] = x[j];
        x[j] = temp1;

        temp1 = sorted_idx[i];
        sorted_idx[i] = sorted_idx[j];
        sorted_idx[j] = temp1;
      }
    }

    temp1 = x[p];
    x[p] = x[j];
    x[j] = temp1;

    temp2 = sorted_idx[p];
    sorted_idx[p] = sorted_idx[j];
    sorted_idx[j] = temp2;

    iquicks(x, sorted_idx, first, j - 1, n_el);
    iquicks(x, sorted_idx, j + 1, last, n_el);
  }
}

void
iquicks_d(double *x, double *sorted_idx, int first, int last, int n_el)
{

  int p, j, i;
  double temp1;
  double temp2;

  if(first < last) {
    p = first;
    i = first;
    j = last;

    while(i < j) {
      while(x[i] <= x[p] && i < last)
        i++;
      while(x[j] > x[p])
        j--;
      if(i < j) {
        temp1 = x[i];
        x[i] = x[j];
        x[j] = temp1;

        temp1 = sorted_idx[i];
        sorted_idx[i] = sorted_idx[j];
        sorted_idx[j] = temp1;
      }
    }

    temp1 = x[p];
    x[p] = x[j];
    x[j] = temp1;

    temp2 = sorted_idx[p];
    sorted_idx[p] = sorted_idx[j];
    sorted_idx[j] = temp2;

    iquicks_d(x, sorted_idx, first, j - 1, n_el);
    iquicks_d(x, sorted_idx, j + 1, last, n_el);
  }
}
