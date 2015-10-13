#include <stdio.h>
#include "quicksort.h"

void
quicksort(double * x, int * sorted_idx,
          int first, int last, int n_el){


  int pivot,j,i;
  double temp1;
  int temp2;

  if(first<last){
    pivot=first;
    i=first;
    j=last;

    while(i<j){
      while(x[i]<=x[pivot]&&i<last)
        i++;
      while(x[j]>x[pivot])
        j--;
      if(i<j){
        temp1=x[i];
        x[i]=x[j];
        x[j]=temp1;

        temp1 = sorted_idx[i];
        sorted_idx[i] = sorted_idx[j];
        sorted_idx[j] = temp1;
      }
    }

    temp1=x[pivot];
    x[pivot]=x[j];
    x[j]=temp1;

    temp2 = sorted_idx[pivot];
    sorted_idx[pivot] = sorted_idx[j];
    sorted_idx[j] = temp2;

    quicksort(x, sorted_idx, first, j-1, n_el);
    quicksort(x, sorted_idx, j+1, last, n_el);
  }
}

void
quicksort_d(double * x, double * sorted_idx,
          int first, int last, int n_el){


  int pivot,j,i;
  double temp1;
  double temp2;

  if(first<last){
    pivot=first;
    i=first;
    j=last;

    while(i<j){
      while(x[i]<=x[pivot]&&i<last)
        i++;
      while(x[j]>x[pivot])
        j--;
      if(i<j){
        temp1=x[i];
        x[i]=x[j];
        x[j]=temp1;

        temp1 = sorted_idx[i];
        sorted_idx[i] = sorted_idx[j];
        sorted_idx[j] = temp1;
      }
    }

    temp1=x[pivot];
    x[pivot]=x[j];
    x[j]=temp1;

    temp2 = sorted_idx[pivot];
    sorted_idx[pivot] = sorted_idx[j];
    sorted_idx[j] = temp2;

    quicksort_d(x, sorted_idx, first, j-1, n_el);
    quicksort_d(x, sorted_idx, j+1, last, n_el);
  }
}
