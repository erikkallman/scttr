/* Copyright (C) 2015 Erik Källman */
/* This file is part of iquicks. */

/* iquicks is free software: you can redistribute it and/or modify */
/* it under the terms of the GNU Lesser General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or */
/* (at your option) any later version. */

/* iquicks is distributed in the hope that it will be useful, */
/* but without any warranty; without even the implied warranty of */
/* merchantability or fitness for a particular purpose. See the */
/* GNU Lesser General Public License for more details. */

/* You should have received a copy of the GNU Lesser General Public License */
/* along with iquicks, found in the "license" subdirectory of the root */
/* directory of any program using the iquicks library. */
/*   If not, see <http://www.gnu.org/licenses/>. */
/**
   * @file iquicks.c
   * @author Erik Källman
   * @date November 2015
   * @brief This file contains the implementation of functions in the iquicks
   * library.
   */
#include <stdio.h>
#include <stdlib.h>
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
  if ((n_el > 0) && ((first > n_el) || (last > n_el))) {
    fprintf(stderr, "\n\n file iquicks.c, function iquicks_d: arguments are not in an acceptable range: first = %d, last = %d, n_el = %d\n", first, last, n_el);
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

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
