/* Copyright (C) 2015 Erik Källman */
/* This file is part of std_num_ops. */

/* std_num_ops is free software: you can redistribute it and/or modify */
/* it under the terms of the GNU Lesser General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or */
/* (at your option) any later version. */

/* std_num_ops is distributed in the hope that it will be useful, */
/* but without any warranty; without even the implied warranty of */
/* merchantability or fitness for a particular purpose. See the */
/* GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License */
/* along with std_num_ops, found in the "license" subdirectory of the root */
/* directory of any program using the std_num_ops library. */
/*   If not, see <http://www.gnu.org/licenses/>. */
/**
   * @file std_num_ops.c
   * @author Erik Källman
   * @date November 2015
   * @brief This file contains implementations of functions
   * used in the std_num_ops library to perform various standard numberical
   * operations.
   */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include "std_num_ops.h"
#include "sci_const.h"

double
fabsc (double complex c1)
{
  double complex c1_r = creall(c1);
  double complex c1_i = cimagl(c1);

  return sqrt(c1_r * c1_r + c1_i * c1_i);
}

int
inrange(double v, double r1, double r2)
{

  const float thrsh = 0.0001;
  double d1 = fabs(v - r1);
  double d2 = fabs(v - r1);

  if (((v > r1) || (d1 < thrsh))
    && ((v < r2) || (d2 < thrsh))) {
    return 1;
  } else {
    return 0;
  }
}

double
get_wi (double *vals, int *idxs_map, int wanted_idx, int n_idxs) {

  int j;
  for (j = 0; j < n_idxs; j++) {
    if (idxs_map[j] == wanted_idx) {
      return vals[j];
    }
  }
  return -1;
}

int
fwdsplice (double **from, double **into, int start, int end,
           int s, int n_dims)
{

  int j,k; /* looping variables */

  if (end < start) {
    fprintf(stderr, "\n\nERROR: std_num_ops.c, function fwdsplice: there is not enough room in the @into array for the specified data in @from to fit. end = %d < start = %d, shift = %d.\n\n", end, start, s);
    printf( "program terminating due to the previous error.\n");
    exit(EXIT_FAILURE);
  }
  else if(s > (end - start)) {
    end += s;
  }
  else if((start < 0) || (end < 0) || (s < 0)) {
    fprintf(stderr, "\n\nERROR: std_num_ops.c, function fwdsplice: negatigve input arguments: start = %d, end = %d, s = %d .\n\n", start, end, s);
    printf( "program terminating due to the previous error.\n");
    exit(EXIT_FAILURE);
  }

  /* make space for the new values by shifting all values forward */
  /*  we know how many values are already in the into array */
  for (j = end; j >= start; j--) {
    for (k = 0; k < n_dims; k++) {

      if (memcpy(&into[k][j + s], &into[k][j],
                 sizeof(into[0][0])) != &into[k][j + s]) {
        fprintf(stderr, "parse_input.c, function sharr_fwd: the double %le stored at memory location %p cannot be copied to location %p.\n",into[j][k], &into[j][k + s], &into[j][k]);
        printf( "program terminating due to the previous error.\n");
        exit(1);
      }
    }
  }

  for (j = 0; j < s; j++) {
    for (k = 0; k < n_dims; k++) {
      into[k][j + start] = from[k][j];
    }
  }

  return 1;
}

double
get_maxl (double *v, int n)
{
  int j;
  double m = v[0];

  for (j = 0; j < n; j++)
    if (v[j] > m)
      m = v[j];

  return m;
}

double
get_minl (double *v, int n)
{
  int j;
  double m = v[0];

  for (j = 0; j < n; j++)
    if (v[j] < m)
      m = v[j];

  return m;
}

int
power(int base, int exp)
{
  int power;
  power=1;
  while(exp-- > 0)
    power *=base;

  return power;
}

double
powerl(double base, int exp)
{
  double power;
  power = 1;
  while(exp-- > 0)
    power *= base;

  return power;
}

int
intinint (int *a, int num, int n_el)
{
  int j; /* looping variables */

  for (j = 0; j < n_el; j++) {
    if (a[j] == num) {
      return j;
    }
  }

  return -1;
}

double
get_boltzw (double e_val)
{
  return exp(-e_val / ((8.6173324e-05) * TEXP));
}

double
lorz (double x_diff, double fwhm)
{
  return (1 / PI) * ((0.5 * fwhm) / ((x_diff * x_diff) + (0.25 * fwhm * fwhm)));
}