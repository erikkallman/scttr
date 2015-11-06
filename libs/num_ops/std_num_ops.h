/* This file is part of std num ops. */

/* std num ops is free software: you can redistribute it and/or modify */
/* it under the terms of the GNU Lesser General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or */
/* (at your option) any later version. */

/* std num ops is distributed in the hope that it will be useful, */
/* but without any warranty; without even the implied warranty of */
/* merchantability or fitness for a particular purpose. See the */
/* GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License */
/* along with std num ops, found in the "license" subdirectory of the root */
/* directory of any program using the std num ops library. */
/*   If not, see <http://www.gnu.org/licenses/>. */
#ifndef STD_NUM_OPS_H
#define STD_NUM_OPS_H
#include <complex.h>

double
fabsc (double complex c1);

/* function

   * synopsis:
   Check if value @v is inside of the closed range bound by @r1 and @r2
   */
int
inrange(double v, double r1, double r2);

/* function get_wi

   * synopsis:
   get long type value from array, using an input array of integer indices and a requested integer

   * side-effects: none

   */
double
get_wi (double *vals, int *idxs_map, int wanted_idx, int n_idxs);

/* function get_wi
 * synopsis:
 splice the values from the array @from, into the array @into by shifting the
 latter forwards, n_vals elements. assumes that from is of the same
 dimensionality as into

*/

int
fwdsplice (double **from, double **into, int start, int end, int s, int n_dims);

double
get_maxl (double *v, int n);

double
get_minl (double *v, int n);

int
power(int base, int exp);

double
powerl(double base, int exp);

int
intinint (int *a, int num, int n_el);

double
get_rbdist (double e_rel, double e_val);

double
lorz (double x_diff, double fwhm);
#endif /* STD_NUM_OPS_H */
