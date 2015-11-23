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
   * @file std_num_ops.h
   * @author Erik Källman
   * @date November 2015
   * @brief This file contains the public interface of functions defined  in the
   * std_num_ops library.
   */
#ifndef STD_NUM_OPS_H
#define STD_NUM_OPS_H
#include <complex.h>

/**
   * @brief Calculate the absolute value of a complex number @p c1.
   *
   * @param c1 A complex number.
   * @returns The absolute value of the complex number @p c1.
   */
double
fabsc (double complex c1);

/**
   * @brief Check if the value @p v is in the closed range of [r1,r2].
   *
   * @param v The value to be checked if in [r1,r2]
   * @param r1 Lower bound of the range [r1,r2]
   * @param f2 Upper bound of the range [r1,r2]
   * @returns 1 if @p is in [r1,r2], otherwise 0.
   */
int
inrange(double v, double r1, double r2);

/**
   * @brief For a given array @p vals, and a number of indices stored in
   * @p idxs_map that maps to @p vals, this function uses @p idxs_map to
   * retrieve the corresponding value in @p vals.
   *
   * @param vals An array of @p double numbers.
   * @param idxs_map An array of indices of values in @p vals.
   * @param wanted_idx The index to be searched for in @p idxs_map
   * @param n_idxs The total number of elements in @p vals and @p idxs_map.
   * @returns If successful, the value found at the right index @p j, otherwise
   * -1.
   */
double
get_wi (double *vals, int *idxs_map, int wanted_idx, int n_idxs);

/**
   * @brief The fwdsplice() function splices an array into another by shifting
   * the arrays of the latter downward enough elements to accomodate for the
   * former.
   *
   * @param from The array to be spliced into @p into
   * @param into The array into which the @p from array will be spliced.
   * @param start The index of the element in @p into at which to start the
   * splicing.
   * @param end The index of the element in @p into at which to end the
   * splicing.
   * @param s The number of elements to shift the elements in the @p into
   * array in order to fit the @p from array.
   * @param n_dims The number of dimensions of @p from over which to perform
   * the splicing.
   */
int
fwdsplice (double **from, double **into, int start, int end, int s, int n_dims);

/**
   * @brief Get the maximum value from an array of @p double pointer type.
   *
   * @param v An array of @p double values.
   * @param n The number of elements in @p v.
   * @return @p the maximum value in @p v.
   */
double
get_maxl (double *v, int n);

/**
   * @brief Get the minimum value from an array of @p double pointer type.
   *
   * @param v An array of @p double values.
   * @param n The number of elements in @p v.
   * @return @p the minimum value in @p v.
   */
double
get_minl (double *v, int n);

/**
   * @brief using a given a base @p base of @p int type, and an exponential @p
   * exp , this function calculates the number resulting from combining them.
   *
   * @param base The base of the value.
   * @param exp The exponent of the value.
   * @returns The number resulting from combining the @p base with @p exp.
   */
int
power(int base, int exp);

/**
   * @brief using a given a base @p base of @p double type, and an exponential
   * @p exp , this function calculates the number resulting from combining them.
   *
   * @param base The base of the value.
   * @param exp The exponent of the value.
   * @returns The number resulting from combining the @p base with @p exp.
   */
double
powerl(double base, int exp);

/**
   * @brief Check if the integer @p num is found in the array @p a.
   *
   * @param a The array to be searched for @p num.
   * @param num The integer number to look for inside @p a.
   * @param n_el The number of elements in @p a.
   * @returns If @p num is found inside @p a, returns the index of @p num in
   * @p a, otherwise returns -1.
   */
int
intinint (int *a, int num, int n_el);

/**
   * @brief A function used to generate the boltzmann weight for a
   * given energy @p e_val.
   *
   * @param e_val An energy value.
   * @returns The relative boltzmann weight calculated from @p e_val.
   */
double
get_boltzw (double e_val);

/**
   * @brief A function used to calculate the lorenztian function value
   * given a specific @p x_diff and full-width half-maximum @p fwhm.
   *
   * @param x_diff An energy difference value.
   * @param fwhm The full-width half-maximum of the lorenztian function.
   * @returns The relative boltzmann weight calculated from @p e_val.
   */
double
lorz (double x_diff, double fwhm);
#endif /* STD_NUM_OPS_H */
