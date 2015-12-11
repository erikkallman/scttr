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
   * @file iquicks.h
   * @author Erik Källman
   * @date November 2015
   * @brief This file contains the public interface to all functions in the
   * iquicks library.
   */
#ifndef QUICKSORT_H
#define QUICKSORT_H

/**
   * @brief For a given array @p x of @p double pointers, this function sorts
   * the a given range of the array, defined from between @p frist and @p last
   * , with a recursive quick sort algorithm that also keeps track of how the
   * indices of the elements in @p x are reordered. The values are sorted in
   * ascending order.
   *
   * @param x The array of values to be sorted.
   * @param sorted_idx The indices of the values in x.
   * @param first The index defining the start of the range to be sorted.
   * @param last The index defining the end of the range to be sorted.
   * @param n_el The number of elements in @p x and @p sorted_idx.
   * @returns void.
   * @note Side effects: re-orders the elements in @p sorted_idx according to
   * the order of values in @p x.
   */
void
iquicks(double *x, int *sorted_idx, int first, int last, int n_el);

/**
   * @brief For a given array @p x of @p double pointers, this function sorts
   * the a given range of the array, defined from between @p frist and @p last
   * , with a recursive quick sort algorithm that also keeps track of how the
   * indices of the elements in @p x are reordered. The values are sorted in
   * ascending order.
   *
   * @param x The array of values to be sorted.
   * @param sorted_idx The indices of the values in x.
   * @param first The index defining the start of the range to be sorted.
   * @param last The index defining the end of the range to be sorted.
   * @param n_el The number of elements in @p x and @p sorted_idx.
   * @returns void.
   * @note Side effects: re-orders the elements in @p sorted_idx according to
   * the order of values in @p x.
   */
void
iquicks_d(double *x, double *sorted_idx, int first, int last, int n_el);

#endif /* QUICKSORT_H */
