/* Copyright (C) 2015 Erik Källman */
/* This file is part of dyn_array, the dynamic array library.*/

/* dyn_array is free software: you can redistribute it and/or modify */
/* it under the terms of the GNU Lesser General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or */
/* (at your option) any later version. */

/* dyn_array is distributed in the hope that it will be useful, */
/* but without any warranty; without even the implied warranty of */
/* merchantability or fitness for a particular purpose. See the */
/* GNU Lesser General Public License for more details. */

/* You should have received a copy of the GNU Lesser General Public License */
/* along with dyn_array, found in the "license" subdirectory of the root */
/* directory of any program using the dyn_array library.*/
/* If not, see <http://www.gnu.org/licenses/>. */
/**
   * @file dyn_array.h
   * @author Erik Källman
   * @date November 2015
   * @brief This file contains the public interface to all functions in the
   * dyn_array library.
   */
#ifndef DYN_ARRAY_H
#define DYN_ARRAY_H
#define DEFAULT_CAP 10 /**< The default capacity for a newly initialized dynamic array. */
#define DEFAULT_INC 10 /**< The default increment in number of elements. */

struct da_s;
typedef struct da_s *da;
/**
   * @brief The da_s struct is fully abstracted (all manipulations of the struct
   * is done through the functions defined in the dyn_array library), so it has
   * a separate type defined for it: the @p da type. The struct only contains
   * some variables to handle the dynamic growth of the arrays, should
   * adding a value to it not be possible since it would incur writing to
   * memory that is not allocated.
   *
   * Currently, array growth and shrinking is handled through element copying
   * and array re-allocation. This is a primitive and naive implementation,
   * but is easily improvable.
   */
struct da_s {
  int n_el; /**< Number of elements containing data in the array. */
  int cap; /**< Data capacity in number of elements. */
  int inc; /**< Increment upon extension. */
  int *a; /**< The array containing the data. */
};

/**
   * @brief Whenever da_set() is called with a filled dynamic array, this
   * function extends that array before the next value is added to it.
   *
   * @param src The dynamic array to be extended.
   * @returns 1 if successful.
   * @note side effect: increases the capacity of @p src according to how many
   * elements @p src was extended.
   */
int
da_extend (da src);

/**
   * @brief This function shrinks the size of the dynamic array @p src by half.
   *
   * @param src The dynamic array to be extended.
   * @returns 1 if successful.
   * @note Side effect: decreases the capacity of @p src, and redefines @p inc
   (see the documentation for the da_s struct).
   * @returns 1 if successful.
   */
int
da_shrink (da src);

/**
   * @brief Delete the value at the element @p el in @p ar, by shifting all
   * values below @p el upwards by one step.
   *
   * @param ar The dynamic array from which the element at @p el is to be
   * deleted.
   * @param el The index of the element to be deleted.
   * @returns 1 if successful.
   */
int
da_del_us (da ar, int el);

/**
   * @brief Sets the value of element @p el in the dynamic array @p ar to @p val
   *
   * @param ar The array in which to set the value.
   * @param el Index of the element to be set.
   * @param val Integer value to be used when setting the @p el element of @p ar
   * @returns 1 if successful.
   */
int
da_set (da ar, int el, int val);

/**
   * @brief A function used to append a value @p val to a given dynamic array
   * @p ar.
   *
   * @param ar The array to @p val will be appended
   * @param val The value to be appended to @p ar.
   * @returns 1 if successful.
   * @note Side effect: if trying to append outside the range of @p ar, the
   * da_extend() function is called prior to appending the value @p val
   */
int
da_append (da ar, int val);

/**
   * @brief Get the value in the array @p a, at element @p el.
   *
   * @param a The array from which to get the value at element @p el.
   * @param el The element in @p ar at which to get a value.
   * @returns The value in @p ar, at element @p el
   * @note Exits with @p EXIT_FAILURE if trying to get a value outside the
   * range of the array.
   */
int
da_get (da a, int el);

/**
   * @brief Get the last value in the array @p a4.
   *
   * @param ar The array from which to get the last value.
   * @returns The value in @p ar, at the last element.
   * @note Exits with @p EXIT_FAILURE if trying to get a value from an empty
   * array.
   */
int
da_getlast (da ar);

/**
   * @brief This function initializes a dynamic array.
   *
   * @param cap The initial value of the @p cap variable.
   * @param inc The initial value of the @p inc variable.
   * @returns @p new_da, the newly initialized dynamic array.
   * @note Side effects: exits with EXIT_FAILURE upon failing malloc call.
   */
da
da_init (int cap, int inc);

#endif /* DYN_ARRAY_H */
