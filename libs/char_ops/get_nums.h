/* Copyright (C) 2015 Erik Källman */
/* This file is part of get_nums. */

/* get_nums is free software: you can redistribute it and/or modify */
/* it under the terms of the GNU Lesser General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or */
/* (at your option) any later version. */

/* get_nums is distributed in the hope that it will be useful, */
/* but without any warranty; without even the implied warranty of */
/* merchantability or fitness for a particular purpose. See the */
/* GNU Lesser General Public License for more details. */

/* You should have received a copy of the GNU Lesser General Public License */
/* along with get_nums, found in the "license" subdirectory of the root */
/* directory of any program using the get_nums library.*/
/*   If not, see <http://www.gnu.org/licenses/>. */
/**
   * @file get_nums.h
   * @author Erik Källman
   * @date November 2015
   * @brief This file contains public interface to the get_nums library.
   */
#ifndef GET_NUMSL_H
#define GET_NUMSL_H
#define BUF_SIZE 256
#define BIN_FLIP(x) ((x) == 1 ? 0 : 1) /**< flip the value of x. */
/**
   * @brief This function extracts the nth number in a string. The string
   * is stored in a buffer, that can be cast into any given type using,
   * for instance, the standard UNIX functions atoi or atof. A "number" is any
   * series of digits possibly containing, but not ending with, a dot. any other
   * character in the string is considered a separator between two numbers.
   *
   * @param s the string in which to search for a number
   * @param buf the buffer in which to store the located number
   * @param idx the integer defining which of n numbers in the string to
   * store in the @p buf variable
   * @param str_len the length (in characters) of the @p s variable
   * @returns l the length (in digits) of the extracted number
   */
int
get_numinstr (char *s, char *buf, int idx, int str_len);

/**
   * @brief get_nums() is a function used to, for any given string containing
   * numbers, extract the numbers and use them to define the values of
   * whatever variables the user provides as input.
   *
   * The user can provide any number of variadic arguments to the function
   * as long as the number of variables provided as arguments matches
   * the length (in elements) of @p idxs_out, the number of numbers in the
   * @p str variable, and the value of the @p n_idxs variable.
   *
   * @param str the string containing the numbers to be extracted.
   * @param idxs_out indices of the numbers found in the @p str string.
   * @param str_len length (in characters) of the @p str variable.
   * @param n_idxs the number of indices to be extracted from @p str.
   * @param ... the @p n_idxs number of variables to be defined from the numbers
   * found in the @p str variable.
   * @returns 1 if successful.
   * @note the variadic arguments are expected to be pointers to memory
   * addresses. example call for the definition of two variables @p v1 and
   * @p v2:  get_nums(str_buffer, idxs_buffer, len_buff, n_idxs, &v1, &v2);
   */
int
get_nums (char *str, int *idxs_out, int str_len, int n_idxs,...);

#endif /* GET_NUMSL_H */
