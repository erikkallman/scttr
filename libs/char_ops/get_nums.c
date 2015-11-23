/* Copyright (C) 2015 Erik Källman */
/* This file is part of get_nums. */

/* get_nums is free software: you can redistribute it and/or modify */
/* it under the terms of the GNU Lesser General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or */
/* (at your option) any later version. */

/* get_nums is distributed in the hope that it will be useful, */
/* but without any warranty; without even the implied warranty of */
/* merchantability or fitness for a particular purpose. See the */
/* GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License */
/* along with get_nums, found in the "license" subdirectory of the root */
/* directory of any program using the get_nums library.*/
/*   If not, see <http://www.gnu.org/licenses/>. */
/**
   * @file get_nums.c
   * @author Erik Källman
   * @date November 2015
   * @brief This file contains implementations of functions
   * used to extract numbers from strings.
   */
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include "std_char_ops.h"
#include "get_nums.h"

int
get_numinstr (char *s, char *buf, int idx, int str_len)
{

  char c = NULL;
  char last_c = NULL;

  int j, k, l;
  int mode = 0;  /* mode = 1 => reading lines of only digits,
                    mode = 0 to reading anything else */

  char num_key[5] = {'-','.','E','+','\0'};

  l = 0;
  k = 0;

  for (j = 0; j <= str_len; c = s[j], j++) {
    /* check if we're still reading a number and avoid dashed lines */
    if (((charinstr(num_key,c) != 0) || (isdigit(c) != 0)) \
        && ((isalpha(last_c) == 0 ) || (last_c == 'E'))) {

      /* check if we're reading the right digit, else ignore the result */
      if (k == idx) {
        buf[l++] = c;
      }

      if (mode == 0) {
        /* now proceeding to read a digit */
        mode = BIN_FLIP(mode);
      }
    }
    /* if we read a non-number character and we're in read mode, we have
     read the end of a digit. */
    else if(mode == 1) {
      if (++k > idx) { /* increase the counter for read numbers */
        break;
      }
      mode = BIN_FLIP(mode); /* go back to reading non numbers */
    }
    last_c = c;
  }
  buf[l] = '\0';
  return l;
}

int
get_nums (char *str, int *idxs_out, int str_len, int n_idxs, ...)
{
  int j, k, l;
  va_list argv;

  char *num_buf;
  double *tmp_num;

  /* store the memory address to the variadic input arguments so that we can
     assign them locally */
  va_start(argv, n_idxs);

  /* a number in the string can maximally be as long as the string */
  num_buf = malloc(str_len);

  for (j = 0; j < n_idxs; j++) { /* loop over indexes */

    tmp_num = va_arg(argv, double*); /* grab the next vararg */

    /* find the correctly indexed number in the string and store it. */
    l = get_numinstr(str, num_buf, idxs_out[j], str_len);

    for (k = 0; k < l; k++) {}
    *tmp_num = satof(num_buf, k + 1); /* extract the next memory location */
  }

  free(num_buf);

  va_end(argv);
  return 1;
}
