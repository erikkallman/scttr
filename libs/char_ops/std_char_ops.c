/* Copyright (C) 2015 Erik Källman */
/* This file is part of std_char_ops. */

/* std_char_ops is free software: you can redistribute it and/or modify */
/* it under the terms of the GNU Lesser General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or */
/* (at your option) any later version. */

/* std_char_ops is distributed in the hope that it will be useful, */
/* but without any warranty; without even the implied warranty of */
/* merchantability or fitness for a particular purpose. See the */
/* GNU Lesser General Public License for more details. */

/* You should have received a copy of the GNU Lesser General Public License */
/* along with std_char_ops, found in the "license" subdirectory of the root */
/* directory of any program using the std_char_ops library.*/
/*   If not, see <http://www.gnu.org/licenses/>. */
/**
   * @file std_char_ops.c
   * @author Erik Källman
   * @date November 2015
   * @brief This file contains the implementation of the functions in the
   * std_char_opt library. It cointains some handy functions for manipulating
   * strings.
   */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <stdarg.h>
#include <inttypes.h>
#include <errno.h>
#include "limits.h"
#include "std_num_ops.h"
#include "std_char_ops.h"


#define ISDDASH(x,y) ((((x) == '-') && ((y) == '-')) ? 1 : 0) /**< check if both x and y are dash characters */
#define MAX_POWERL 100000 /**< the maximally allowed power of a number */

int
env2int (const char * name, int *v)
{
  char *s = getenv(name);

  if (s == NULL) {
    return 1;
  } else {
    *v = str2int(s);
  }

  return 0;
}

int
str2int (const char *str)
{
  uintmax_t num = strtoumax(str, NULL, 10);
  if (num == UINTMAX_MAX && errno == ERANGE) {
    fprintf(stderr, "unable to convert \n");
    printf( "program terminating due to the previous error.\n");
    exit(EXIT_FAILURE);
  }
  return (int)num;
}

char *
a2str (int * a, int n, char delim)
{
  int j = 0;
  int k = 0;

  char *s = malloc((n*2+1)*sizeof(char));

  while (k < n) {
    s[j] = a[k++] + '0';
    s[j + 1] = delim;
    j += 2;
  }
  s[j-1] = '\0';

  return s;
}

char *
concs (int n_args, ...) {

  int j,k,l;
  int slen_tot = 0;

  char *tmp_str;
  char *path;
  char **s = malloc(n_args * sizeof(char *));

  int *slen = malloc(n_args * sizeof(int));

  va_list argv;
  va_start(argv, n_args);

  for (j = 0; j < n_args; j++) {
    tmp_str = va_arg(argv, char *);
    slen[j] = strlen(tmp_str);

    s[j] = malloc(slen[j]);
    s[j] = tmp_str;
    slen_tot += slen[j];
  }

  path = malloc(slen_tot + 1);
  for(l = 0, j = 0, k = 0;; j++, l++) {
    if (j == slen[k]) {
      if (++k == n_args) {
        break;
      }
      j = 0;
    }
    path[l] = s[k][j];
  }

  path[l] = '\0';
  free(s);
  free(slen);
  va_end(argv);

  return path;
}

int
strscmp (const char *str1, const char ** str2, int n_str)
{
  int j;

  for (j = 0; j < n_str; j++)
    if (strstr(str1, str2[j]))
      return j;

  return -1;
}

int
isanyalpha (char *s, int len)
{
  int j;

  char num_key[6] = {'-','.','E','+','\0'};

  for (j = 0; j < len; j++)
    if ((isalpha(s[j]) && (charinstr(num_key,s[j]) == 0))
        || (s[j] == '#')
        || (s[j] == '*'))
      return 1;

  return 0;
}

int
charinstr (char *str, char c)
{
  int j;

  for (j = 0; str[j] != '\0' ; j++)
    if (str[j] == c)
      return 1;

  return 0;
}

int
satopow(char *s, int len)
{

  int p= -1;
  int j = len - 1;
  int psign = 0;

  double pval = 0;

  /* extract the power of the number from the back of the string*/
  while (j > 0) {
    if (s[j] == '+') {
      psign = 1;
      break;
    }
    else if(s[j] == '-') {
      psign = -1;
      break;
    }
    else if(s[j] == '.') {
      psign = 1;
      p = 0;
      break;
    }
    p++;
    j--;
  }

  if (j == 0)
    return 0;

  /* check so that we're not trying to store a number that is larger than
     what can be contained in a double */
  if (p > MAX_POWERL) {
    fprintf(stderr, "std_char_ops.c, function satof: attempted to store a number of power greater than what can be held in a double (MAX_POWERL = %d < power = %d).\n"
            ,MAX_POWERL, p);
    printf( "program terminating due to the previous error.\n");
    exit(EXIT_FAILURE);
  }

  while((p--) > 0)
    pval = pval * 10 + (s[++j] - '0');

  return pval*psign;
}

double
satof(char *s, int len)
{

  int j;
  int nsign; /* sign of the number */

  double v = 0;
  double pval = satopow(s, len);

  if (s[0] == '-') {
    nsign = -1;
    j = 1;
  } else {
    nsign = 1;
    j = 0;
  }

  /* locate the decimal point, if there is one */
  while(isdigit(s[j]) && (j < len))
    v = v * 10 + (double)(s[j++] - '0');

  if (s[j++] == '.')
    /* account for all powers below the decimal */
    while(isdigit(s[j]) && (s[j] != '\0')) {
      v = v * 10 + (double)(s[j] - '0');
      pval = pval - 1;
      j++;
    }

  if (pval < 0) {
    pval = 1 / powerl(10, -pval);
  } else {
    pval = powerl(10, pval);
  }

  return ((double)nsign) * v * pval;
}

double
sci_atof(char *s)
{

  int esign = 0;
  int sign, j, exp;
  int power(int base, int exp);

  double val, pow;

  for(j = 0; isspace(s[j]); j++){}

  sign = (s[j] == '-') ? -1 : 1;

  if(s[j] == '+' || s[j] == '-')
    j++;

  for(val = 0.0; isdigit(s[j]); j++)
    val = (double)10.0 * val + (s[j] - '0');

  if(s[j] == '.')
    j++;

  for(pow = 1.0; isdigit(s[j]); j++)
    {
      val = 10.0 * val + (s[j] - '0');
      pow *= 10.0;
    }

  if(s[j] == 'e' || s[j] == 'E')
    j++;

  if(s[j] == '+' || s[j] == '-') {
      esign = s[j];
      j++;
    }

  for(exp = 0; isdigit(s[j]); j++)
    exp = 10.0 * exp + (s[j] - '0');

  if( esign == '-')
    return (double)sign * (val / pow) / (double)power(10, exp);
  else
    return (double)sign * (val / pow) * (double)power(10, exp);
}

int
isdashes (char *s, int len)
{
  char c, last_c;

  int j;

  for (j = 0; j < len; j++) {
    c = s[j];

    if (ISDDASH(c, last_c))
      return 1;/* found at least one double dash */

    last_c = c;
  }
  /* we got through the string without finding a double dash" */
  return 0;
}

int
isempty (char *s, int len)
{
  int j;

  static char key[5] = { ' ','\n','\t','-','\0' };

  for (j = 0; j < len; j++)
    if (strchr(key, s[j]) == NULL)
      return 0;/* match found, it isnt empty */

  /* we got through the string without finding a match to key.
     the string is "empty" */
  return 1;
}

int
isdigitin (char *s, int len)
{
  int j;

  for (j = 0; j < len; j++)
    if (isdigit(s[j]))
      return 1;

  return 0;
}
