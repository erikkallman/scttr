/* This file is part of dyn array, the dynamic array library.*/

/* Dyn array is free software: you can redistribute it and/or modify */
/* it under the terms of the GNU Lesser General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or */
/* (at your option) any later version. */

/* dyn array is distributed in the hope that it will be useful, */
/* but without any warranty; without even the implied warranty of */
/* merchantability or fitness for a particular purpose. See the */
/* GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License */
/* along with dyn array, found in the "license" subdirectory of the root */
/* directory of any program using the dyn array library.*/
/* If not, see <http://www.gnu.org/licenses/>. */
#ifndef DYN_ARRAY_H
#define DYN_ARRAY_H
#define DEFAULT_CAP 10
#define DEFAULT_INC 10

struct da_s;
typedef struct da_s *da;

struct da_s {
  int n_el; /* number of elements */
  int cap; /* data capacity in number of elements */
  int inc; /* increment upon extension */
  int *a;
};

int
da_extend (da src);

int
da_shrink (da src);

int
da_del_us (da ar, int el);

int
da_set (da ar, int el, int val);

int
da_append (da ar, int val);

int
da_get (da a, int el);

int
da_getlast (da ar);

da
da_init (int cap, int inc);

#endif /* DYN_ARRAY_H */
