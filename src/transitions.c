/* Copyright (C) 2015 Erik Källman */
/* This file is part of the scttr program. */

/* scttr is free software: you can redistribute it and/or modify */
/* it under the terms of the GNU Lesser General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or */
/* (at your option) any later version. */

/* scttr is distributed in the hope that it will be useful, */
/* but without any warranty; without even the implied warranty of */
/* merchantability or fitness for a particular purpose. See the */
/* GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License */
/* along with scttr, found in the "license" subdirectory of the root */
/* directory of the scttr program. */
/* If not, see <http://www.gnu.org/licenses/>. */
/**
   * @file transitions.c
   * @author Erik Källman
   * @date November 2015
   * @brief This file contains the implementation of all functions defined for
   * the @p trs variable in the input node struct (see @p inp_node).
   */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "iquicks.h"
#include "transitions.h"
#include "std_num_ops.h"
#include "sci_const.h"
#include "inp_node_s.h"

int
get_erange (struct inp_node *inp, double e)
{
  int j;

  for (j = 1; j < inp -> md -> state_er[0]; j += 2 ) {
    if (inrange((e - inp -> e0) * AUTOEV
                , inp -> md -> state_er[j]
                , inp -> md -> state_er[j + 1])) {
        return j - 1;
    }
  }

  return -1;

}

int
get_trs (int from, double **trs)
{

  int last_i = (int)trs[0][0];
  int j = 0;

  while (last_i != -1) {
    if ((int)trs[0][j] == from) {
      return j;
    }
    j++;
    last_i = (int)trs[0][j];
  }

  return -1;
}

int
get_inext (double **trs, int from)
{

  int j = 0;

  while ((int)trs[0][j] != -1) {
    if ((int)trs[0][j] == from) {
      break;
    }
    j++;
  }

  if ((int)trs[0][j] != from) {
    return (int)trs[0][j];
  }

  while((int)trs[0][j] == from) {
    j++;
  }

  return j;
}


void
count_states (struct inp_node *inp)
{

  int j = 0; /* looping variables */
  int last_i = -2;
  inp -> n_gfs = 0;
  inp -> n_is = 0;
  inp -> n_tmax = 0;

  int s_idx;
  int t_max;

  int n_proc = 0;
  int *t;
  int *proc_is;

  double *state_er = inp -> md -> state_er;

  if((proc_is = malloc(inp -> n_trans*sizeof(int))) == NULL ) {
    fprintf(stderr, "transitions.c, count_states: malloc: failed to allocate memory for \"proc_is\"\n");
    printf("program terminating due to the previous error.\n");
    exit(EXIT_FAILURE);
  }

  if((t = malloc(inp -> n_trans * sizeof(int))) == NULL ) {
    fprintf(stderr, "transitions.c, count_states: malloc: failed to allocate memory for \"t\"\n");
    printf("program terminating due to the previous error.\n");
    exit(EXIT_FAILURE);
  }

  for (j=0; j <inp -> n_trans; j++) {
    t[j] = 0;
  }

  j = 0;

  /* gount all ground and final states */
  while((int)(inp -> trs[0][j]) != -1) {
    if(((inp -> trs[2][j]-inp -> e0)*AUTOEV >= state_er[1])
       && ((inp -> trs[2][j]-inp -> e0)*AUTOEV <= state_er[2])) {
      if (last_i != (int)inp -> trs[0][j]) {
        inp -> n_gfs++;
      }
      if(((inp -> trs[3][j]-inp -> e0)*AUTOEV >= state_er[3])
         && ((inp -> trs[3][j]-inp -> e0)*AUTOEV <= state_er[4])) {
        if ((s_idx = intinint(proc_is,(int)inp -> trs[1][j], n_proc)) == -1) {
          proc_is[n_proc] = (int)inp -> trs[1][j];
          t[n_proc] += 1;
          n_proc++;
        }
        else{
          t[s_idx] += 1;
        }
        inp -> n_is++;
      }
    }

    last_i = (int)inp -> trs[0][j++];
  }

  t_max = 0;
  for (j=0; j < n_proc; j++) {
    if (t[j] > t_max) {
      t_max = t[j];
    }
  }

  inp -> n_tmax = t_max;
  free(proc_is);
  free(t);
}

int
add_eltrans (struct inp_node *inp)
{

  printf("      adding elastic transitions ..\n");

  int j, k, l;
  int n_proc, nb;
  int tmp_idx, next_to;
  int last_i;
  int nt_el = inp -> n_trans;
  int sz_buf = inp -> n_tmax;

  /* in the worst case, the number of transitions are doubled when adding
     elastic transitions to the trs matrix */
  long int sz_el = (inp -> n_trans * 2) + 1;

  if (nt_el > (int)sz_el) {
    fprintf(stderr, "parse_input.c, function add_eltrans: input buffer writing outside its memory. nt_el = %d >= sz_el = %d.\n"
            , nt_el, (int)(sz_el));
    printf("program terminating due to the previous error.\n");
    exit(1);
  }

  double e0 = inp -> e0;
  double *state_er = inp -> md -> state_er;

  /* an array to keep track of which states that have had their
   elastic transitions accounted for */
  int *proc_st;

  double **trs_buf;

  /* transition matrix with enough space to acommodate the elastic
     transitions */
  double **trs_el;

  if((trs_buf = malloc(6 * sizeof(double *))) == NULL ) {
    fprintf(stderr, "transitions.c, function add_eltrans: failed to allocate memory for \"trs_buf\"\n");
    printf("program terminating due to the previous error.\n");
    exit(1);
  }

  for (j = 0; j < 6; j++) {
    if((trs_buf[j] = malloc(sz_buf * sizeof(double))) == NULL ) {
      fprintf(stderr, "transitions.c, function add_eltrans: failed to allocate memory for pointers in \"trs_buf[%d]\"\n"
              ,j);
      printf("program terminating due to the previous error.\n");
      exit(1);
    }
  }

  if((trs_el = malloc(6 * sizeof(double *))) == NULL ) {
    fprintf(stderr, "transitions.c, function add_eltrans: failed to allocate me mory for \"trs_el\"\n");
    printf("program terminating due to the previous error.\n");
    exit(1);
  }

  for (j=0; j < 6; j++) {
    if((trs_el[j] = malloc((sz_el + 1) * sizeof(double))) == NULL ) {
      fprintf(stderr, "transitions.c, function add_eltrans: failed to allocate memory for pointers in \"trs_el[%d]\". sz_el = %d\n"
              ,j, (int)sz_el);
      printf("program terminating due to the previous error.\n");
      exit(1);
    }
  }

  if((proc_st = malloc(inp -> n_trans * sizeof(int))) == NULL ) {
    fprintf(stderr, "transitions.c , function add_eltrans, malloc: failed to allocate memory for \"proc_st\"\n");
    printf("program terminating due to the previous error.\n");
    exit(1);
  }

  /* copy the old trans data to trs_el */
  for (l = 0; l < inp -> n_trans; l++) {
    trs_el[0][l] = inp -> trs[0][l];
    trs_el[1][l] = inp -> trs[1][l];
    trs_el[2][l] = inp -> trs[2][l];
    trs_el[3][l] = inp -> trs[3][l];
    trs_el[4][l] = inp -> trs[4][l];
    trs_el[5][l] = inp -> trs[5][l];
    next_to = (int)inp -> trs[0][l];
  }

  trs_el[0][l] = inp -> trs[0][l] = -1;

  j = nb = n_proc = 0;
  while ((int)inp -> trs[0][j] > 0) {

    printf("        %.2f%%\r", (((float)j / (float)inp -> n_trans) * 100));
    if ((intinint(proc_st, (int)inp -> trs[1][j], n_proc) == -1)
        && (((inp -> trs[2][j] - e0) * AUTOEV >= state_er[1])
            && ((inp -> trs[2][j] - e0) * AUTOEV <= state_er[2]))
        && (((inp -> trs[3][j] - e0) * AUTOEV >= state_er[3])
            && ((inp -> trs[3][j] - e0) * AUTOEV <= state_er[4]))) {

      next_to = (int)inp -> trs[1][j];

      /* found a "to" state that has not had its elastic transitions
       added yet. loop over the trs matrix and check if there are transitions
       from other states that need to be taken into account as well. */
      nb = 0;
      l = j;
      while ((int)inp -> trs[0][l] != -1) {
        if (((int)inp -> trs[1][l] == next_to)
            && (((inp -> trs[2][l] - e0) * AUTOEV >= state_er[1])
                && ((inp -> trs[2][l] - e0) * AUTOEV <= state_er[2]))
            ) {

          /* transitions from another state */
          trs_buf[0][nb] = inp -> trs[1][l];
          trs_buf[1][nb] = inp -> trs[0][l];
          trs_buf[2][nb] = inp -> trs[3][l];
          trs_buf[3][nb] = inp -> trs[2][l];
          trs_buf[4][nb] = inp -> trs[4][l];
          trs_buf[5][nb] = inp -> trs[5][l];

          nb++;

          /* jump to the next "from" state, since any given state can only have
           one transition to another specific state */
          last_i = inp -> trs[0][l];
          while((int)inp -> trs[0][l] != -1) {
            if ((int)inp -> trs[0][l] != last_i) {
              break;
            }
            l++;
          }
        }
        else {
          l++;
        }
      }
      if ((nt_el + nb + 1) > sz_el) {
      fprintf(stderr, "transitions.c, function add_eltrans: input buffer writing outside its memory. nt_el+nb+1 = %d >= sz_el = %ld.\n"
              ,nt_el+nb+1,sz_el);
        printf("program terminating due to the previous error.\n");
        exit(1);
      }

      /* if the from state cant be found in trs, just store the data in the last
         available place in trs_el */
      /* otherwise use the fwdsplice function to add it to trs_el */
      if (get_trs(next_to, trs_el) == -1) {
        for (k = 0; k < nb; nt_el++, k++) {
          trs_el[0][nt_el] = trs_buf[0][k];
          trs_el[1][nt_el] = trs_buf[1][k];
          trs_el[2][nt_el] = trs_buf[2][k];
          trs_el[3][nt_el] = trs_buf[3][k];
          trs_el[4][nt_el] = trs_buf[4][k];
          trs_el[5][nt_el] = trs_buf[5][k];
        }
        trs_el[0][nt_el + 1] = -1;
      }
      else {
        tmp_idx = get_inext(trs_el, next_to);
        fwdsplice(trs_buf, trs_el, tmp_idx, nt_el, nb, 6);
        trs_el[0][nt_el + nb + 1] = -1;
        nt_el += nb;
      }
      proc_st[n_proc++] = next_to;
    }
    j++;
  }

  trs_el[0][nt_el] = -1;

  for (j = 0; j < 6; j++) {
    free(inp -> trs[j]);
  }
  free(inp -> trs);
  inp -> trs = NULL;

  /* allocate new space for trs, now that we know the total size */
  if((inp -> trs = malloc(6*sizeof(double *))) == NULL ) {
    fprintf(stderr, "transitions.c, function add_eltrans: failed to allocate memory for \"input_data\"\n");
    printf("program terminating due to the previous error.\n");
    exit(1);
  }

  for (j = 0; j < 6; j++) {
    if((inp -> trs[j] = malloc((nt_el + 1) * sizeof(double))) == NULL ) {
      fprintf(stderr, "transitions.c, function add_eltrans: failed to allocate memory for pointe rs in \"input_data\"\n");
      printf("program terminating due to the previous error.\n");
      exit(1);
    }
  }
  for (l = 0; l < 6; l++) {
    for (j = 0; j <= nt_el; j++) {
        inp -> trs[l][j] = trs_el[l][j];
    }
  }

  inp -> trs[0][nt_el] = -1;

  for (j = 0; j < 6; j++) {
    free(trs_buf[j]);
  }
  free(trs_buf);

  for (j = 0; j < 6; j++) {
    free(trs_el[j]);
  }
  free(trs_el);
  free(proc_st);

  inp -> n_trans = nt_el;

  printf("\n      done.\n");
  return 1;
}
