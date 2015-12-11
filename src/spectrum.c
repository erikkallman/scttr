/* Copyright (C) 2015 Erik Källman */
/* This file is part of the scttr program. */

/* scttr is free software: you can redistribute it and/or modify */
/* it under the terms of the GNU Lesser General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or */
/* (at your option) any later version. */

/* scttr is distributed in the hope that it will be useful, */
/* but without any warranty; without even the implied warranty of */
/* merchantability or fitness for a particular purpose. See the */
/* GNU Lesser General Public License for more details. */

/* You should have received a copy of the GNU Lesser General Public License */
/* along with scttr, found in the "license" subdirectory of the root */
/* directory of the scttr program. */
/* If not, see <http://www.gnu.org/licenses/>. */
/**
   * @file spectrum.c
   * @author Erik Källman
   * @date November 2015
   * @brief This file contains all functions defined for the spectrum struct.
   */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "spectrum.h"
#include "sci_const.h"
#include "spectrum_s.h"
#include "inp_node_s.h"
#include "transitions.h"
#include "std_num_ops.h"
#include "iquicks.h"

struct spectrum *
get_spec (struct inp_node *inp, int idx)
{
  struct spectrum *spec = inp -> root_spec;

  while ((spec -> idx ) != idx) {
    spec = spec -> next_spec;
    if (spec == NULL) {
      fprintf(stderr, "spectrum.c, function get_spec: spectrum of index %d not found in list.\n"
              ,idx);
      printf("program terminating due to the previous error.\n");
      exit(EXIT_FAILURE);
    }
  }

  return spec;
}

struct spectrum *
init_spec (struct inp_node *inp, int cap, int inc)
{

  struct spectrum *spec = malloc(sizeof(struct spectrum));

  spec -> emin_x = (inp -> md -> state_er[4] - 2) / AUTOEV;
  spec -> emax_x = (inp -> md -> state_er[3] + 2) / AUTOEV;
  spec -> emin_y = (inp -> md -> state_er[6] - 2) / AUTOEV;
  spec -> emax_y = (inp -> md -> state_er[5] + 2) / AUTOEV;

  spec -> is2fs = da_init(cap, inc);
  spec -> is_idxs = da_init(cap, inc);
  spec -> gs2is = da_init(cap, inc);
  spec -> ii_start = da_init(cap, inc);
  spec -> n_st = 0;
  spec -> s_mat = NULL;
  spec -> next_spec = NULL;
  spec -> last_spec = NULL;
  spec -> idx = 0;
  spec -> sfac = -0.1;

  return spec;
}

int
set_root_spec (struct inp_node *inp)
{

  int j, k, l;

  int r1, r2;
  int is_num = 0;
  int is_idx = 0;
  int n_is_tmp;
  int nt = 0; /* number of transitions fround in the provided energy range */
  int wh = 0; /* write head for the gs2is array */

  struct spectrum *root_spec = init_spec(inp, 10, 10);

  float bw_thrsh = inp -> md -> state_t[2];
  double bw;
  double e_gs, e_is, e_fs;

  /* keep check on what states that have been processed to make sure that no
     "to" state appears in the pi matrix more than once  */
  int *proc_st;
  int n_proc = 0;
  int is_proc;

  inp -> n_states = 0;

  if((proc_st = malloc(inp -> n_trans * sizeof(double))) == NULL ) {
    fprintf(stderr, "spectrum.c, function set_root_spec: failed to allocate memory for \"proc_st\"\n");
    printf("program terminating due to the previous error.\n");
    exit(1);
  }

  if((inp -> idx_map = malloc((inp -> n_trans + 1)
                              *  sizeof(int))) == NULL ) {
    fprintf(stderr, "spectrum.c, function set_root_spec: failed to allocate memory for \"inp->idx_map\"\n");
    printf("program terminating due to the previous error.\n");
    exit(1);
  }

  for (j = 0; j < inp -> n_trans; j++) {
    inp -> idx_map[j] = -1;
  }

  inp -> idx_map[inp -> n_trans] = -2; /* mark the end of the list */
  inp -> idx_map[0] = 0;

  printf("      input matrix integrity check and initial screening ..");
  fflush(stdout);
  /* check so that every state in PI can be reached with the get_i function */
  j = 0;

  inp -> tmax_d = inp -> tmax_q = -1;

  while((int)inp->trs[0][j] != -1) {

    /* make sure the transition is not taking place between states
       in the same energy range interval */
    r1 = get_erange(inp, inp -> trs[2][j]);
    r2 = get_erange(inp, inp -> trs[3][j]);

    if ((r1 != -1) && (r2 != -1)
      && (get_erange(inp, inp -> trs[2][j])
            != get_erange(inp, inp -> trs[3][j]))) {

      if ((int)inp -> trs[5][j] == 1) {
        if (inp -> trs[4][j] > inp -> tmax_d) {
          inp -> tmax_d = inp -> trs[4][j];
        }
      }
      else{
        if (inp -> trs[4][j] > inp -> tmax_q) {
          inp -> tmax_q = inp -> trs[4][j];
        }
      }
    }

    e_fs = inp -> trs[3][j];
    e_is = inp -> trs[2][j];

    /* is this a transition from an intermediate to final state? */
    if (inrange((e_fs-inp -> e0) * AUTOEV, inp -> md -> state_er[5]
                , inp -> md -> state_er[6])
        && inrange((e_is-inp -> e0) * AUTOEV, inp -> md -> state_er[3]
                   , inp -> md -> state_er[4])) {

      da_append(root_spec -> is2fs, j);
      is_num = inp -> trs[0][j];

      /* have the gs-is transitions for this particular is already
         been processed? */
      is_proc = intinint(proc_st, is_num, n_proc);

      if (is_proc == -1) {

        /* loop over the list of transitions looking for corresponding
           g.state->i-state transitions */
        k = 0;
        n_is_tmp = root_spec -> gs2is -> n_el;
        while((int)inp -> trs[0][k] != -1) {
          e_gs = inp -> trs[2][k];
          bw = get_boltzw((e_gs - inp -> e0) * (double)AUTOEV);
          if ((((int)inp -> trs[1][k]) == is_num)
              && inrange((e_gs - inp -> e0) * AUTOEV
                         , inp -> md -> state_er[1]
                         , inp -> md -> state_er[2])
              && (bw > bw_thrsh)
              ) {

            da_append(root_spec -> is_idxs, root_spec -> gs2is -> n_el);
            da_append(root_spec -> gs2is, k);
            nt++;
          }
          k++;
        }
        if (n_is_tmp == root_spec -> gs2is -> n_el) {

          /* no gs2is transitions found for that is2fs transition */
          root_spec -> is2fs -> n_el--;
        }
        else{

          /* gs2is transitions were added */
          da_append(root_spec -> is_idxs, -1);
          da_append(root_spec -> ii_start, wh);

          wh = root_spec -> is_idxs -> n_el;
          proc_st[n_proc++] = is_num;
        }
      }
      else {

        /* find the intermediate state in the gs2trs matrix */
        for (l=0; l < root_spec -> is_idxs -> n_el; l++) {
          is_idx = root_spec -> is_idxs -> a[l];
          if (is_idx != -1) {
            if ((int)inp -> trs[1][root_spec -> gs2is -> a[is_idx]]
                == is_num) {
              break;
            }
          }
        }
        da_append(root_spec -> ii_start, l);
        while(root_spec -> is_idxs -> a[l++] != -1)
          nt++;
      }
    }
    j++;
  }

  root_spec -> n_st = nt;
  add_spec(inp, root_spec);
  inp -> n_states = n_proc;

  if ((root_spec -> is2fs -> n_el == 0)
      || (root_spec -> gs2is -> n_el == 0)){
    fprintf(stderr, "spectrum.c, function set_root_spec: no intermediate or final states states were found in the energy range you provided. (spec -> is2fs -> n_el = %d, root_spec -> gs2is -> n_el = %d)\n"
            ,root_spec -> is2fs -> n_el,root_spec -> gs2is -> n_el);
    printf("program terminating due to the previous error.\n");
    exit(EXIT_FAILURE);
  }

  printf(" done.\n");
  free(proc_st);

  return 0;
}

int
set_spec (struct inp_node *inp)
{

  int j, k, l, m;

  int n_fs, n_is;
  int wh;
  int is_proc;
  int is_pos = 0;

  int is_idx = 0;
  int n_is_tmp = 0;
  int nt = inp -> root_spec -> n_st;

  struct spectrum *root_spec = inp -> root_spec;
  struct spectrum *spec = init_spec(inp, 10, 10);

  int n_sfs = root_spec -> is2fs -> n_el;

  int * r_gi = root_spec -> gs2is -> a;
  int * r_fi = root_spec -> is2fs -> a;
  int * r_ii = root_spec -> is_idxs -> a;

  double int_thrsh = inp -> md -> state_t[1];

  double tmom_gi,tmom_if;
  double bw;
  double int_tot,int_scr;

  double **int_dist  = malloc(sizeof(double *) * 2);
  double **tmp_evals = malloc(sizeof(double *) * 2);

  for (j = 0; j < 2; j++) {
    if((int_dist[j] = malloc(nt * sizeof(double))) == NULL ) {
      fprintf(stderr, "spectrum.c, function set_spec: failed to allocate memory for \"int_dist[%d]\"\n"
              ,j);
      printf("program terminating due to the previous error.\n");
      exit(1);
    }
    if((tmp_evals[j] = malloc(root_spec -> is2fs -> n_el * sizeof(double)))
       == NULL ) {
      fprintf(stderr, "spectrum.c, function set_spec: failed to allocate memory for \"tmp_evals[%d]\"\n"
              ,j);
      printf("program terminating due to the previous error.\n");
      exit(1);
    }
  }

  for (int_tot = 0, l = 0, j = 0; j<n_sfs; j++) {

    is_pos = root_spec -> ii_start -> a[j];
    is_idx = r_ii[is_pos];
    tmom_if = inp -> trs[4][r_fi[j]];

    for (k = is_pos; ((is_idx = r_ii[k]) != -1); k++) {

      tmom_gi = inp -> trs[4][r_gi[is_idx]];
      bw = get_boltzw((inp -> trs[2][r_gi[is_idx]] - inp -> e0) * (double)AUTOEV);
      int_dist[0][l] = l;

      int_dist[1][l] = tmom_if * tmom_gi * bw;
      int_tot += int_dist[1][l];
      l++;
    }
  }

  /* sort ascendingly according to energy */
  iquicks_d(int_dist[1], int_dist[0], 0, nt - 1, nt);

  /* screen until retaining x% of the intensity of all transitions,
   set all "sceened" states to the n_trans+1 */
  int_scr = 0;

  for (j = 0; j < nt; j++) {
    if ((int_scr / int_tot) <= int_thrsh) {
      int_scr += int_dist[1][j];
    }
    else {
      j--;
      int_dist[0][j] = -1;
      break;
    }
  }
  iquicks_d(int_dist[0], int_dist[1], j, nt - 1, nt);

  /* finally, generate the new screen by only including those states that were not screened out above */
  n_fs = n_is = 0;
  l = 0;
  m = j + 1;
  wh = 0;

  for (j = 0; j < n_sfs; j++) {

    is_pos = root_spec -> ii_start -> a[j]; /* position at is_idx where gs2is
                                     idx will start */
    is_idx = r_ii[is_pos];
    n_is_tmp = spec -> is_idxs -> n_el;

    for (k = is_pos; ((is_idx = root_spec->is_idxs->a[k]) != -1); k++) {
      if (l == (int)int_dist[0][m]) {

        /* is there an index in spec that already maps to
           this gs2is transition?*/
        if ((is_proc = intinint(spec -> gs2is -> a, root_spec -> gs2is
                                -> a[is_idx],spec -> gs2is -> n_el)) != -1) {
          da_append(spec -> is_idxs, is_proc);
        }
        else {
          da_append(spec -> is_idxs, spec -> gs2is -> n_el);
          da_append(spec -> gs2is, r_gi[is_idx]);
        }
        m++; /* start looking for the next state in int_dist */
        if (m == nt) {
          break;
        }
      }
      l++;
    }
    /* check if any gs2is transitions were added */
    if (n_is_tmp != spec -> is_idxs -> n_el) {
      da_append(spec -> is_idxs, -1);
      da_append(spec -> is2fs, r_fi[j]);
      da_append(spec -> ii_start, wh);

      wh = spec -> is_idxs -> n_el;
      tmp_evals[0][n_is++] = (inp -> trs[3]
                              [r_gi[r_ii[is_pos]]] - inp -> e0) * AUTOEV;
      tmp_evals[1][n_fs++] = (inp -> trs[3][r_fi[j]] - inp -> e0) * AUTOEV;
    }
    if (m == nt) {
      break;
    }
  }

  spec -> emin_x = floor((get_minl(tmp_evals[0], n_is) - 1)) / AUTOEV;
  spec -> emax_x = ceil((get_maxl(tmp_evals[0], n_is) + 1)) / AUTOEV;
  spec -> emin_y = floor((get_minl(tmp_evals[1], n_fs) - 1)) / AUTOEV;
  spec -> emax_y = ceil((get_maxl(tmp_evals[1], n_fs) + 1)) / AUTOEV;

  add_spec(inp, spec);

  free(tmp_evals[0]);
  free(tmp_evals[1]);
  free(tmp_evals);

  return 1;
}

int
add_spec (struct inp_node *inp, struct spectrum *spec)
{

  struct spectrum *tmp_spec;

  if ((int)sizeof(spec) != (int)sizeof(struct spectrum *)) {
    fprintf(stderr, "\n\nspectrum.c, function add_spec: variable struct spectrum * provided by callee is not initialized (sizeof(spec) = %d != %d = sizeof(struct struct spectrum)).\n"
            , (int)sizeof(spec), (int)sizeof(struct spectrum));
    printf("program terminating due to the previous error.\n");
    exit(EXIT_FAILURE);
  }

  spec -> idx = 1;
  if (inp -> root_spec == NULL) {
    inp -> root_spec = spec;
  }
  else {
    tmp_spec = inp -> root_spec;

    spec -> idx++;
    while (tmp_spec -> next_spec != NULL) {
      (spec -> idx)++;
      tmp_spec = tmp_spec -> next_spec;
    }

    spec -> last_spec = tmp_spec -> last_spec;
    tmp_spec -> next_spec = spec;
  }

  return 1;
}

int
set_trs_red (struct inp_node *inp, int spec_idx)
{

  struct spectrum *spec = get_spec(inp, spec_idx);

  if (spec -> is2fs -> n_el == 0) {
    fprintf(stderr, "spectrum.c, function set_trs_red: the provided spectrum has not had its spectrum indexing arrays defined. check your function call sequence and read the documentation on how to use the functions set_spec() and set_trs_red().\n");
    printf( "program terminating due to the previous error.\n");
    exit(EXIT_FAILURE);
  }

  int j, k, l;

  int is_idx, is_pos;
  int n_el = 0;
  /* create local copies and pointers to variables whos name would
     otherwise clutter the code loop */
  int n_sfs = spec -> is2fs -> n_el;
  int *is2fs = spec -> is2fs -> a;
  int *gs2is = spec -> gs2is -> a;
  int *is_idxs = spec -> is_idxs -> a;
  int *ii_start = spec -> ii_start -> a;

  double **tr;

  /* figure out size of tr */
  for (n_el = 0, j = 0; j < n_sfs; j++) {

    is_pos = ii_start[j];
    is_idx = is_idxs[j];
    /* printf("%le \n", inp -> trs[4][is2fs[j]]); */
    for (n_el++, k = is_pos; ((is_idx = is_idxs[k]) != -1); k++,n_el++) {    /* printf("  %le\n", inp -> trs[4][gs2is[is_idx]]); */}
  }
  /* printf("\n\n" ); */
  if((tr = malloc((n_el+1) * sizeof(double*))) == NULL ) {
    fprintf(stderr, "spectrum.c, function set_trs_red: failed to allocate memory for \"tr\"\n");
    printf("program terminating due to the previous error.\n");
    exit(EXIT_FAILURE);
  }

  for (j = 0; j < n_el+1; j++) {
    if((tr[j] = malloc(4 * sizeof(double))) == NULL ) {
      fprintf(stderr, "spectrum.c, function set_trs_red: failed to allocate memory for \"tr[%d]\"\n"
              , j);
      printf("program terminating due to the previous error.\n");
      exit(EXIT_FAILURE);
    }
  }
  spec -> n_st = n_el;

  for (l = 0, k = 0, j = 0; j < n_sfs; j++) {
    /* printf("%d %d %d %d\n", j, l, n_sfs, spec -> n_st); */
    /*       fflush(stdout); */

    tr[l][0] = inp -> trs[0][is2fs[j]];
    tr[l][1] = inp -> trs[1][is2fs[j]];
    tr[l][2] = inp -> trs[3][is2fs[j]] - inp -> trs[2][is2fs[j]];
    tr[l][3] = inp -> trs[4][is2fs[j]];
    /* printf("%le \n", tr[l][3]); */
    is_pos = ii_start[j];
    is_idx = is_idxs[j];

    for (l++, k = is_pos; ((is_idx = is_idxs[k]) != -1); k++, l++) {
      /* printf("%d %d %d %d %le\n", is_idx, gs2is[is_idx], gs2is[is_idx], inp -> n_trans, inp -> trs[2][gs2is[is_idx]]); */
      /* fflush(stdout); */
      tr[l][0] = get_boltzw((inp -> trs[2][gs2is[is_idx]]
                       - inp -> e0) * (double)AUTOEV);
      tr[l][1] = 0;
      tr[l][2] = inp -> trs[3][gs2is[is_idx]] - inp
        -> trs[2][gs2is[is_idx]];
      tr[l][3] = inp -> trs[4][gs2is[is_idx]];
      /* printf(" %le \n", tr[l][3]); */
    }
    /* for (m = 0; m < l; m++) { */
    /*   printf("tr[%d][0] = %d\n", m, (int)tr[m][1] ); */
    /*   fflush(stdout); */
    /* }  */
  }
  /* printf("\n\n" ); */

  tr[l][1] = -1;

  spec -> trs_red = tr;
  return 0;
}

int
free_spec (struct spectrum *spec)
{
  int j;
  struct spectrum *ls = spec -> last_spec;
  free(spec -> is2fs);
  free(spec -> is_idxs);
  free(spec -> gs2is);
  free(spec -> ii_start);



  if (spec -> idx > 1) {
    for (j = 0; j < spec -> n_elx; j++) {
      free(spec -> omega_x[j]);
      free(spec -> omega_y[j]);
      free(spec -> s_mat[j]);
    }
    free(spec -> s_mat);
    free(spec -> omega_x);
    free(spec -> omega_y);
  }

  /* maintain the structure of the list of spectra */
  if (spec -> last_spec != NULL)  {
    ls -> next_spec = spec -> next_spec;
  }
  if (spec -> next_spec != NULL) {
    spec -> next_spec -> last_spec = ls;
  }

  free(spec);

  return 1;
}

int
free_all_specs (struct inp_node *inp)
{

  struct spectrum *spec;
  struct spectrum *ns = inp -> root_spec;

  while (ns != NULL) {
    spec = ns;
    ns = spec -> next_spec;
    free_spec(spec);
  }

  return 1;
}
