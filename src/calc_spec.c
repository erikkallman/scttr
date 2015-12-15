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
   * @file calc_spec.c
   * @author Erik Källman
   * @date November 2015
   * @brief This file contains the implementation of the calc_spec() function.
   * It calculates a spectrum, stored in the s_mat variable in a
   * spectrum node.
   */
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <signal.h>
#include <time.h>
#include <complex.h>
#include <omp.h>
#include <string.h>
#include <time.h>
#include "calc_spec.h"
#include "dyn_array.h"
#include "scttr_io.h"
#include "sci_const.h"
#include "std_num_ops.h"
#include "std_char_ops.h"
#include "spectrum_s.h"
#include "spectrum.h"
#include "inp_node_s.h"
#include "cache_opt.h"
#include "cpu_opt.h"

int
intf_0( struct inp_node *inp, struct spectrum *spec, struct metadata *md)
{

  int j, k, l, x, y; /* iteration variables */
  int nth = env2int("OMP_NUM_THREADS");

  /* create local copies and pointers to variables whos name would
     otherwise clutter loops */
  int n_st = spec -> n_st;
  int rchunk,rem;
  int * rchunks = malloc(nth * sizeof(int));
  memset(rchunks, 0, nth*sizeof(nth));

  double ** sm_th0; /* thread-shared pointer to the spectrum matrix of thread 1 */

  double ** tr = spec -> trs_red;

  double o_x, o_y; /* energies at each point in the spectrum matrix  */

  double fwhm_x, fwhm_y; /* full-width half-maxima */

  /* excitation energy (x-axis) broadening parameters */
  double grms_x; /* gaussian RMS value */
  double gvar_x; /* square root of the variance of the gaussian */

  /* transfer energy (y-axis) broadening parameters */
  double grms_y; /* gaussian RMS value */
  double gvar_y; /* square root of the variance of the gaussian */

  double emin_x = spec -> emin_x;
  double emin_y = spec -> emin_y;

  /* element to element energy difference in the s_mat matrix */
  double de_x = md -> res[0] / AUTOEV;
  double de_y = md -> res[1] / AUTOEV;

  double wt; /* wall-time counter */

  fwhm_x = md -> fwhm[0] / AUTOEV;
  fwhm_y = md -> fwhm[1] / AUTOEV;
  grms_x = 2.0 * powerl((fwhm_x / (2 * sqrt(2 * log(2)))), 2);
  gvar_x = fwhm_x / (2.0 * sqrt(2.0 * log(2))) * sqrt(2.0 * 3.1415927);

  grms_y = 2.0 * powerl((fwhm_y / (2 * sqrt(2 * log(2)))), 2);
  gvar_y = fwhm_y / (2.0 * sqrt(2.0 * log(2))) * sqrt(2.0 * 3.1415927);

  /* how many rows in tr to assign to each thread due to cache constraints */
  spec -> prsz = cache_cfg -> l_sz / sizeof(double);
  spec -> npr = (int)floorf(spec -> n_ely / spec -> prsz) + (spec -> n_ely % spec -> prsz  > 1 ? 1 : 0);
  spec -> npr_tot = spec -> n_elx * spec -> npr;
  /* printf("x = %d, %d %d\n", x, spec->npr_tot / spec -> npr, spec -> n_elx); */
  /* printf("j = %d, %d %d\n", j, spec->npr * spec->prsz, spec -> n_ely); */
  /* fprintf(stderr, "\n\n=======Valgrind eject point=======\n\n"); */
  /* exit(1); */

  if ( spec-> npr * spec->prsz < spec -> n_ely) {
    fprintf(stderr, "calc_spec.c, function intf_0: matrix incorrectly partitioned (spec-> npr * spec->prsz < spec -> n_ely)\n");
    printf( "program terminating due to the previous error.\n");
    exit(EXIT_FAILURE);
  }

  rchunk = get_row_chunk(4, sizeof(double), cache_cfg);
  if (rchunk > n_st/nth) {
    /* if the amount of data we can assign to a given thread due to cache constraints is larger than the total number of rows divided on threads, reset the rchunk, since otherwise, one thread might carry a significantly larger portion of the work load. also, we know that n_st/nth is smaller than rchunk at this point, meaning that it will fit the cache constraints. */
    rchunk = n_st/nth;
  }

  rem = n_st % rchunk;

  /* if there is a remainder, split it up among the threads as well as possible */
  if (rem > 0) {
    if (rem % nth > 0) {
      for (j = 0; rem > 0; rem--, j++) {
        rchunks[j] += 1;
      }
    }
    else {
      for (j = 0; j < nth; j++) {
        rchunks[j] += rem / nth;
      }
    }
  }

  for (j = 0; j < nth; j++) {
    rchunks[j] += rchunk;
    /* printf("pre chu%d\n", rchunks[j]); */
    /* fflush(stdout); */
  }

  /* divide up the transitions in the tr matrix  */
  for (j = 0, k = 0, l = 0; j < n_st; j++) {

    /* printf("%le\n", tr[j][3]); */
    if (k >= rchunks[l]) {
      l++;
      rchunks[l] = k;
      k = 0;
    }
    while((int)tr[++j][1] == 0){
      /* printf("  %le\n", tr[j][3]); */
    }
    j--;
  }

  /* for (j = 0; j < nth; j++) { */
  /*   printf("post chu%d\n", rchunks[j]); */
  /*   fflush(stdout); */
  /* } */

  /* printf("%d %d %d %d\n", rchunk, n_st/nth, rchunk%n_st/nth, rem); */
  /* fprintf(stderr, "\n\n=======Valgrind eject point=======\n\n"); */
  /* exit(1); */
  for (x = 0; x < spec -> n_elx; x++) {
    o_x = emin_x + (x * de_x);
    for (y = 0; y < spec -> n_ely; y++) {
      o_y = emin_y + (y * de_y);
      spec -> omega_x[x][y] = o_x;
      spec -> omega_y[x][y] = o_y;
    }
  }

  /* printf("%le to %le, %le to %le \n", spec -> omega_x[0][0] * AUTOEV,spec -> omega_x[x-1][y-1] * AUTOEV, spec -> omega_y[0][0]* AUTOEV,  spec -> omega_y[x-1][y-1]* AUTOEV); */
  /* printf("%d %d %d\n", spec -> prsz, spec -> npr, spec -> npr_tot / spec -> npr); */
  /* printf("%d %d %d\n", spec -> n_elx, spec -> n_ely, spec -> npr_tot / spec -> npr); */
  /* fprintf(stderr, "\n\n=======Valgrind eject point=======\n\n"); */
  /* exit(1); */
  wt = omp_get_wtime();
#pragma omp parallel num_threads(env2int("OMP_NUM_THREADS"))
  {
    /* looping variables for the partitioned spectrum matrix  */
    int j; /* row index */
    int k; /* column index */

    int l; /* row index for tr matrix */

    int x, y; /* row and column index for the plotted spectrum */
    int ith = omp_get_thread_num();

    int l_st = 0;

    for (j = 0; j != ith; j++) {
      l_st += rchunks[j];
    }

    /* thread-local start and finish index in the tr matrix */
    int l_fn = l_st + rchunks[ith];

    /* energies and transition moments along the axis .. */
    double ediff_x, tmom_gi; /* .. of excitiation (x-axis)*/
    double ediff_y, tmom_if; /* .. of energy transfer (y-axis)*/
    double de_gi, de_if; /* energy eigenvalue differences */
    double bw; /* boltzmann weight */
    double omega_x;

    double complex tmp;  /* accumulator used in the Kramers-Heisenberg formula */

    double ** sm;

    if((sm = malloc(spec -> npr_tot * sizeof(double *))) == NULL ) {
      fprintf(stderr, "calc_spec.c, function intf_0: failed to allocate memory for \"sm\"\n");
      printf("program terminating due to the previous error.\n");
      exit(EXIT_FAILURE);
    }

    for (j = 0; j < spec -> npr_tot; j++) {
      if((sm[j] = malloc(spec -> prsz * sizeof(double))) == NULL ) {
        fprintf(stderr, "calc_spec.c, function intf_0: failed to allocate memory for \"sm[%d]\"\n",j);
        printf("program terminating due to the previous error.\n");
        exit(EXIT_FAILURE);
      }
      memset(sm[j], 0, spec -> prsz * sizeof(double));
    }

    if (ith == 0) {
      sm_th0 = sm;
    }

    for (j = 0, x = 0, y = 0; x < spec -> n_elx/* j < spec -> npr_tot */; j++) {
      for (k = 0; (k < spec -> prsz) && (x < spec -> n_elx); k++, y++) {
        omega_x = emin_x + (x * de_x);
        for (l = l_st; l < l_fn; l++) {
          de_if = tr[l][2];
          tmom_if = tr[l][3];
          tmp = 0 + 0*I;
          while((int)tr[++l][1] == 0) { /* loop over ground to intermediate
                                          transitions */
            bw = tr[l][0];
            de_gi = tr[l][2];
            tmom_gi = tr[l][3];
            /* printf("%d %d %d %d\n", x, spec -> n_elx, y, spec -> n_ely); */
            /* fflush(stdout); */
            ediff_x = omega_x - de_gi;
            ediff_y = (emin_y + (y * de_y)) - de_gi - de_if;

            tmp += tmom_gi * tmom_if * bw / (-de_gi + omega_x
                                             - (fwhm_x / 2)*I);

            tmp *= (exp(-(powerl(ediff_x, 2)) / grms_x) / gvar_x * de_x)
              * (exp(-(powerl(ediff_y, 2)) / grms_y) / gvar_y * de_y);
          }
          tmp = fabsc(tmp);
          tmp *= tmp;
          sm[j][k] += creal(tmp);
          l--;
        }

        if (y == spec -> n_ely-1) {
          /* we have traversed one row in the spectrum */
          /* printf("%d %d %d\n", j*spec->prsz+k, y, spec -> n_ely); */
          /* fprintf(stderr, "\n\n=======Valgrind eject point=======\n\n"); */
          /* exit(1); */
          x++;
          y = 0;
        }
      }
    }

    /* printf("%d %d\n",x,y); */
    /* printf("x = %d, %d %d\n", x, spec->npr_tot / spec -> prsz, spec -> n_elx); */
    /* printf("j = %d, %d %d\n", j, spec->npr * spec->prsz, spec -> n_ely); */
    /* printf("tot = %d, %d %d\n", spec->n_elx * spec->n_ely, j * spec -> prsz, spec->npr_tot * spec -> prsz); */
    /* printf("%le to %le, %le to %le \n", spec -> omega_x[0][0] * AUTOEV,spec -> omega_x[x-1][y-1] * AUTOEV, spec -> omega_y[0][0]* AUTOEV,  spec -> omega_y[x-1][y-1]* AUTOEV); */
    /* fflush(stdout); */
    /* fprintf(stderr, "\n\n=======Valgrind eject point=======\n\n"); */
    /* exit(1); */
    /* printf("\n%d \n%d \n%d \n%d \n%d \n%d \n%d \n%d\n", j, spec -> npr_tot, k, spec -> prsz, y, spec -> n_ely, x, spec -> n_elx ); */
    /* fflush(stdout); */
    /* fprintf(stderr, "\n\n=======Valgrind1 eject point=======\n\n"); */
    /* exit(1); */

#pragma omp barrier

#pragma omp critical
    {
      if (ith != 0) {
        /* printf("\n\nthread %d summing\n\n", ith); */
        /* fflush(stdout); */
        for (j = 0, l = 0; j < spec -> npr_tot; j++) {
          for (k = 0; k < spec -> prsz; k++, l++) {

            sm_th0[j][k] += sm[j][k];

            /* if (l == n_elx-1) { */
            /*   /\* we have traversed one row in the spectrum *\/ */
            /*   fprintf(stderr, "\n\n=======Valgrind eject point=======\n\n"); */
            /*   exit(1); */
            /*   l = 0; */
            /*   break; */
            /* } */
          }
        }
        for (j = 0; j < spec -> npr_tot; j++) {
          free(sm[j]);
        }
        free(sm);
      }
    }
  }
  wt = omp_get_wtime() - wt;

  FILE *cpu_out;
  cpu_out = fopen("./cpu.dat", "a");
  fprintf(cpu_out,"%le\n",wt);
  fclose(cpu_out);

  spec -> s_mat = sm_th0;
  free(rchunks);
  return 0;
}


int
intf_0_old( struct inp_node *inp, struct spectrum *spec, struct metadata *md)
{

  int j, k, x, y;

  /* create local copies and pointers to variables whos name would
     otherwise clutter the code loop */
  int n_sfs = spec -> is2fs -> n_el;
  int *is2fs = spec -> is2fs -> a;
  int *gs2is = spec -> gs2is -> a;
  int *is_idxs = spec -> is_idxs -> a;
  int *ii_start = spec -> ii_start -> a;
  int is_idx, is_pos;

  /* energies and transition moments along the axis .. */
  double ediff_x, tmom_gi; /* .. of excitiation (x-axis)*/
  double ediff_y, tmom_if; /* .. of energy transfer (y-axis)*/
  double de_gi, de_if; /* energy eigenvalue differences */
  double o_x, o_y; /* energies at each point in the spectrum matrix  */
  double bw; /* boltzmann weight */
  double fwhm_x, fwhm_y; /* full-width half-maxima */
  double complex tmp;  /* accumulator used in the Kramers-Heisenberg formula */

  /* excitation energy (x-axis) broadening parameters */
  double grms_x; /* gaussian RMS value */
  double gvar_x; /* square root of the variance of the gaussian */

  /* transfer energy (y-axis) broadening parameters */
  double grms_y; /* gaussian RMS value */
  double gvar_y; /* square root of the variance of the gaussian */

  double emin_x = spec -> emin_x;
  double emin_y = spec -> emin_y;

  /* element to element energy difference in the s_mat matrix */
  double de_x = md -> res[0] / AUTOEV;
  double de_y = md -> res[1] / AUTOEV;

  fwhm_x = md -> fwhm[0] / AUTOEV;
  fwhm_y = md -> fwhm[1] / AUTOEV;
  grms_x = 2.0 * powerl((fwhm_x / (2 * sqrt(2 * log(2)))), 2);
  gvar_x = fwhm_x / (2.0 * sqrt(2.0 * log(2))) * sqrt(2.0 * 3.1415927);

  grms_y = 2.0 * powerl((fwhm_y / (2 * sqrt(2 * log(2)))), 2);
  gvar_y = fwhm_y / (2.0 * sqrt(2.0 * log(2))) * sqrt(2.0 * 3.1415927);
  for (x = 0; x < spec -> n_elx; x++) {

    o_x = emin_x + (x * de_x);
    for (y = 0; y < spec -> n_ely; y++) {
      spec -> s_mat[x][y] = 0;
      o_y = emin_y + (y * de_y);
      spec -> omega_x[x][y] = o_x;
      spec -> omega_y[x][y] = o_y;

      for (k = 0, j = 0; j < n_sfs; j++) {
        is_pos = ii_start[j];
        tmom_if = inp -> trs[4][is2fs[j]];
        de_if = inp -> trs[3][is2fs[j]] - inp -> trs[2][is2fs[j]];

        is_idx = is_idxs[j];
        tmp = 0 + 0*I;

        for (k = is_pos; ((is_idx = is_idxs[k]) != -1); k++) {

          tmom_gi = inp -> trs[4][gs2is[is_idx]];
          de_gi = inp -> trs[3][gs2is[is_idx]] - inp
            -> trs[2][gs2is[is_idx]];

          bw = get_boltzw((inp -> trs[2][gs2is[is_idx]]
                           - inp -> e0) * (double)AUTOEV);

          ediff_x = spec -> omega_x[x][y] - de_gi;
          ediff_y = spec -> omega_y[x][y] - de_gi - de_if;

          tmp += tmom_gi * tmom_if * bw / (-de_gi + spec -> omega_x[x][y]
                                           - (fwhm_x / 2)*I);

          tmp *= (exp(-(powerl(ediff_x, 2)) / grms_x) / gvar_x * de_x)
            * (exp(-(powerl(ediff_y, 2)) / grms_y) / gvar_y * de_y);
        }
        tmp = fabsc(tmp);
        tmp *= tmp;
        spec -> s_mat[x][y] += creal(tmp);
      }
    }
  }

  return 0;
}


int
calc_spec (struct inp_node *inp, int spec_idx)
{
  /* int x, y; */

  struct metadata *md = inp -> md;
  struct spectrum *spec = get_spec(inp, spec_idx);

  int j, k;

  double emin_x = spec -> emin_x;
  double emax_x = spec -> emax_x;
  double emin_y = spec -> emin_y;
  double emax_y = spec -> emax_y;

  spec -> n_elx = (int)((emax_x - emin_x) / (md -> res[0] / AUTOEV));
  spec -> n_ely = (int)((emax_y - emin_y) / (md -> res[1] / AUTOEV));

  if((spec -> omega_x = malloc(spec -> n_elx * sizeof(double *))) == NULL ) {
    fprintf(stderr, "calc_spec.c, function calc_spec: failed to allocate memory for \"spec -> omega_x\"\n");
    printf("program terminating due to the previous error.\n");
    exit(1);
  }

  if((spec -> omega_y = malloc(spec -> n_elx * sizeof(double *))) == NULL ) {
    fprintf(stderr, "calc_spec.c, function calc_spec: failed to allocate memory for \"spec -> omega_y\"\n");
    printf("program terminating due to the previous error.\n");
    exit(1);
  }

  if((spec -> s_mat = malloc(spec -> n_elx * sizeof(double *))) == NULL ) {
    fprintf(stderr, "calc_spec.c, function calc_spec: failed to allocate memory for \"s_mat\"\n");
    printf("program terminating due to the previous error.\n");
    exit(1);
  }

  for (j=0; j<spec -> n_elx; j++) {
    if((spec -> s_mat[j] = malloc(spec -> n_ely * sizeof(double))) == NULL ) {
      fprintf(stderr, "calc_spec.c, function calc_spec: failed to allocate memory for \"spec -> s_mat[%d]\"\n"
              , j);
      printf("program terminating due to the previous error.\n");
      exit(1);
    }

    if((spec -> omega_x[j] = malloc(spec -> n_ely * sizeof(double))) == NULL ) {
      fprintf(stderr, "calc_spec.c, function calc_spec: failed to allocate memory for \"spec -> omega_x[%d]\"\n"
              , j);
      printf("program terminating due to the previous error.\n");
      exit(1);
    }

    if((spec -> omega_y[j] = malloc(spec -> n_ely * sizeof(double))) == NULL ) {
      fprintf(stderr, "calc_spec.c, function calc_spec: failed to allocate memory for \"spec -> omega_y[%d]\"\n"
              , j);
      printf("program terminating due to the previous error.\n");
      exit(1);
    }
  }

  printf("  - calculating the RIXS map, using");
  switch(inp -> md -> intf_mode)
    {

    case 1:
      printf(" a constructive interference model .. \n");
      intf_0(inp, spec, md);

      for (j = 0; j < spec -> npr_tot; j++) {
        for (k = 0; k < spec -> prsz; k++) {
          if (spec -> s_mat[j][k] > spec -> sfac) {
            spec -> sfac = spec -> s_mat[j][k];
          }
        }
      }

      write_spec(inp, get_spec(inp,2));
      break;

    default:
      printf(" a constructive interference model, since you didnt specify this parameter .. \n");
      intf_0_old(inp, spec, md);

      for (j = 0; j < spec -> npr_tot; j++) {
        for (k = 0; k < spec -> prsz; k++) {
          if (spec -> s_mat[j][k] > spec -> sfac) {
            spec -> sfac = spec -> s_mat[j][k];
          }
        }
      }

      write_spec_old(inp, get_spec(inp,2));
    }

  printf("    .. done.\n");
  printf("  - finding the normalization constant of the  calculated RIXS map ..");

  /* for (j = 0; j < spec -> npr_tot; j++) { */
  /*   for (k = 0; k < spec -> prsz; k++) { */
  /*     if (spec -> s_mat[j][k] > spec -> sfac) { */
  /*       spec -> sfac = spec -> s_mat[j][k]; */
  /*     } */
  /*   } */
  /* } */

  /* for (x = 0; x < spec -> n_elx; x++) { */
  /*   for (y = 0; y < spec -> n_ely; y++) { */
  /*     if (spec -> s_mat[x][y] > spec -> sfac) { */
  /*       spec -> sfac = spec -> s_mat[x][y]; */
  /*     } */
  /*   } */
  /* } */

  printf(" done.\n");

  return 0;
}
