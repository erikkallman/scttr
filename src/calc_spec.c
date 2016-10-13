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
#include "timing.h"
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
#include "glob_time.h"

double para_t;

int
intf_0( struct inp_node *inp, struct spectrum *spec, struct metadata *md)
{
  int x, y;
  int j, k, l; /* iteration variables */
  int nth = 0;
  env2int("OMP_NUM_THREADS", &nth);

  /* variables used in the inner loop of thelorentzian convolution*/
  int x_in, y_in; /* row and column index */
  int j_in, k_in;

  int lx_i, ly_i; /* indices of the currently used broadening
                                   values */
  lx_i = ly_i = - 1;

  /* create local copies and pointers to variables whos name would
     otherwise clutter loops */
  int n_st = spec -> n_st;
  int rchunk,rem;
  int * rchunks = malloc(nth * sizeof(int));
  memset(rchunks, 0, nth*sizeof(nth));

  char *ltime = malloc(20);

  double tmp_int; /* variable used to accumulate intensites */
  double omega_x, omega_y;
  double omega_x_in, omega_y_in;
  double lx_t, ly_t;
  double eu_lx, eu_ly;
  double emin_x = spec -> emin_x;
  double emin_y = spec -> emin_y;
  eu_lx = -fabs(emin_x * 2);
  eu_ly = -fabs(emin_y * 2);
  double ediff_x, ediff_y;

  double lx_mfac, lx_afac;
  double ly_mfac, ly_afac;

  double ** sm_th0; /* thread-shared pointer to the spectrum matrix of thread 1 */

  double ** tr = spec -> trs_red;

  /* element to element energy difference in the s_mat matrix */
  double de_x = md -> res[0] / AUTOEV;
  double de_y = md -> res[1] / AUTOEV;

  double wt; /* wall-time counter */

  /* how many rows in tr to assign to each thread due to cache constraints? */
  spec -> prsz = cache_cfg -> l_sz / sizeof(double);
  spec -> npr = (int)floorf(spec -> n_ely / spec -> prsz) + (spec -> n_ely % spec -> prsz  > 1 ? 1 : 0);
  spec -> npr_tot = spec -> n_elx * spec -> npr;

  if ( spec-> npr * spec->prsz < spec -> n_ely) {
    fprintf(stderr, "calc_spec.c, function intf_0: matrix incorrectly partitioned (spec-> npr * spec->prsz < spec -> n_ely)\n");
    printf( "program terminating due to the previous error.\n");
    exit(EXIT_FAILURE);
  }

  rchunk = get_row_chunk(4, sizeof(double), cache_cfg);
  goto threading;
 threading:
  memset(rchunks, 0, nth*sizeof(nth));
  if (rchunk > n_st/nth) {
    /* if the amount of data we can assign to a given thread due to cache constraints is larger than the total number of rows divided on threads, reset the rchunk, since otherwise, one thread might carry a significantly larger portion of the work load. also, we know that n_st/nth is smaller than rchunk at this point, meaning that it will fit the cache constraints. */
    rchunk = n_st/nth;
  }
  else {
    /* since the memory chunk size is smaller, is it even possible to use it? */
    if (n_st/rchunk > nth) {
      /* no, increase the chunk size so that the total number of transitions */
      rchunk = n_st / nth;
    }
  }

  rem = n_st % rchunk;

  /* if there is a remainder, split it up among the threads as well as possible */
  /* printf("pre rem = %d\n", rem); */
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
  /* printf("post rem = %d\n", rem); */

  /*   printf("\nCHUNKS PRE0 %d %d\n\n", n_st, k); */
  /* for (j = 0; j < nth; j++) { */
  /*   printf("rchunks[%d] = %d\n",j, rchunks[j]); */
  /*   fflush(stdout); */
  /* } */

  k = 0;
  for (j = 0; j < nth; j++) {
    rchunks[j] += rchunk;
    k += rchunks[j];
  }

  /* printf("\nCHUNKS PRE1 %d %d\n\n", n_st, k); */
  /* for (j = 0; j < nth; j++) { */
  /*   printf("rchunks[%d] = %d\n",j, rchunks[j]); */
  /*   fflush(stdout); */
  /* } */

  /* divide up the transitions in the tr matrix  */
  j = k = l = 0;
  while (l != nth) {
    if (j >= n_st){
      rchunks[l] = k;
      l++;
      k = 0;
      break;
    }
    else if (k >= rchunks[l]) {
      rchunks[l] = k;
      l++;
      k = 0;
    }
    for (k++; (++j < n_st) && ((int)tr[j][1] == 0); k++){};
  }

  /* distribute the remainder */
  rem = n_st - j - 1;
  if (rem > 0) {
    /* if there is a remainder, the last control loop exited since l == nth */
    rchunks[l-1] += rem;
  }

  if (l < nth) {
    printf("========Does not scale beyond %d (nth = %d)\n", l ,nth );
    nth = l;
    goto threading;
  }

  /* printf("\nCHUNKS POST\n\n"); */
  /* for (j = 0; j < nth; j++) { */
  /*   printf("rchunks[%d] = %d\n",j, rchunks[j]); */
  /*   fflush(stdout); */
  /* } */

  /* j = k = l = 0; */
  /* while (j<n_st) { */
  /*   printf("tr[%d] = %le\n",j , tr[j][3]); */
  /*   while((int)tr[++j][1] == 0){ */
  /*     printf("  tr[%d] = %le\n", j, tr[j][3]); */
  /*   } */
  /* } */

  /* printf("%d %d %d %d\n", nth, n_st, rchunk%n_st/nth, rem); */
  /* fflush(stdout); */

  /* printf("%le to %le, %le to %le \n", spec -> omega_x[0][0] * AUTOEV,spec -> omega_x[x-1][y-1] * AUTOEV, spec -> omega_y[0][0]* AUTOEV,  spec -> omega_y[x-1][y-1]* AUTOEV); */
  /* printf("%d %d %d\n", spec -> prsz, spec -> npr, spec -> npr_tot / spec -> npr); */
  /* printf("%d %d %d\n", spec -> n_elx, spec -> n_ely, spec -> npr_tot / spec -> npr); */
  /* printf("==== THREAD INFO START ==== \n\n"); */
  /* fflush(stdout); */
  if((spec -> s_mat = malloc(spec -> npr_tot * sizeof(double *))) == NULL ) {
    fprintf(stderr, "calc_spec.c, function intf_0: failed to allocate memory for \"spec -> s_mat\"\n");
    printf("program terminating due to the previous error.\n");
    exit(EXIT_FAILURE);
  }

  for (j = 0; j < spec -> npr_tot; j++) {
    if((spec -> s_mat[j] = malloc(spec -> prsz * sizeof(double))) == NULL ) {
      fprintf(stderr, "calc_spec.c, function intf_0: failed to allocate memory for \"spec -> s_mat[%d]\"\n",j);
      printf("program terminating due to the previous error.\n");
      exit(EXIT_FAILURE);
    }
    memset(spec -> s_mat[j], 0, spec -> prsz * sizeof(double));
  }

  wt = omp_get_wtime();

#pragma omp parallel num_threads(nth)
  {
    /* looping variables for the partitioned spectrum matrix  */
    int j; /* row index */
    int k; /* column index */
    int l; /* row index for tr matrix */

    int x, y; /* row and column index for the plotted spectrum */
    int gx_i, gy_i; /* indices of the currently used broadening
                       values */
    gx_i = gy_i = -1;
    int ith = omp_get_thread_num();

    /* thread 0 carries an array of pointers to the spectru matrices
     of all other threads, to later perform summation of the individual
    matrices. */
    int l_st = 0;

    for (j = 0; j != ith; j++) {
      l_st += rchunks[j];
    }

    /* thread-local start and finish index in the tr matrix */
    int l_fn = l_st + rchunks[ith];
    /* printf("started calculating: thread %d, start %d, end %d\n", ith, l_st, l_fn); */
    /* fflush(stdout); */

    /* energies and transition moments along the axis .. */
    double ediff_x, tmom_gi; /* .. of excitiation (x-axis)*/
    double ediff_y, tmom_if; /* .. of energy transfer (y-axis)*/
    double de_gi, de_if; /* energy eigenvalue differences */
    double bw; /* boltzmann weight */
    double omega_x, omega_y;

  /* excitation energy (x-axis) and transfer energy (y-axis) broadening
     parameters */
    double gx_t, gy_t; /* temporary full-width half-maximum
                                      broadening values */
    double eu_gx, eu_gy; /* upper energy bound for a given broadening value */

    /* initialize them so that the broadening value update will trigger
    at the first iteration */
    eu_gx = -fabs(emin_x * 2);
    eu_gy = -fabs(emin_y * 2);

    /* factors to reduce the calculation of the gaussian and lorentzian
    distributions */
    double gy_stddev, gx_stddev;
    double gy_var, gx_var;
    double gy_mfac, gx_mfac;

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
    omega_x = emin_x;
    for (j = 0, x = 0, y = 0; x < spec -> n_elx/* j < spec -> npr_tot */; j++) {
      for (k = 0; (k < spec -> prsz) && (x < spec -> n_elx); k++, y++) {
        omega_y = emin_y + (y * de_y);

        /* update the broadenings, if needed */
        if (omega_x > eu_gx) {
          /* printf("thread %d set gx since %le > %le\n", ith, omega_x, eu_gx); */
          /* fflush(stdout); */
          gx_i += 2;
          gx_t = md -> gx[gx_i] / AUTOEV;

          gx_stddev = gx_t / (2 * sqrt(2 * log(2)));
          gx_var = 2.0 * powerl(gx_stddev, 2);
          gx_mfac = de_x / gx_stddev * sqrt(2.0 * PI);

          eu_gx = (md -> gx)[gx_i+1] / AUTOEV;
          /* printf("thread %d set eu_gx to %le\n", ith, eu_gx); */
          /* fflush(stdout); */
        }
        if (omega_y > eu_gy) {
          /* printf("thread %d set gy since %le > %le\n", ith, omega_y, eu_gy); */
          /* fflush(stdout); */
          gy_i += 2;
          gy_t = md -> gy[gy_i] / AUTOEV;

          gy_stddev = gy_t / (2 * sqrt(2 * log(2)));
          gy_var = 2.0 * powerl(gy_stddev, 2);
          gy_mfac = de_y / gy_stddev * sqrt(2.0 * PI);

          eu_gy = (md -> gy)[gy_i+1] / AUTOEV;
          /* printf("thread %d set eu_gy to %le\n\n", ith, eu_gy); */
          /* fflush(stdout); */
        }

        for (l = l_st; l < l_fn;) {
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
            ediff_y = omega_y - de_gi - de_if;

            tmp += tmom_gi * tmom_if * bw / (-de_gi + omega_x
                                             - (gx_t / 2)*I);

            tmp *= (exp(-(powerl(ediff_x, 2)) / gx_var) * gx_mfac)
              * (exp(-(powerl(ediff_y, 2)) / gy_var) * gy_mfac);
            /* tmp = fabsc(tmp); */
            /* sm[j][k] += creal(tmp); */
          }
          tmp = fabsc(tmp);
          tmp *= tmp;
          sm[j][k] += creal(tmp);
        }

        if (y == spec -> n_ely-1) {
          /* we have traversed one row in the spectrum */
          /* printf("%d %d %d\n", j*spec->prsz+k, y, spec -> n_ely); */

          x++;
          omega_x = emin_x + (x * de_x);
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
    /* printf("finished calculating: thread %d, start %d, end %d\n", ith, l_st, l_fn); */

#pragma omp barrier
    if (ith == 0) {
      printf("\n      summing up the thread-local spectrum layers.. (%s)",get_loctime(ltime));
      fflush(stdout);
    }
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
      }
    } /* implicit barrier */
  }

  printf(" done (%s).", get_loctime(ltime));
  fflush(stdout);
  para_t = omp_get_wtime() - wt;
  printf("\n      applying lorentzian boadening.. (%s)",get_loctime(ltime));
  fflush(stdout);

  for (j = 0, x = 0, y = 0; x < spec -> n_elx; j++) {
    for (k = 0; (k < spec -> prsz) && (x < spec -> n_elx); k++, y++) {
      omega_x = emin_x + (x * de_x);
      omega_y = emin_y + (y * de_y);

      /* update the broadenings, if needed */
      if (omega_x > eu_lx) {
        lx_i += 2;
        lx_t = md -> lx[lx_i] / AUTOEV;

        lx_mfac = 0.5 * lx_t / PI;
        lx_afac = (0.25 * lx_t * lx_t);

        eu_lx = (md -> lx)[lx_i+1] / AUTOEV;
      }
      if (omega_y > eu_ly) {
        ly_i += 2;
        ly_t = md -> ly[ly_i] / AUTOEV;

        ly_mfac = 0.5 * ly_t / PI;
        ly_afac = (0.25 * ly_t * ly_t);

        eu_ly = (md -> ly)[ly_i+1] / AUTOEV;
      }

      tmp_int = 0;
      for (j_in = 0, x_in = 0, y_in = 0; x_in < spec -> n_elx; j_in++) {
        for (k_in = 0; (k_in < spec -> prsz) && (x_in < spec -> n_elx)
               ;k_in++, y_in++) {
          omega_x_in = emin_x + (x_in * de_x);
          omega_y_in = emin_y + (y_in * de_y);

          /* calculate the lorentzian contribution to this intensity_in */
          ediff_x = omega_x - omega_x_in;
          ediff_y = omega_y - omega_y_in;

          tmp_int += sm_th0[j_in][k_in]
            * (lx_mfac / ((ediff_x * ediff_x) + lx_afac))
            * (ly_mfac / ((ediff_y * ediff_y) + ly_afac));

          if (y_in == spec -> n_ely-1) {
            x_in++;
            y_in = 0;
          }
        }
      }
      spec -> s_mat[j][k] += tmp_int;

      if (y == spec -> n_ely-1) {
        /* we have traversed one row in the spectrum */
        x++;
        y = 0;
      }
    }
  }
  printf(" done (%s).", get_loctime(ltime));

  /* spec -> s_mat = sm_th0; */
  /* printf("\n\n==== THREAD INFO END ==== \n\n"); */
  /* fflush(stdout); */
  /* for (j = 0; j < spec -> npr_tot; j++) { */
  /*   free(sm_th0[j]); */
  /* } */
  /* free(sm_th0); */
  free(ltime);
  free(rchunks);
  return 0;
}

int
intf_0_old( struct inp_node *inp, struct spectrum *spec, struct metadata *md)
{
  int x, y;
  int j, k, l; /* iteration variables */
  int nth = 0;
  env2int("OMP_NUM_THREADS", &nth);

  /* variables used in the inner loop of thelorentzian convolution*/
  int x_in, y_in; /* row and column index */
  int j_in, k_in;

  int lx_i, ly_i; /* indices of the currently used broadening
                                   values */
  lx_i = ly_i = - 1;

  /* create local copies and pointers to variables whos name would
     otherwise clutter loops */
  int n_st = spec -> n_st;
  int rchunk,rem;
  int * rchunks = malloc(nth * sizeof(int));
  memset(rchunks, 0, nth*sizeof(nth));

  char *ltime = malloc(20);

  double tmp_int; /* variable used to accumulate intensites */
  double omega_x, omega_y;
  double omega_x_in, omega_y_in;
  double lx_t, ly_t;
  double eu_lx, eu_ly;
  double emin_x = spec -> emin_x;
  double emin_y = spec -> emin_y;
  eu_lx = -fabs(emin_x * 2);
  eu_ly = -fabs(emin_y * 2);
  double ediff_x, ediff_y;

  double lx_mfac, lx_afac;
  double ly_mfac, ly_afac;

  double ** sm_th0; /* thread-shared pointer to the spectrum matrix of thread 1 */

  double ** tr = spec -> trs_red;

  /* element to element energy difference in the s_mat matrix */
  double de_x = md -> res[0] / AUTOEV;
  double de_y = md -> res[1] / AUTOEV;

  double wt; /* wall-time counter */

  /* how many rows in tr to assign to each thread due to cache constraints? */
  spec -> prsz = cache_cfg -> l_sz / sizeof(double);
  spec -> npr = (int)floorf(spec -> n_ely / spec -> prsz) + (spec -> n_ely % spec -> prsz  > 1 ? 1 : 0);
  spec -> npr_tot = spec -> n_elx * spec -> npr;

  if (spec-> npr * spec->prsz < spec -> n_ely) {
    fprintf(stderr, "calc_spec.c, function intf_0_old: matrix incorrectly partitioned (spec-> npr * spec->prsz < spec -> n_ely)\n");
    printf( "program terminating due to the previous error.\n");
    exit(EXIT_FAILURE);
  }

  rchunk = get_row_chunk(4, sizeof(double), cache_cfg);
  goto threading;
 threading:
  memset(rchunks, 0, nth*sizeof(nth));
  if (rchunk > n_st/nth) {
    /* if the amount of data we can assign to a given thread due to cache constraints is larger than the total number of rows divided on threads, reset the rchunk, since otherwise, one thread might carry a significantly larger portion of the work load. also, we know that n_st/nth is smaller than rchunk at this point, meaning that it will fit the cache constraints. */
    rchunk = n_st/nth;
  }
  else {
    /* since the memory chunk size is smaller, is it even possible to use it? */
    if (n_st/rchunk > nth) {
      /* no, increase the chunk size so that the total number of transitions */
      rchunk = n_st / nth;
    }
  }

  rem = n_st % rchunk;

  /* if there is a remainder, split it up among the threads as well as possible */
  /* printf("pre rem = %d\n", rem); */
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
  /* printf("post rem = %d\n", rem); */

  /*   printf("\nCHUNKS PRE0 %d %d\n\n", n_st, k); */
  /* for (j = 0; j < nth; j++) { */
  /*   printf("rchunks[%d] = %d\n",j, rchunks[j]); */
  /*   fflush(stdout); */
  /* } */

  k = 0;
  for (j = 0; j < nth; j++) {
    rchunks[j] += rchunk;
    k += rchunks[j];
  }

  /* printf("\nCHUNKS PRE1 %d %d\n\n", n_st, k); */
  /* for (j = 0; j < nth; j++) { */
  /*   printf("rchunks[%d] = %d\n",j, rchunks[j]); */
  /*   fflush(stdout); */
  /* } */

  /* divide up the transitions in the tr matrix  */
  j = k = l = 0;
  while (l != nth) {
    if (j >= n_st){
      rchunks[l] = k;
      l++;
      k = 0;
      break;
    }
    else if (k >= rchunks[l]) {
      rchunks[l] = k;
      l++;
      k = 0;
    }
    for (k++; (++j < n_st) && ((int)tr[j][1] == 0); k++){};
  }

  /* distribute the remainder */
  rem = n_st - j - 1;
  if (rem > 0) {
    /* if there is a remainder, the last control loop exited since l == nth */
    rchunks[l-1] += rem;
  }

  if (l < nth) {
    printf("\n    - the calculation does not scale beyond %d threads (current nth = %d). Reducing the number of threads used to %d.\n", l , nth, l );
    fflush(stdout);
    nth = l;
    goto threading;
  }

  /* printf("\nCHUNKS POST\n\n"); */
  /* for (j = 0; j < nth; j++) { */
  /*   printf("rchunks[%d] = %d\n",j, rchunks[j]); */
  /*   fflush(stdout); */
  /* } */

  /* j = k = l = 0; */
  /* while (j<n_st) { */
  /*   printf("tr[%d] = %le\n",j , tr[j][3]); */
  /*   while((int)tr[++j][1] == 0){ */
  /*     printf("  tr[%d] = %le\n", j, tr[j][3]); */
  /*   } */
  /* } */

  /* printf("%d %d %d %d\n", nth, n_st, rchunk%n_st/nth, rem); */
  /* fflush(stdout); */

  /* printf("%le to %le, %le to %le \n", spec -> omega_x[0][0] * AUTOEV,spec -> omega_x[x-1][y-1] * AUTOEV, spec -> omega_y[0][0]* AUTOEV,  spec -> omega_y[x-1][y-1]* AUTOEV); */
  /* printf("%d %d %d\n", spec -> prsz, spec -> npr, spec -> npr_tot / spec -> npr); */
  /* printf("%d %d %d\n", spec -> n_elx, spec -> n_ely, spec -> npr_tot / spec -> npr); */
  /* printf("==== THREAD INFO START ==== \n\n"); */
  /* fflush(stdout); */
  if((spec -> s_mat = malloc(spec -> npr_tot * sizeof(double *))) == NULL ) {
    fprintf(stderr, "calc_spec.c, function intf_0: failed to allocate memory for \"spec -> s_mat\"\n");
    printf("program terminating due to the previous error.\n");
    exit(EXIT_FAILURE);
  }

  for (j = 0; j < spec -> npr_tot; j++) {
    if((spec -> s_mat[j] = malloc(spec -> prsz * sizeof(double))) == NULL ) {
      fprintf(stderr, "calc_spec.c, function intf_0: failed to allocate memory for \"spec -> s_mat[%d]\"\n",j);
      printf("program terminating due to the previous error.\n");
      exit(EXIT_FAILURE);
    }
    memset(spec -> s_mat[j], 0, spec -> prsz * sizeof(double));
  }

  wt = omp_get_wtime();

#pragma omp parallel num_threads(nth)
  {
    /* looping variables for the partitioned spectrum matrix  */
    int j; /* row index */
    int k; /* column index */
    int l; /* row index for tr matrix */

    int x, y; /* row and column index for the plotted spectrum */
    int gx_i, gy_i; /* indices of the currently used broadening
                       values */
    gx_i = gy_i = -1;
    int ith = omp_get_thread_num();

    /* thread 0 carries an array of pointers to the spectru matrices
     of all other threads, to later perform summation of the individual
    matrices. */
    int l_st = 0;

    for (j = 0; j != ith; j++) {
      l_st += rchunks[j];
    }

    /* thread-local start and finish index in the tr matrix */
    int l_fn = l_st + rchunks[ith];
    /* printf("started calculating: thread %d, start %d, end %d\n", ith, l_st, l_fn); */
    /* fflush(stdout); */

    /* energies and transition moments along the axis .. */
    double ediff_x, tmom_gi; /* .. of excitiation (x-axis)*/
    double ediff_y, tmom_if; /* .. of energy transfer (y-axis)*/
    double de_gi, de_if; /* energy eigenvalue differences */
    double bw; /* boltzmann weight */
    double omega_x, omega_y;

  /* excitation energy (x-axis) and transfer energy (y-axis) broadening
     parameters */
    double gx_t, gy_t; /* temporary full-width half-maximum
                                      broadening values */
    double eu_gx, eu_gy; /* upper energy bound for a given broadening value */

    /* initialize them so that the broadening value update will trigger
    at the first iteration */
    eu_gx = -fabs(emin_x * 2);
    eu_gy = -fabs(emin_y * 2);

    /* factors to reduce the calculation of the gaussian and lorentzian
    distributions */
    double gy_stddev, gx_stddev;
    double gy_var, gx_var;
    double gy_mfac, gx_mfac;

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

    omega_x = emin_x;
    for (j = 0, x = 0, y = 0; x < spec -> n_elx/* j < spec -> npr_tot */; j++) {
      for (k = 0; (k < spec -> prsz) && (x < spec -> n_elx); k++, y++) {
        omega_y = emin_y + (y * de_y);

        /* update the broadenings, if needed */
        if (omega_x > eu_gx) {
          /* printf("thread %d set gx since %le > %le\n", ith, omega_x, eu_gx); */
          /* fflush(stdout); */
          gx_i += 2;
          gx_t = md -> gx[gx_i] / AUTOEV;

          gx_stddev = gx_t / (2 * sqrt(2 * log(2)));
          gx_var = 2.0 * powerl(gx_stddev, 2);
          gx_mfac = de_x / gx_stddev * sqrt(2.0 * PI);

          eu_gx = (md -> gx)[gx_i+1] / AUTOEV;
          /* printf("thread %d set eu_gx to %le\n", ith, eu_gx); */
          /* fflush(stdout); */
        }
        if (omega_y > eu_gy) {
          /* printf("thread %d set gy since %le > %le\n", ith, omega_y, eu_gy); */
          /* fflush(stdout); */
          gy_i += 2;
          gy_t = md -> gy[gy_i] / AUTOEV;

          gy_stddev = gy_t / (2 * sqrt(2 * log(2)));
          gy_var = 2.0 * powerl(gy_stddev, 2);
          gy_mfac = de_y / gy_stddev * sqrt(2.0 * PI);

          eu_gy = (md -> gy)[gy_i+1] / AUTOEV;
          /* printf("thread %d set eu_gy to %le\n\n", ith, eu_gy); */
          /* fflush(stdout); */
        }

        for (l = l_st; l < l_fn;) {
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
            ediff_y = omega_y - de_gi - de_if;

            tmp += tmom_gi * tmom_if * bw / (ediff_x
                                             - (gx_t / 2)*I);

            tmp *= (exp(-(powerl(ediff_x, 2)) / gx_var) * gx_mfac)
              * (exp(-(powerl(ediff_y, 2)) / gy_var) * gy_mfac);
            tmp = fabsc(tmp);
            sm[j][k] += creal(tmp);
            /* printf("%le %le %le %le %le\n ", sm[j][k], tmom_gi, tmom_if, ediff_x, ediff_y); */
            /* printf("%le %le %le %le %le\n ", sm[j][k], tmom_gi, tmom_if, exp(-(powerl(ediff_x, 2)) / gx_var) * gx_mfac, exp(-(powerl(ediff_y, 2)) / gy_var) * gy_mfac); */

            /* printf("x %d %d:%le %le %le %le %le\n ",tr[0][1]j,tr[l][2], creal(tmp), ediff_x, gx_var, gx_mfac, exp(-(powerl(ediff_x, 2)) / gx_var)); */
            /* printf("y %d %d:%le %le %le %le %le\n ",j,k, creal(tmp), ediff_y, gy_var, gy_mfac, exp(-(powerl(ediff_y, 2)) / gy_var)); */
          }
          /* fprintf(stderr, "\n\n=======Valgrind eject point=======\n\n"); */
          /* exit(1); */
          /* tmp = fabsc(tmp); */
          /* tmp *= tmp; */
          /* sm[j][k] += creal(tmp); */
        }

        if (y == spec -> n_ely-1) {
          /* we have traversed one row in the spectrum */
          /* printf("%d %d %d\n", j*spec->prsz+k, y, spec -> n_ely); */
          x++;
          omega_x = emin_x + (x * de_x);
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
    /* printf("finished calculating: thread %d, start %d, end %d\n", ith, l_st, l_fn); */

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
      }
    } /* implicit barrier */
  }
  /* printf(" done (%s).", get_loctime(ltime)); */
  fflush(stdout);
  para_t = omp_get_wtime() - wt;

  if (inp -> md -> lorz) {

    /* printf("\n      applying lorentzian boadening.. (%s)",get_loctime(ltime)); */
    /* fflush(stdout); */

    for (j = 0, x = 0, y = 0; x < spec -> n_elx; j++) {
      for (k = 0; (k < spec -> prsz) && (x < spec -> n_elx); k++, y++) {
        omega_x = emin_x + (x * de_x);
        omega_y = emin_y + (y * de_y);

        /* update the broadenings, if needed */
        if (omega_x > eu_lx) {
          /* printf("thread %d at i= %d set lx since %le > %le\n", ith,lx_i, omega_x, eu_lx); */
          /* fflush(stdout); */
          lx_i += 2;
          lx_t = md -> lx[lx_i] / AUTOEV;

          lx_mfac = 0.5 * lx_t / PI;
          lx_afac = (0.25 * lx_t * lx_t);

          eu_lx = (md -> lx)[lx_i+1] / AUTOEV;
        }
        if (omega_y > eu_ly) {
          ly_i += 2;
          ly_t = md -> ly[ly_i] / AUTOEV;

          ly_mfac = 0.5 * ly_t / PI;
          ly_afac = (0.25 * ly_t * ly_t);

          eu_ly = (md -> ly)[ly_i+1] / AUTOEV;
        }

        tmp_int = 0;
        for (j_in = 0, x_in = 0, y_in = 0; x_in < spec -> n_elx; j_in++) {
          for (k_in = 0; (k_in < spec -> prsz) && (x_in < spec -> n_elx)
                 ;k_in++, y_in++) {
            omega_x_in = emin_x + (x_in * de_x);
            omega_y_in = emin_y + (y_in * de_y);

            /* calculate the lorentzian contribution to this intensity */
            ediff_x = omega_x - omega_x_in;
            ediff_y = omega_y - omega_y_in;

            tmp_int += sm_th0[j_in][k_in]
              * (lx_mfac / ((ediff_x * ediff_x) + lx_afac))
              * (ly_mfac / ((ediff_y * ediff_y) + ly_afac));

            if (y_in == spec -> n_ely-1) {
              x_in++;
              y_in = 0;
            }
          }
        }

        spec -> s_mat[j][k] += tmp_int;

        if (y == spec -> n_ely-1) {
          /* we have traversed one row in the spectrum */
          x++;
          y = 0;
        }
      }
    }
    /* printf(" done (%s).", get_loctime(ltime)); */
  }
  else {
    spec -> s_mat = sm_th0;
  }


  /* printf("\n\n==== THREAD INFO END ==== \n\n"); */
  /* fflush(stdout); */
  /* for (j = 0; j < spec -> npr_tot; j++) { */
  /*   free(sm_th0[j]); */
  /* } */
  /* free(sm_th0); */
  free(ltime);
  free(rchunks);
  return 0;
}

int
calc_spec (struct inp_node *inp, int spec_idx)
{
  /* int x, y; */
  char *ltime = malloc(20);

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
      printf(" a constructive interference model (%s) ..", get_loctime(ltime));
      fflush(stdout);
      intf_0(inp, spec, md);
      break;

    default:
      printf(" no interference model (%s) ..", get_loctime(ltime));
      fflush(stdout);
      intf_0_old(inp, spec, md);
    }

  printf("    done (%s).\n", get_loctime(ltime));

  for (j = 0; j < spec -> npr_tot; j++) {
    for (k = 0; k < spec -> prsz; k++) {
      if (spec -> s_mat[j][k] > spec -> sfac) {
        spec -> sfac = spec -> s_mat[j][k];
      }
    }
  }

  free(ltime);
  return 0;
}
