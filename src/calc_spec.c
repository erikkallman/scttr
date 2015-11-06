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
#include "calc_spec.h"
#include "dyn_array.h"
#include "scttr_io.h"
#include "sci_const.h"
#include "std_num_ops.h"
#include "std_char_ops.h"
#include "spectrum_s.h"
#include "spectrum.h"
#include "inp_node_s.h"

int
calc_spec (struct inp_node *inp)
{

  struct metadata *md = inp -> md;
  struct spectrum *spec = get_spec(inp, 2);

  int j, k, x, y;
  int is_idx, is_pos;

  /* create local copies and pointers to variables whos name would
     otherwise clutter the code loop */
  int n_sfs = spec -> is2fs -> n_el;
  int *is2fs = spec -> is2fs -> a;
  int *gs2is = spec -> gs2is -> a;
  int *is_idxs = spec -> is_idxs -> a;
  int *ii_start = spec -> ii_start -> a;

  double de_gi, de_if; /* energy eigenvalue differences */
  double o_x, o_y; /* energies at each point in the spectrum matrix  */
  double bw; /* boltzmann weight */

  /* energies and transition moments along the axis .. */
  double ediff_x, tmom_gi; /* .. of excitiation (x-axis)*/
  double ediff_y, tmom_if; /* .. of energy transfer (y-axis)*/
  double emin_x = spec -> emin_x;
  double emax_x = spec -> emax_x;
  double emin_y = spec -> emin_y;
  double emax_y = spec -> emax_y;
  double de_x, de_y; /* element to element energy difference in the s_mat
                        matrix */
  double fwhm_x, fwhm_y; /* full-width half-maxima */
  double complex tmp;  /* accumulator used in the Kramers-Heisenberg formula */

  /* excitation energy (x-axis) broadening parameters */
  double grms_x; /* gaussian RMS value */
  double gvar_x; /* square root of the variance of the gaussian */

  /* transfer energy (y-axis) broadening parameters */
  double grms_y; /* gaussian RMS value */
  double gvar_y; /* square root of the variance of the gaussian */

  fwhm_x = md -> fwhm[0] / AUTOEV;
  fwhm_y = md -> fwhm[1] / AUTOEV;

  de_x = md -> res[0] / AUTOEV;
  de_y = md -> res[1] / AUTOEV;

  spec -> n_elx = (int)((emax_x - emin_x) / de_x);
  spec -> n_ely = (int)((emax_y - emin_y) / de_y);

  if((spec -> omega_x = malloc(spec -> n_elx * sizeof(double *))) == NULL ) {
    fprintf(stderr, "calc_spec.c:function calc_spec, malloc: failed to allocate memory for \"spec -> omega_x\"\n");
    printf("program terminating due to the previous error.\n");
    exit(1);
  }

  if((spec -> omega_y = malloc(spec -> n_elx * sizeof(double *))) == NULL ) {
    fprintf(stderr, "calc_spec.c:function calc_spec, malloc: failed to allocate memory for \"spec -> omega_y\"\n");
    printf("program terminating due to the previous error.\n");
    exit(1);
  }

  if((spec -> s_mat = malloc(spec -> n_elx * sizeof(double *))) == NULL ) {
    fprintf(stderr, "calc_spec.c:function calc_spec, malloc: failed to allocate memory for \"s_mat\"\n");
    printf("program terminating due to the previous error.\n");
    exit(1);
  }

  for (j=0; j<spec -> n_elx; j++) {
    if((spec -> s_mat[j] = malloc(spec -> n_ely * sizeof(double))) == NULL ) {
      fprintf(stderr, "calc_spec.c:function calc_spec, malloc: failed to allocate memory for \"spec -> s_mat[%d]\"\n"
              , j);
      printf("program terminating due to the previous error.\n");
      exit(1);
    }

    if((spec -> omega_x[j] = malloc(spec -> n_ely * sizeof(double))) == NULL ) {
      fprintf(stderr, "calc_spec.c:function calc_spec, malloc: failed to allocate memory for \"spec -> omega_x[%d]\"\n"
              , j);
      printf("program terminating due to the previous error.\n");
      exit(1);
    }

    if((spec -> omega_y[j] = malloc(spec -> n_ely * sizeof(double))) == NULL ) {
      fprintf(stderr, "calc_spec.c:function calc_spec, malloc: failed to allocate memory for \"spec -> omega_y[%d]\"\n"
              , j);
      printf("program terminating due to the previous error.\n");
      exit(1);
    }
  }

  grms_x = 2.0 * powerl((fwhm_x / (2 * sqrt(2 * log(2)))), 2);
  gvar_x = fwhm_x / (2.0 * sqrt(2.0 * log(2))) * sqrt(2.0 * 3.1415927);

  grms_y = 2.0 * powerl((fwhm_y / (2 * sqrt(2 * log(2)))), 2);
  gvar_y = fwhm_y / (2.0 * sqrt(2.0 * log(2))) * sqrt(2.0 * 3.1415927);

  printf("  - calculating the RIXS map .. \n");

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

          bw = get_rbdist(inp -> e0,inp -> trs[2][gs2is[is_idx]]);

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

  printf("    .. done.\n");
  printf("  - finding the normalization constant of the  calculated RIXS map ..");
  for (x = 0; x < spec -> n_elx; x++) {
    for (y = 0; y < spec -> n_ely; y++) {
      if (spec -> s_mat[x][y] > spec -> sfac) {
        spec -> sfac = spec -> s_mat[x][y];
      }
    }
  }

  printf(" done.\n");
  spec -> sfac = spec -> sfac;
  return EXIT_SUCCESS;
}
