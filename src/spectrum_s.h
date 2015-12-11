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
   * @file spectrum_s.h
   * @author Erik Källman
   * @date November 2015
   * @brief This file contains the implementation of the spectrum struct.
   */
#ifndef SPECTRUM_S_H
#define SPECTRUM_S_H
#include "dyn_array.h"

/**
   * @brief For a given energy range and screening values provided by the user
   * the spectrum struct stores the data resulting from applying that user input
   *  to the @p trs variable in a given inp_node struct
   *
   * The spectrum struct stores a map of the only transitions in the @p trs
   * variable that is accepted after screening. This data is later used
   * in calc_spec() to avoid branching of the most instruction-dense
   * loop structure, as well as not having to iterate over the @p trs matrix in
   * order to find the transitions that survived screening.
   *
   *  As a result of a call to the set_spec function, two dynamic arrays
   * (see the da struct) are defined, that each respectively contain indices
   * of elements in the @p trs matrix of a screened transition from intermediate
   * to a final state (is2fs) and from a ground state to an intermediate state
   * (gs2is).
   *
   * Two arrays, @p is_idxs and @p ii_start map elements from is2fs to gs2is.
   * For a given intermediate (I) to final state transition in is2fs,
   * transitions from a given ground state to that intermediate state I are
   * located in the gs2is matrix, starting from the index is_idxs[j], where
   * j is the index of the corresponding final state transition in is2fs.
   * is_idxs is a mapping between transitions in is2fs and gs2is which is used
   * to avoid having to search through gs2is for the right intermediate state index.
   * Essentially, the four arrays form a set of indexed linked lists.
   *
   * @note The default values of the @p emin and @p emax variables, used in
   * init_spec() are defined from the energy ranges provided by command
   * line input. The variables contained in this struct are defined and
   * initialized at various places in the code, and are always mentioned as
   * notes to functions where they are.
   */
struct spectrum
{

  int idx; /**< index of this spectrum in the linked list of spectra */
  int n_st; /**< number of transitions kept after screening */
  int n_elx, n_ely; /**< Number of elements in x and y (row and column) direction of the calculated spectrum matrix. */

  /* variables used to initialize the partitioned spectrum matrix, passed to
     each thread in the calculation. */
  int prsz; /**< Partitioned row size in number of elements */
  int npr; /**< Number of partitioned rows for each row in the spectrum matrix.  */
  int npr_tot; /**< Total number of partitioned rows */

  double emin_x, emax_x; /**< maximum and minimum energies in the range of
                            final states */
  double emin_y, emax_y; /**< maximum and minimum energies in the range of
                            intermediate states. */

  double sfac; /**< scaling factor used to normalized the spectrum matrix */

  double **omega_x; /**< energy range in x dimension */
  double **omega_y; /**< energy range in x dimension */
  double **s_mat; /**< the matrix containing the calculated spectrum,
                   resulting from calling the calc_spec() function */

  double ** trs_red; /**< The transition matrix reduced to the optimal size for the current cache architecture. */

  da is2fs; /**< index of a transition from intermediate to final state */
  da gs2is; /**< index of a transition from ground to intermediate state */
  da is_idxs; /**< intermediate state transition indices in gs2is */
  da ii_start; /**< intermediate state transition index start in gs2is */

  struct spectrum *last_spec; /**< pointer to the "left" spectrum in the list */
  struct spectrum *next_spec; /**< pointer to the "right" spectrum in the list */

};

#endif /* SPECTRUM_S_H */
