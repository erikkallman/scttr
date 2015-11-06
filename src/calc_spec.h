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
   * @file calc_spec.h
   * @author Erik Källman
   * @date November 2015
   * @brief This file contains the interface for the calc_spec() function,
   * which calculates a spectrum, stored in the s_mat variable in a
   * spectrum node.
   */
#ifndef CALC_SPEC_H
#define CALC_SPEC_H
#include "inp_node_s.h"
/**
   * @brief This function uses the data contained in a given spectrum attached
   * to the @p inp provided as argument, and the @p trs variable to define a
   * matrix containing the calculated spectrum.
   *
   * To calculate the spectrum, stored in the @p s_mat variable of the spectrum
   * struct, the reduced Kramers-Heisenberg formula (M. Lundberg et. al. 2013)
   * is used.
   * Suitable functions are used to broaden the intensities summed in
   * each matrix element of the @p s_mat matrix.
   * Variables with either and @p _x or @p _y suffix refers to properties of
   * the excitiation and energy transfer axis respectively.
   *
   * @param inp the input node containing the data from which a spectrum matrix will be calculated
   * @returns EXIT_SUCCESS if successful, otherwise EXIT_FAILURE
   * @note side-effects: defines and initializes the following variables of a
   * spectrum struct: @p omega_x, @p omega_y, @p s_mat, @p sfac.
   * See the spectrum struct for more information.
   * @note the calculated spectrum matrix is not normalized in this function.
   * This is done in the write_spec() function.
   */
int
calc_spec (struct inp_node *inp);

#endif /* CALC_SPEC_H */
