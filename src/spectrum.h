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
   * @file spectrum.h
   * @author Erik Källman
   * @date November 2015
   * @brief This file contains the interface to functions defined for the spectrum struct.
   */
#ifndef SPECTRUM_H
#define SPECTRUM_H
#include "spectrum_s.h"
#include "inp_node_s.h"

/**
 * @brief For a given index number @p idx, this function searches the linked
 * list of spectra attached to the input node @p inp for a spectrum with the
 * corresponding index, and then returns it.
 *
 * @param inp the input node inside which a spectrum of index @p idx will
 * be searched for
 * @param idx the integer number defining the index of what spectrum to
 * search for inside the input node.
 * @returns the spectrum struct @p spec if it is found in the list.
 * @note upon failed memory allocation the program terminates prematurily
 * with the EXIT_FAILURE macro.
 *
 */
struct spectrum *
get_spec (struct inp_node *inp, int idx);

/**
   * @brief The init_spec() function initializes a spectrum struct and returns
   * it to the callee.
   *
   * @param inp the input node from which the
   * @param cap the initial capacity (in terms of elements) of the arrays
   * contained in the spectrum struct.
   * @param inc the initial increment (in terms of elements) of the arrays
   * contained in the spectrum struct (see the dyn_array struct for more
   * information on how this variable is used).
   * @returns spec the newly initialized spectrum struct.
   * @note this function does not attach the newly initialized spectrum to the
   *  input node @p inp, provided in the argument list.
   */
struct spectrum *
init_spec (struct inp_node *inp, int cap, int inc);

/**
   * @brief Using the boltzmann weight of the energies the ground states
   * found in the @p trs matrix of the input node @p inp, this function
   * screens out all ground states below the user-defined threshold. It uses
   * all transitions from states that were kept after screening, along with the
   * energy ranges provided by the user, to initialize the dynamic arrays
   * of @p root_spec (see struct spectrum for more details).
   *
   * @param inp the input node to which the @root_spec spectrum will be added
   * @returns 1 if successfully executed, 0 if screening failed.
   * @note side effects: initializes the @p root_spec variable with init_spec(),
   * then appends the spectrum @p root_spec to @p inp. exits with EXIT_FAILURE
   * upon failed malloc call.
   */
int
set_root_spec (struct inp_node *inp);

int
set_root_spec_el (struct inp_node *inp);

/**
 * @brief The set_spec() function adds a spectrum to the list currently on
 * @p inp, based on its root spectrum (see the inp_node struct)
 *
 * In addition to including the boltzmann weight screening results already
 * in the root spectrum of @p inp, set_spec() also uses the total intensity
 * threshold value provided by the user, to sceen out transitions not
 * intense enough to be included. If the user requests 98% of the total
 * intensity to be kept after screening, set_spec() will screen out, starting
 * from the least intense, transitions until that threshold is reached.
 * @param inp the input node for which a spectrum will defined, and then
 * added
 * @returns 1 if successful.
 * @note side effects: initializes a spectrum struct (see init_spec()),
 * defines values in that struct and also appends the struct to the input node
 * @p inp provided as input argument. exit with EXIT_FAILURE upon failed
 * malloc call.
 */
int
set_spec (struct inp_node *inp);

/**
   * @brief if the input argument @p spec is a previously initialized spectrum
   * struct, this function adds that @p spec to the input node @p inp.
   * @param inp the input node to which the callee wants to add a spectrum
   * @param spec the spectrum to be appended to the input node @p inp
   * @returns 1 if successful, otherwise, if the @p spec struct is
   * not initialized, exits with EXIT_FAILURE.
   * @note side effects: sets the @p idx variable of @p spec
   */
int
add_spec (struct inp_node *inp, struct spectrum * spec);

/**
   * @brief In order to form the matrix of transitions, which is used in the
   * calc_spec() function to generate the spectrum, this function iterates over
   * the indexing arrays of the provided spectrum @p spec.
   *
   * @param inp the input node containing the spectrum of index @p spec_idx
   *
   * @returns 1 if successful.
   * @note side effects: exits with @p EXIT_FAILURE upon failed @p malloc call
   * , and if the spectrum has not had its indexing arrays defined.
   */
int
set_trs_red (struct inp_node *inp, int spec_idx);

/**
   * @brief Frees up the memory previously allocated for a given spectrum
   * struct, while keeping the linked list structure intact.
   * @param spec the spectrum struct to have its memory freed.
   * @returns 1 if successful.
   * @note side effects: exits with @p EXIT_FAILURE upon failed @p malloc call.
   */
int
free_spec (struct spectrum * spec);


/**
   * @brief Frees up the memory for all previously added spectra for a given input
   * input node struct, while keeping the linked list structure intact.
   * @param inp the input node struct to have the memory for its entire linked
   * list of spectrum structs freed.
   * @returns 1 if successful.
   * @note exits with @p EXIT_FAILURE upon failed malloc call.
   */
int
free_all_specs (struct inp_node *inp);

void
trs2str (struct spectrum *spec);

void
strs2str (struct inp_node *inp, struct spectrum *spec);

void
tot_strs2str (struct inp_node *inp, struct spectrum *spec);

#endif /* SPECTRUM_H */
