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
   * @file transitions.h
   * @author Erik Källman
   * @date November 2015
   * @brief This file defines the interface to all functions defined for
   * the @p trs variable in the input node struct (see @p inp_node).
   */
#ifndef TRS_C
#define TRS_C
#include "inp_node_s.h"

/**
 * @brief provided with the energy value @p e this function determines which of
 * the energy ranges (see the metadata struct), provided by the user, that @p e
 * belongs to.
 *
 * @param inp the input node that contains the energy ranges
 * @param e the energy value for which the energy range is to be determined
 * @returns the integer j, defining which of the energy ranges @p e is inside of
 * , or -1 if @p e is not inside any of the user-defined energy ranges.
   */
int
get_erange (struct inp_node *inp, double e);

/**
   * @brief This function searches the @p trs matrix for the index of the
   * @p from state and provides the user with its index inside the trs matrix.
   *
   * @param from the index of the state indentifying a specific energy eigenstate
   * for a given transition from state x to y
   * @param trs the transition matrix inside which the @p from index is to
   * be searched for.
   * @returns the integer @p j if the @p from state is found, otherwise -1
   */
int
get_trs (int from, double **trs);

/**
   * @brief The @p trs matrix contains blocks of transitions from a state x
   * to multiple states y1, y2,...yn. To skip ahead to the next block of
   * transitions, this function looks up the index of the first transition
   * in a new block of transitions.
   *
   * @param trs the @p trs matrix in which to look for the state of index
   * @p from.
   * @param from the index of the eigenstate from which to start searching
   * for the next block of transitions.
   * @returns the index @p j of the position of the element in @p trs that
   * starts the new block of transitions.
   */
int
get_inext (double **trs, int from);

/**
   * @brief In order to add elastic transitions to the trs matrix, this
   * functions allows the program to figure out how many states fit into
   * each energy range provided by the user (see main.c, @p state_er,
   * also the metadata struct).
   *
   * @param inp the input node that contains the @p trs matrix that will
   * be searched for states inside of the energy ranges provided by the user.
   * @returns void
   * @note side effects: exits with EXIT_FAILURE upon failed memory allocation.
   */
void
count_states (struct inp_node *inp);

/**
   * @brief Specifically for input data derived from Molcas calculations,
   * elastic transitions (from state x to y and then from y back to x) are not
   * present. This functions adds those transitions to a give @p trs matrix.
   *
   * If the energy range provided by the user specifies that all
   * transitions are from ground states and to final states in the same energy
   * range this function is used to add elastic transitions to the @p trs matrix
   * of the provided input node struct @p inp.
   *
   * @param inp the input node whos @p trs matrix should be expanded to hold
   * elastic transitions.
   * @returns 1 if successful
   * @note side effects: exits with EXIT_FAILURE if memory allocation fails.
   * allocates enough memory to the @p trs matrix to acommodate all possible
   * elastic transitions.
   */
int
add_eltrans (struct inp_node *inp);

#endif /* TRS_C */
