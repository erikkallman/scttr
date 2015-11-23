/* This file is part of Scatter. */

/* Scatter is free software: you can redistribute it and/or modify */
/* it under the terms of the GNU Lesser General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or */
/* (at your option) any later version. */

/* Scatter is distributed in the hope that it will be useful, */
/* but without any warranty; without even the implied warranty of */
/* merchantability or fitness for a particular purpose. See the */
/* GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License */
/* along with Scatter, found in the "license" subdirectory of the root */
/* directory of the Scatter program. If not, see <http://www.gnu.org/licenses/>. */
/**
   * @file inp_node_s.h
   * @author Erik Källman
   * @date November 2015
   * @brief This file contains the implementation of the inp_node (input node)
   * struct.
   */
#ifndef SCTTR_INPUT_S_H
#define SCTTR_INPUT_S_H
#include "metadata_s.h"
#include "spectrum_s.h"

int n_inp; /* total number of input nodes */
struct inp_node *root_inp;

/**
 * @brief Each input file loaded into the program has an associated
 * node containing all data about, and extracted from, that
 * specific input file. That node is stored in a inp_node_s struct.
 *
 * In the global scope, root_inp provides an entry point for any function
 * to search the linked list of input nodes for a specific node,
 * corresponding to a specific input file read into the data stored in
 * and used by the program. In addition to basic datatypes, each input node
 * also contains @p metadata defined from the input, as well as a linked list
 * of spectra, derived from the screening and energy-range parameters
 * provided by the user.
 */
struct inp_node
{

  int idx; /**< index of node in the list value */

  int n_states; /**< number of electronic states */
  int n_trans; /**< number of transitions between those states */
  int n_spec; /**< number of spectra contained in this node */
  int n_tmax; /**< maximum number of intermediate state transitions from
               any given electronic state */

  int n_gfs; /**< number of ground and final electronic states */
  int n_is; /**< number of intermediate electronic states */

  int *idx_map; /**< a mapping between a state of a given index, to the position
                  of its transitions stored in the trs variable (see below) */

  /**< sum of the boltzmann weights of all states in the system*/
  double bw_sum; /**< sum of boltzmann weights  or all ground states */
  double tmax_q,tmax_d; /**< maximum transition moment of all quadrupole (q) and
                           dipole (tmax_d) transitions */
  double e0; /**< lowest energy eigenstate */

  double **trs;  /**< the trs variable contains an array of six rows and n_trans
                    columns. for a transition from state x to state yeach
                    row stores the following data:
                    [0]: index of state x
                    [1]: index of state y
                    [2]: energy eigenvalue of state x
                    [3]: energy eigenvalue of state y
                    [4]: transition moment for the transition from x to y
                    [5]: type of transition, 1=dipole, 2=quadrupole.
                 */

  struct metadata *md; /**< spectrum information node metadata */
  struct spectrum *root_spec; /**< spectrum derived from the input */

  struct inp_node *next;
  struct inp_node *last;
};

#endif /* SCTTR_INPUT_S_H */
