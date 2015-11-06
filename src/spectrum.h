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

#ifndef SPECTRUM_H
#define SPECTRUM_H
#include "spectrum_s.h"
#include "inp_node_s.h"

struct spectrum *
get_spec (struct inp_node *inp, int idx);

struct spectrum *
init_spec (struct inp_node *inp, int cap, int inc);

int
set_root_spec (struct inp_node *inp);

/**
 * @brief The set_spec() function adds a spectrum to the list currently on @p inp, based on its root spectrum (see the inp_node struct)
 *
 * In addition to including the boltzmann weight screening results already
 * in the root spectrum of @p inp, set_spec() also uses the total intensity
 * threshold value provided by the user. If the user requests 98% of the total
 * intensity to be kept after screening, set_spec() will screen out, starting
 * from the least intense, transitions until that threshold is reached.
 * @param inp the input node to which a spectrum will set, and appended
 */
int
set_spec (struct inp_node *inp);

int
add_spec (struct inp_node *inp, struct spectrum * spec);

int
free_spec (struct spectrum * spec);

int
free_all_specs (struct inp_node *inp);

#endif /* SPECTRUM_H */
