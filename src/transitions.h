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
#ifndef TRS_C
#define TRS_C
#include "inp_node_s.h"

int
get_erange (struct inp_node *inp, double e);

int
get_i (struct inp_node *inp, int from);

int
get_il (struct inp_node *inp, int from);

int
get_trs (int from, double **trs);

int
get_inext (struct inp_node *inp, int from);

int
get_ilnext (struct inp_node *inp, int from);

int
get_trsnext (double **trs, int from);

void
count_states (struct inp_node *inp);

int
add_sym (struct inp_node *inp);

#endif /* TRS_C */
