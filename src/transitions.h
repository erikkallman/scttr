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
#ifndef TRS_C
#define TRS_C
#include "sctr_input.h"

int
get_erange (sctr_input s_inp,
            double e);

int
get_i (sctr_input s_inp,
       int from);

int
get_il (sctr_input s_inp,
        int from);

int
get_trs (int from,
         double ** trs);

int
get_inext (sctr_input s_inp,
           int from);

int
get_ilnext (sctr_input s_inp,
            int from);

int
get_trsnext (double ** trs,
             int from);

int
eval_trs (sctr_input s_inp);

void
count_states (sctr_input s_inp);

int
add_sym (sctr_input s_inp);

#endif /* TRS_C */
