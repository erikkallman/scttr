#ifndef TRS_C
#define TRS_C
#include "spectrum_info.h"

int
get_erange (spec_info s,
            double e);

int
get_i (spec_info s,
       int from);

int
get_il (spec_info s,
        int from);

int
get_trs (int from,
         double ** trs);

int
get_inext (spec_info s,
           int from);

int
get_ilnext (spec_info s,
            int from);

int
get_trsnext (double ** trs,
             int from);

int
eval_trs (spec_info s);

int
eval_trs_mult (spec_info s);

void
count_states (spec_info s);

int
add_sym (spec_info s);

#endif /* TRS_C */
