#ifndef SMAP_H
#define SMAP_H
#include "dynarray.h"
/* function calc_smap

   * synopsis:

   * algorithm:

   * input:

   * output:

   * side-effects:

   */
void
calc_smap_m (char * fn_relpath,
             char * fn_infile,
             int len_infile,
             double * state_er,
             double * state_t,
             double * res
             );

void
calc_smap_dbg(double * state_er,
             double * state_t
             );

void
write_log (double * state_er,
           double * state_t,
           double * res,
           char * fn_relpath,
           char * fn_infile,
           int len_infile,
           int n_max
           );

#endif /* SMAP_H */
