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
calc_smap_m (char * method,
           char * inode_id,
           mdda_s * idxs
           );

void
calc_smap_dbg (char * method,
           char * inode_id,
           mdda_s * idxs
           );
#endif /* SMAP_H */
