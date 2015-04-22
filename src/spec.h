#ifndef SPEC_H
#define SPEC_H
#include "structs.h"

/* function get_spec
   Searches the list of specta connected to the root spec @root_s of a given
   info node for a node of index == @idx, then returns that node. if the
   requested spec is not in the list, the last spec in the list is
   returned instead.

   * synopsis:

   * algorithm:

   * input:

   * output:

   * side-effects:

   */
spec
get_spec (spec root_s,
          int idx,
          int layer);

int
free_spec (spec s,
           int idx,
           int layer);

void
free_spec_stack (int idx);

spec
get_last_layer (spec root_s );

/* function append_spec_layer
   append spec @s with index @s->idx to the spec of the same index
   as a layer. called by broadening functions (yet to be implemented).
   * synopsis:

   * algorithm:

   * input:

   * output:

   * side-effects:

   */
void
append_spec_layer (spec s,
                   spec root_s);

/* function append_spec

   * synopsis:
   appends the spec @s of index @idx to the list of spectra connected to the
   root spec @root_s of an info node.

   * algorithm:

   * input:

   * output:

   * side-effects:

   */

void
append_spec (spec s,
             spec root_s
             );


void
init_spec (info_node inode,
           double ** s_data,
           int s_idx,
           int ly,
           int h,
           int l
           );

#endif /* SPEC_H */
