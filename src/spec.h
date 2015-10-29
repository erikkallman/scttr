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
#ifndef SPEC_H
#define SPEC_H

struct spec_s;
typedef struct spec_s * spec;

struct spec_s{

  int idx;
  int layer;
  int n_layers;
  int height;
  int length;

  double ** sdat; /* the actual spec data */

  /* next spec for the information node */
  spec next;
  spec last;

  /* next layer of this spec */
  spec next_l;
  spec last_l;
  spec root_l;
};


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
free_spec_stack (sctr_input inode,
                 spec root_s,
                 int idx);

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
init_spec (sctr_input inode,
           double ** s_data,
           int s_idx,
           int ly,
           int h,
           int l
           );

#endif /* SPEC_H */
