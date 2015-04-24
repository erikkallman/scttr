#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "spec.h"
#include "structs.h"

spec
get_spec (spec root_s,
          int idx,
          int layer) {

  spec next_spec;
  spec curr_spec = root_s;

  while((curr_spec == NULL) && ((curr_spec -> idx) != idx)){
    next_spec = curr_spec -> next;
    curr_spec = next_spec;
/*     if (curr_spec == NULL ) { */
/*       fprintf(stderr, "\n\nERROR: spec.c, function get_spec, tried to locate \ */
/* layer %d of spec %d, but failed. spec not found in list.\n",idx,layer); */
/*       printf( "program terminating due to the previous error.\n\n"); */
/*       exit(1); */
/*     } */
  }

  if (layer > 1) { /* locate the right layer of the spec */
    curr_spec = curr_spec -> next_l;
    while((curr_spec -> layer) != layer){
      next_spec = curr_spec -> next_l;
      curr_spec = next_spec;
    }
  }

  return curr_spec;
}

int
free_spec (spec s,
           int idx,
           int layer){
  int j;

  int l = s -> length;
  int h = s -> height;
  /* printf( "%le %d %d\n", (s->sdat)[50][50], l, h); */
  /* fprintf(stderr, "\n\n=======Valgrind eject point=======\n\n"); */
  /* exit(1); */
  for (j=0; j<h; j++) {
    free((s->sdat)[j]);
  }
  free(s->sdat);

  if ((layer == 1) && (s -> last != NULL)) { /* connect the adjecent spectra */
    s -> last -> next = s -> next;
    s -> next -> last = s -> last;
  }
  else if((s -> last == NULL) && (s -> next == NULL)){
    /* the last of all spectra */

  }
  else { /* connect the adjecent layers */
    s -> last_l -> next_l = s -> next_l;
    s -> next_l -> last_l = s -> last_l;
  }

  free(s);

}

void
free_spec_stack (info_node inode,
                 spec root_s,
                 int idx){

  int j;
  int n_l; /* number of specta */

  spec s_last;

  /* get the root spec of the stack on idex idx */
  spec s = get_spec(root_s, idx, 1);

  n_l = s -> n_layers;

  if (s -> n_layers > 1) {
    s = get_spec(s, idx, n_l);
  }

  for (j=n_l; j>0; j--) {
    s_last = s -> last_l;
    free_spec(s, idx, j);
    inode -> n_spec--;
    s = s_last;
  }
}

spec
get_last_layer (spec root_s) {

  spec next_spec = root_s;
  spec curr_spec;

  while(next_spec != NULL){
    curr_spec = next_spec;
    next_spec = curr_spec -> next_l;
  }

  return curr_spec;
}


void
append_spec_layer (spec s,
                   spec root_s
                       ) {

  /* spec to append to */
  spec s_a = get_spec(root_s, s->idx, root_s ->  n_layers);

  /* the last layer appended to that spec */
  spec s_l;
  s -> root_l = root_s;

  if ((s_a -> idx) != (s -> idx)) {
    fprintf(stderr, "spec.c, function append_spec_layer: tried to add \
spec of idx %d as a layer to %d, but the zero-layer spec for that index was\
 not found in the list. \n", s_a -> idx, s -> idx);
    printf( "program terminating due to the previous error.\n");
    exit(1);
  } else {
    s_a -> root_l -> n_layers++;
    s_l = get_last_layer(s_a);
    s_l -> next_l = s;
    s -> last = s_l;
    s -> next = NULL;
  }
}

void
append_spec (spec s,
             spec root_s
                 ) {

  /* spec to append to */
  spec s_a = get_spec(root_s, s->idx,0);

  s_a -> next = s;
  s -> last = s_a;
  s -> next = NULL;

}

void
init_spec (info_node inode,
           double ** s_data,
           int s_idx,
           int ly,
           int h,
           int l
           ){

  int j;

  spec new_spec;

  if((new_spec = malloc(sizeof(struct spec_s))) == NULL ){
    fprintf(stderr, "spec.c:function init_spec, malloc: failed \
to allocate memory for \"new_spec\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  inode -> n_spec++;
  new_spec -> idx = inode -> n_spec;
  new_spec -> layer = ly;
  new_spec -> n_layers = 1;
  new_spec -> height = h;
  new_spec -> length = l;

  new_spec -> sdat = s_data;

  new_spec -> next = NULL;

  if (ly == 1) {
    if (inode -> n_spec == 1) { /* there is no root info node defined  */

      new_spec -> last = NULL;
      inode -> root_spec = new_spec;

      inode -> root_spec -> root_l = new_spec;
      inode -> root_spec -> last_l = NULL;
      inode -> root_spec -> next_l = NULL;

    }
    else{ /* append the new spec to the list */
      append_spec(new_spec, inode->root_spec);
    }
  }
  else{ /* if the function got called with a non-zero layer, append
         the spec as a layer instead */
    inode -> n_spec--;
    new_spec -> idx = s_idx;
    append_spec_layer(new_spec, inode->root_spec);
  }
}

void
specs2s (info_node inode) {

  printf( "\n  -printing the content of the spec list for the inode \
to file %s:\n", inode -> str_id);

  spec rs, cs, ns, cs_l, ns_l;

  if (inode->n_spec > 0) {
    rs  = inode -> root_spec;
    ns = rs -> next;
    cs = rs;

    while(cs != NULL){
      printf( "  index %d/%d layer %d\n", cs->idx, cs->n_layers, cs->layer);
      cs_l = cs -> root_l;

      while(cs_l != NULL){
        printf( "    layer %d\n", ns_l -> layer);
        ns_l = cs_l -> next_l;
        cs_l = ns_l;
      }
      ns = cs -> next;
      cs = ns;
    }
  } else {
    printf( "  the list is empty.\n");
  }

}
