#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "dyn_array.h"

int
da_extend (da src) {

  int j;
  int * a_dest = malloc(sizeof(int)*(src->cap+src->inc));

  for (j=0; j < src->n_el; j++) {
    a_dest[j] = src->a[j];
  }

  free(src->a);
  src-> a = a_dest;
  src-> cap = src->cap+src->inc;

  return(EXIT_SUCCESS);
}

int
da_shrink (da src) {

  int j;
  src -> inc = src-> n_el/2;
  src -> cap = src->n_el + src -> inc;
  int * a_dest = malloc(sizeof(int)*(src->cap));

  for (j=0; j < src->n_el; j++) {
    a_dest[j] = src->a[j];
  }

  free(src->a);
  src-> a = a_dest;

  return(EXIT_SUCCESS);
}

/* delete element by shifting all below it up one step */
int
da_del_us (da ar,
           int el) {
  int j;
  int val;
  int n_el = ar -> n_el;

  for (j=el; j<n_el-2; j++) {
    val =  ar -> a[j+1];
    ar -> a[j] = val;
  }

  if ( (--ar->n_el)/ar->cap <= 0.25) {
    da_shrink(ar);
  }

  ar->n_el--;

  return EXIT_SUCCESS;
}

int
da_set (da ar,
        int el,
        int val) {

  if (el >= ar->cap ) {
    da_extend(ar);
  }

  ar -> a[el] = val;

  return EXIT_SUCCESS;
}

int
da_append (da ar,
           int val) {

  if (ar->n_el >= ar->cap) {
    da_extend(ar);
  }

  ar -> a[ar->n_el] = val;
  ar -> n_el++;

  return EXIT_SUCCESS;
}

int
da_get (da ar,
        int el) {

  if (el > ar -> cap) {
    fprintf(stderr, "\n\ndyn_array.c, function da_get: trying to get value if index outside of array range.\n");
    printf( "program terminating due to the previous error.\n");
    exit(EXIT_FAILURE);
  }
  return ar -> a[el];
}

int
da_getlast (da ar) {

  if (ar->n_el == 0) {
    fprintf(stderr, "\n\ndyn_array.c, function da_getlast: the array in a is empty.\n");
    printf( "program terminating due to the previous error.\n");
    exit(EXIT_FAILURE);
  }

  return ar -> a[(int)(ar -> n_el)-1];
}

da
da_init (int cap,
         int inc) {

  int * a;
  da new_da;

  if((new_da = malloc(sizeof(struct da_s))) == NULL ){
    fprintf(stderr, "da.c:function da_init, malloc: failed \
to allocate memory for \"new_da\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  if((a = malloc(sizeof(int)*cap)) == NULL ){
    fprintf(stderr, "da.c:function da_init, malloc: failed \
to allocate memory for \"a\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  if (cap == 0) {
    cap = DEFAULT_CAP;
  }
  if (inc == 0) {
    inc = DEFAULT_INC;
  }

  new_da -> n_el = 0;
  new_da -> cap = cap;
  new_da -> inc = inc;
  new_da -> a = a;

  return new_da;
}
