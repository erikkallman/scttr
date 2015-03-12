#include <stdio.h>
#include <stdlib.h>
#include "std_num_ops.h"
#include "dynarray.h"

void da_init(da_s *da) {
  // initialize size and capacity
  da->size = 1;
  da->capacity = DA_INITIAL_CAPACITY;

  // allocate memory for da->data
  da->data = malloc(sizeof(int) * da->capacity);
}

mdda_s *
mdda_init() {
  int j=0; /* looping variables */
  mdda_s * root_mdda;
  mdda_s * next_mdda;
  mdda_s * curr_mdda;

  /* create the mdda chain root */
  root_mdda = malloc(sizeof(mdda_s));
  root_mdda -> idx = 0;
  root_mdda -> size = MDDA_INITIAL_CAPACITY ; /* used to check if the da needs initialization */
  root_mdda -> mdda_capacity = MDDA_INITIAL_CAPACITY;

  root_mdda->da = malloc(sizeof(da_s));
  (root_mdda->da)->size = 1;
  (root_mdda->da)->capacity = DA_INITIAL_CAPACITY;

  // allocate memory for da->data
  (root_mdda->da)->data = malloc(sizeof(int) * ((root_mdda->da)->capacity));

  root_mdda -> next = malloc(sizeof(mdda_s));
  root_mdda -> root = root_mdda;
  root_mdda -> head = root_mdda;
  curr_mdda = root_mdda;

  while(++j < MDDA_INITIAL_CAPACITY){
    /* printf( "initializing %d %d %d \n", j,  (curr_mdda->root) -> idx, root_mdda -> idx); */
    next_mdda = (curr_mdda->next);
    curr_mdda = next_mdda;
    curr_mdda -> idx = j;
    curr_mdda -> size = 1; /* used to check if the da needs initialization */
    curr_mdda -> mdda_capacity = MDDA_INITIAL_CAPACITY;
    /* printf( "da_initing\n" ); */
    curr_mdda->da = malloc(sizeof(da_s));
    (curr_mdda->da)->size = 1;
    (curr_mdda->da)->capacity = DA_INITIAL_CAPACITY;

    // allocate memory for da->data
    (curr_mdda->da)->data = malloc(sizeof(int) * ((curr_mdda->da)->capacity));

    /* da_init(&(curr_mdda->da)); */
    /* printf( "da_inited\n" ); */
    /* printf( "allocating\n" ); */
    curr_mdda -> next = malloc(sizeof(mdda_s));
    /* printf( "allocated\n" ); */
    curr_mdda -> root = root_mdda;
    curr_mdda -> head = root_mdda;

    /* printf( "initialized %d \n\n", curr_mdda -> idx); */
  }
  curr_mdda -> next = NULL;

  /* each mdda node is just a column containing a pointer to the next column as well as a pointer to a da */
  return root_mdda;
}

void da_append(da_s *da, int value) {
  // make sure there's room to expand into
  da_double_capacity_if_full(da);

  // append the value and increment da->size
  da->data[da->size++] = value;
}

void mdda_append(mdda_s *mdda, int val) {

  da_s * da = ((mdda->head)->da);
  da_append(da,val);

  da->size++;
}

int da_get(da_s *da, int index) {
  if (index >= (da->size+1) || index < 0) {
    printf("da_get: Index %d out of bounds for da of size %d\n", index, da->size);
    exit(1);
  }
  return da->data[index];
}

int mdda_get(mdda_s *mdda, int c, int r) {

  mdda = mdda_get_node(mdda, c);
  return da_get((mdda -> da), r);
}

mdda_s * mdda_get_node(mdda_s *mdda, int c){

  int j;
  mdda_s * next_mdda = (mdda->root);
  int sz = next_mdda -> size;
  if (c >= sz || c < 0) {
    printf("mdda_get_node: Index %d out of bounds for mdda of size %d\n", c, mdda->size);
    exit(1);
  }

  for (j=0; j<=sz; j++) {
    mdda = next_mdda;
    if ((mdda -> idx) == c) {
      break;
    } else {
      next_mdda = (mdda->next);
    }
  }

  return next_mdda;
}


void da_set(da_s *da, int index, int value) {
  // zero fill the da up to the desired index
  while (index >= da->size) {
    da_append(da, 0);
  }

  // set the value at the desired index
  da->data[index] = value;
}


void mdda_set(mdda_s *mdda, int c, int r, int value) {

  mdda = mdda_get_node(mdda, c);
  da_set((mdda->da), r, value);

}

void
mdda_set_branch(mdda_s *t, /* tree mdda */
                mdda_s *b, /* branch mdda */
                int t_idx, /* index of the branch inside its mdda */
                int b_idx /* index of the tree inside its mdda */
                ){

  t = mdda_get_node(t, t_idx);
  b = mdda_get_node(b, b_idx);

  t -> branch = b;
}

void da_double_capacity_if_full(da_s *da) {
  if (da->size >= da->capacity) {
    // double da->capacity and resize the allocated memory accordingly
    da->capacity *= 2;
    da->data = realloc(da->data, sizeof(int) * da->capacity);
  }
}

void mdda2s(mdda_s * mdda){
  int j,k,l,m; /* looping variables */
  int jgrid,kgrid; /* looping variables */

  int n_fs,n_is;
  int n_gs = mdda_get(mdda, 0, 0);
  int gs_idx,is_idx,fs_idx;


  mdda_s * root_mdda = (mdda -> root);
  mdda_s * curr_mdda = mdda;
  mdda_s * next_mdda = (mdda -> next);

  /* get a pointer to the 0th column of iis */

  mdda_s * iis = root_mdda -> branch;

  curr_mdda = next_mdda;
  printf( "  -printing the content of the screened state matrix:");
  for (j=1; j<n_gs+1; j++) {
    next_mdda = (curr_mdda -> next);
    gs_idx = mdda_get(mdda, 0, j);
    printf( "\n   gs[%d/%d] = %d\n", j, n_gs, gs_idx);
    n_is = mdda_get(mdda, j, 0);
    for (k=1; k<n_is+1; k++) {
      is_idx = mdda_get(curr_mdda, j, k);
      printf( "     is[%d/%d] = %d\n", k, n_is, is_idx);


      if (mdda_intinint(iis, is_idx)) {

        /* the IS index is in iis, meaning we can read and print final state
         indices associated with that state*/
        for (l=1; l<n_is+1; l++) {
          if (mdda_get(iis, 0, l) == is_idx) {
            break;
          }
        }

        n_fs = mdda_get(iis, l, 0);
        /* l is the index of the column containing the final state indices
         we're interested in */
        for (m=1; m<n_fs+1; m++) {
          fs_idx = mdda_get(iis, l, m);
          printf( "       fs[%d/%d] = %d\n", m, n_fs, fs_idx);
        }
      }
      /* check if the IS idx is in the iis index array on root*/
      /* print error message if it isnt */
      /* else, repeat the printing process above */
    }
  }

  /* while((next_mdda = (curr_mdda -> next)) != NULL ){ */
  /*   if ((curr_mdda -> idxs) == ) { */
  /*     /\* we found the right intermediate state data *\/ */

  /*   } */

  /*   curr_mdda = next_mdda; */
  printf( "\n\n" );
}


int
mdda_intinint(mdda_s * mdda, int val){
  int j;
  int * dat = ((mdda -> da) -> data);
  int sz = ((mdda -> da) -> size);
  return intinint(dat, val, sz);
}

void da_free(da_s *da) {
  free(da->data);
}

void mdda_free(mdda_s *mdda) {

  mdda_s * last_mdda;
  mdda_s * next_mdda;
  mdda_s * curr_mdda = mdda;
  int * data;
  while((next_mdda = (curr_mdda -> next)) != NULL ){
    last_mdda = curr_mdda;
    curr_mdda = next_mdda;
    data = (last_mdda -> da)->data;

    free(data);
    free(last_mdda);
  }
}
