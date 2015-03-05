#include <stdio.h>
#include <stdlib.h>
#include "dynarray.h"

void da_init(da_s *da) {
  // initialize size and capacity
  da->size = 0;
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
  (root_mdda->da)->size = 0;
  (root_mdda->da)->capacity = DA_INITIAL_CAPACITY;

  // allocate memory for da->data
  (root_mdda->da)->data = malloc(sizeof(int) * ((root_mdda->da)->capacity));

  root_mdda -> next = malloc(sizeof(mdda_s));
  root_mdda -> root = root_mdda;
  root_mdda -> head = root_mdda;
  curr_mdda = root_mdda;

  while(++j < MDDA_INITIAL_CAPACITY){
    printf( "initializing %d %d %d \n", j,  (curr_mdda->root) -> idx, root_mdda -> idx);
    next_mdda = (curr_mdda->next);
    curr_mdda = next_mdda;
    curr_mdda -> idx = j;
    curr_mdda -> size = 0; /* used to check if the da needs initialization */
    curr_mdda -> mdda_capacity = MDDA_INITIAL_CAPACITY;
    printf( "da_initing\n" );
    curr_mdda->da = malloc(sizeof(da_s));
    (curr_mdda->da)->size = 0;
    (curr_mdda->da)->capacity = DA_INITIAL_CAPACITY;

    // allocate memory for da->data
    (curr_mdda->da)->data = malloc(sizeof(int) * ((curr_mdda->da)->capacity));

    /* da_init(&(curr_mdda->da)); */
    printf( "da_inited\n" );
    printf( "allocating\n" );
    curr_mdda -> next = malloc(sizeof(mdda_s));
    printf( "allocated\n" );
    curr_mdda -> root = root_mdda;
    curr_mdda -> head = root_mdda;

    printf( "initialized %d \n\n", curr_mdda -> idx);
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
  if (index >= da->size || index < 0) {
    printf("da_get: Index %d out of bounds for da of size %d\n", index, da->size);
    exit(1);
  }
  return da->data[index];
}

int mdda_get(mdda_s *mdda, int c, int r) {

  int j;

  mdda_s * next = (mdda->root);

  int sz = next -> size;

  if (c >= sz || c < 0) {
    printf("mdda_get: Index %d out of bounds for da of size %d\n", c, mdda->size);
    exit(1);
  }

  /* loop from root and find the right mdda */
  for (j=0; j<=sz; j++) {
    mdda = next;
    if ((mdda -> idx) == c) {
      break;
    } else {
      next = (mdda->next);
    }
  }

  return da_get((mdda -> da), r);
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

  int j;

  mdda_s * next_mdda = (mdda->root);

  int sz = next_mdda -> size;

  /* if (c >= sz || c < 0) { */
  /*   printf("mdda_set: Index %d out of bounds for da of size %d\n", c, mdda->size); */
  /*   exit(1); */
  /* } */

  /* loop from root and find the right mdda */
  for (j=0; j<=sz; j++) {
    mdda = next_mdda;
    if ((mdda -> idx) == c) {
      break;
    } else {
      next_mdda = (mdda->next);
    }
  }

  da_set((mdda->da), r, value);

}

void da_double_capacity_if_full(da_s *da) {
  if (da->size >= da->capacity) {
    // double da->capacity and resize the allocated memory accordingly
    da->capacity *= 2;
    da->data = realloc(da->data, sizeof(int) * da->capacity);
  }
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
