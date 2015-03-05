#ifndef DA_H
#define DA_H

#define DA_INITIAL_CAPACITY 100
#define MDDA_INITIAL_CAPACITY 100
typedef struct mdda_s mdda_s;
typedef struct da_s da_s;

// Define a da type
struct da_s{
  int size;      // slots used so far
  int capacity;  // total available slots
  int *data;     // array of integers we're storing
};

struct mdda_s{
  int size;      // slots used so far
  int idx;
  int mdda_capacity; /* capacity of new mddas in the llist */

  da_s *da;
  mdda_s * root;
  mdda_s * head; /* the last created node */
  mdda_s * next; /* pointer to the next mdda in the array */
}; /* multi-dimensional dynamic array */

void da_init(da_s *da);

mdda_s * mdda_init();

void da_append(da_s *da, int value);

void mdda_append(mdda_s *mdda, int val);

int da_get(da_s *da, int index);

int mdda_get(mdda_s *mdda, int c, int r);

void da_set(da_s *da, int index, int value);

void mdda_set(mdda_s *mdda, int c, int r, int value);

void da_double_capacity_if_full(da_s *da);

/* void mdda_double_capacity_if_full(mdda_s *da); */

void da_free(da_s *da);

void mdda_free(mdda_s *mdda);

#endif /* DA_H */
