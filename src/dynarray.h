#ifndef DA_H
#define DA_H

#define DA_INITIAL_CAPACITY 500
#define MDDA_INITIAL_CAPACITY 500
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

  /* a pointer to an array in a different mdda */
  mdda_s * branch;

  mdda_s * root;
  mdda_s * head; /* the last created node */
  mdda_s * next; /* pointer to the next mdda in the array */
}; /* multi-dimensional dynamic array */

/* function isiniis

   * synopsis:

   * algorithm:

   * input:

   * output:

   * side-effects:

   */
int
isiniis (mdda_s * igs,
         mdda_s * iis
         );

void da_init(da_s *da);

/* any node will always have a size of minimum 1 since the 0th element is
   used to store the number of states in the array */
mdda_s * mdda_init();

void da_append(da_s *da, int value);

void mdda_append(mdda_s *mdda, int val);

int da_get(da_s *da, int index);

int mdda_get(mdda_s *mdda, int c, int r);

mdda_s * mdda_get_node(mdda_s *mdda, int c);

void da_set(da_s *da, int index, int value);


/* function mdda_set

   * synopsis:
   this function sets the element in column @c and row @r to the value @value.

   * algorithm:
   obtain a local reference to the @mdda node, use da_set to add the value to
   the element in the data array of the node, identified by the row number @r.
   use the da_append function to both zero-pad the array (if needed) and set
   the right element to the right number.

   * input:
   mdda : the multi-dimensional dynamic array column
   c : the column index of the mdda
   r : the row index of the column c
   value : the value the [c,r] element will be set to

   * output:    none

   * side-effects:
   the size of the data array of column node c will be increased as new values
   are added to it. also, the data array will be padded with zeros up until the
   [c,r] element if this is not yet allocated.

   */
void mdda_set(mdda_s *mdda, int c, int r, int value);

void
mdda_set_branch(mdda_s *t, /* tree mdda */
                mdda_s *b, /* branch mdda */
                int tree_idx, /* index of the branch inside its mdda */
                int branch_idx /* index of the tree inside its mdda */
                );

void da_double_capacity_if_full(da_s *da);

/* void mdda_double_capacity_if_full(mdda_s *da); */

/* function to debug print the contents of a mdda */
void mdda2s(mdda_s * mdda);

int
mdda_intinint(mdda_s * mdda, int val);

void da_free(da_s *da);

/* function mdda_show

   * synopsis:

   * algorithm:

   * input:

   * output:

   * side-effects:

   */
void
mdda_show (mdda_s * mdda);

void mdda_free(mdda_s *mdda);

#endif /* DA_H */
