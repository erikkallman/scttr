#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "appc.h"
#include "info_node.h"
#include "structs.h"
#include "estate.h"

info_node
get_inode (char * id
           ){

  info_node curr_inode = root_inode;
  info_node next_inode;

    /* locate the info_node corresponding to the file name input argument */
  while(strstr((curr_inode -> str_id),id) == NULL) {

    next_inode = curr_inode -> next;
    curr_inode = next_inode;

    if (curr_inode == NULL) {
      fprintf(stderr, "parse_input.c: get_inode, no info node can be found\
 with a str_id == %s\n", id);
      printf( "program terminating due to the previous error.\n");
      exit(1);
    }
  }

  return curr_inode;
}

info_node
init_inode (char * s,
            int ns,
            int nt
            ){

  int j;

  static info_node last_inode;
  int str_sz = strlen(s);
  int * el_idxs;
  char * si = malloc(strlen(s)+1); /* info node id string */
  info_node new_inode;

  estate * estate_list;

  if((new_inode = malloc(sizeof(struct info_node_s))) == NULL ){
    fprintf(stderr, "info_node.c:function init_inode, malloc: failed \
to allocate memory for \"new_inode\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  if((estate_list = malloc((ns+1)*sizeof(struct e_state_s))) == NULL ){
    fprintf(stderr, "info_node.c:function init_inode, malloc: failed \
to allocate memory for \"estate_list\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  if((el_idxs = malloc((ns+1)*sizeof(int))) == NULL ){
    fprintf(stderr, "info_node.c:function init_inode, malloc: failed \
to allocate memory for \"el_idxs\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  for (j=0; j<str_sz; j++) {
    si[j] = s[j];
  }

  new_inode -> str_id = si;
  new_inode -> n_states = ns;
  new_inode -> n_trans = nt;
  new_inode -> n_gfs = 0;
  new_inode -> n_is = 0;
  new_inode -> n_spec = 0;
  new_inode -> e_states = estate_list;
  new_inode -> ev_idxs = el_idxs;

  if (n_inodes == 0) { /* there is no root info node defined  */
    n_inodes = 1;
    new_inode -> last = NULL;
    root_inode = new_inode;
    last_inode = root_inode;
  }
  else { /* the root info node is already defined */
    n_inodes++;
    new_inode -> last = last_inode;
    last_inode -> next = new_inode;
  }
  new_inode -> next = NULL;
  new_inode -> idx = n_inodes-1;

  return new_inode;
}

int
set_symtrans (info_node inode
              ){
  int j,k; /* looping variables */

  int from,to, n_add, n_proc;
  int tmp_n_tfrom, from_type;
  int n_states = inode -> n_states;

  /* a register of all intermediate states whos symmetric transitions
   have been processed*/
  int * idxs_proc;

  double * tmp_e_vals;
  double * tmp_idxs_to;
  double * tmp_t_moms;

  double ** sym_dat;

  estate bookmark;
  estate next_state;
  estate curr_state;

  if((idxs_proc = malloc((n_states)*sizeof(int))) == NULL ){
    fprintf(stderr, "parse_input:function set_symtrans, malloc: failed \
to allocate memory for \"ixds_proc\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  if((sym_dat = malloc(4*sizeof(double*))) == NULL ){
    fprintf(stderr, "parse_input:function set_symtrans, malloc: failed \
to allocate memory for \"tmp_dat\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  for (j=0; j<3; j++) {
    if((sym_dat[j] = malloc((n_states+1)*sizeof(double))) == NULL ){
      fprintf(stderr, "parse_input:function set_symtrans, malloc: failed \
to allocate memory for \"sym_dat[%d]\"\n",j);
      printf( "program terminating due to the previous error.\n");
      exit(1);
    }
  }

  n_proc = 0;
  curr_state;
  next_state = inode -> root_e_state;
  /* e_statelist2s(inode,1); */
  while(next_state != NULL){
    curr_state = next_state;
    /* only check for symmetric transitions for non IS  */
    /* printf( "====start at state = %d\n", (curr_state -> state_idx)); */
    if ((curr_state -> type != 2) &&                            \
      intinint(idxs_proc,curr_state -> state_idx,n_proc) == 0) {
      /* printf( "====found state = %d\n", (curr_state -> state_idx)); */
      bookmark = curr_state;

      for (j=0; j<(curr_state -> n_tfrom); j++) {

        from = (curr_state -> idxs_to)[j];
        from_type = (inode -> e_states)[(inode -> ev_idxs)[from]] -> type;
        /* printf( "====new from = %d\n", from); */
        /* first off, check if the state is even in the list.*/
        if ((from_type != (bookmark -> type)) && is_state_inlist(inode, from)) {

          /* printf( "from %d to %d\n", from, curr_state -> state_idx); */
          /* check if the state has already been processed */
          if (intinint(idxs_proc, from, n_proc) == 0) {
            n_add = 1;
            sym_dat[0][0] = 0;
            sym_dat[1][0] = from;


            /* locate all non-intermediate states that have transitions to the
               "from" state defined above */
            curr_state = inode -> root_e_state;

            while((next_state = curr_state -> next) != NULL){

              /* skip transitions between states of the same type
                 (fs/gs -> fs/gs, etc.) */
              if (intinint(curr_state -> idxs_to, from, curr_state -> n_tfrom)) {
              /* if (intinint(curr_state -> idxs_to, from, curr_state -> n_tfrom)) { */
                /* printf( "  from %d to %d\n", from, curr_state -> state_idx); */
                sym_dat[0][n_add] = (double)curr_state -> state_idx;
                sym_dat[1][n_add] = get_trans(curr_state, from);
                sym_dat[2][n_add] = curr_state -> e_val;
                sym_dat[0][0] = (++n_add)-1;
                /* printf( "sym_dat[0][0] = %le\n", sym_dat[0][0]); */
                /* printf( "  sym_dat = %d %le %le\n", (int)sym_dat[0][n_add], sym_dat[1][n_add], sym_dat[2][n_add]); */
                /* add the index of this inner state to the matrix */
              }
              curr_state = next_state;
            }
            /* printf( "  assigning new data1\n" ); */
            curr_state = get_state_sil(inode, from);
            idxs_proc[n_proc++] = from;
            /* printf( "  assigning new data\n" ); */
            /* the data we need for all symmetric transitions to the "from" state
               is now stored in sym_dat. update the from state with that data. */
            curr_state = get_state_si(inode, from);

            curr_state -> idxs_to = appc_d((curr_state -> idxs_to), &sym_dat[0][1] \
                                           , (curr_state -> n_tfrom), (int)sym_dat[0][0]);

            curr_state -> t_moms = appc((curr_state -> t_moms), &sym_dat[1][1] \
                                        , (curr_state -> n_tfrom), (int)sym_dat[0][0]);

            curr_state -> e_vals = appc((curr_state -> e_vals), &sym_dat[2][1]\
                                        , (curr_state -> n_tfrom), (int)sym_dat[0][0]);

            curr_state -> n_tfrom += sym_dat[0][0];

            free(curr_state -> ttypes);

            curr_state -> ttypes = malloc((curr_state -> n_tfrom)*sizeof(int));

            /* printf( "sym_dat[0][0] = %d\n", sym_dat[0][0]); */
            /* curr_state -> max_tmom = get_maxl(curr_state ->  t_moms, curr_state -> n_tfrom); */
            curr_state = bookmark;

            /* printf( "====back at state = %d\n", (curr_state -> state_idx)); */
          }
        }
        /* jump back to the bookmark to keep iterating from the last
         position in the linked list */
      }
      curr_state = bookmark;
    }
    next_state = curr_state -> next;
    /* if (curr_state -> state_idx == 6) { */
    /*   0cv */
    /* } */
  }

  /* printf( "proc_idx=" ); */
  /* for (j=0; j<n_proc; j++) { */
  /*   printf( "%d, ", idxs_proc[j]); */
  /* } */
  /* printf( "\n" ); */
  /* printf( "====done\n" ); */

  for (j=0; j<3; j++) {
    free(sym_dat[j]);
  }

  free(sym_dat);
  free(idxs_proc);
}
