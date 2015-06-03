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
            int * mom,
            int ns,
            int nt
            ){

  int j,k;

  static info_node last_inode;
  int str_sz = strlen(s);

  int * el_idxs;
  char * si = malloc(strlen(s)+1); /* info node id string */
  double ** max_t;

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

  if((el_idxs = malloc((ns*2)*sizeof(int))) == NULL ){
    fprintf(stderr, "info_node.c:function init_inode, malloc: failed \
to allocate memory for \"el_idxs\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  /* mat_t[type][ttype (moment type)] */
  if((max_t = malloc(2*sizeof(double *))) == NULL ){
    fprintf(stderr, "info_node.c:function init_inode, malloc: failed \
to allocate memory for \"max_t\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  for (j=0; j<2; j++) {
    if((max_t[j] = malloc(mom[0]*sizeof(double))) == NULL ){
      fprintf(stderr, "info_node.c:function init_inode, malloc: failed \
to allocate memory for \"max_t[%d]\"\n",j);
      printf( "program terminating due to the previous error.\n");
      exit(1);
    }
    for (k=0; k<mom[0]; k++) {
      max_t[j][k] = -1;
    }
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
  new_inode -> mom_types = mom;
  new_inode -> mt = max_t;

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
  int tmp_n_tfrom, from_type, from_ttype;
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

    /* only check for symmetric transitions from ground states to states that
       have not already been processed  */
    /* printf( "====start at state = %d\n", (curr_state -> state_idx)); */
    if (curr_state -> type == 1){
      /* printf( "====found state = %d\n", (curr_state -> state_idx)); */
      /* store a reference to the current ground state to continue iterating
         from it later */
      bookmark = curr_state;

      /* loop over each transition from the current ground state and check .. */
      for (j=0; j<(curr_state -> n_tfrom); j++) {

        from = (curr_state -> idxs_to)[j];
        from_type = (inode -> e_states)[(inode -> ev_idxs)[from]] -> type;
        /* printf( "====new from = %d\n", from); */
        /* .. if the state is in the estate list, that is is a ground satate,
           and that it hasnt already been processed.*/
        if ((intinint(idxs_proc, from, n_proc) == 0) &&
            (from_type != (bookmark -> type)) && is_state_inlist(inode, from)
            /* && (from_type == 1) */
            ) {

          /* printf( "from %d to %d\n", from, curr_state -> state_idx); */
          n_add = 1;
          sym_dat[0][0] = 0;
          sym_dat[1][0] = from;

          /* locate all non-intermediate states that have transitions to the
             "from" state defined above */
          curr_state = inode -> root_e_state;

          while(curr_state != NULL){

            /* if there is a transition from another ground state G2 to the
               intermediate state "from" add the G2 state to the from
               transition list as well as a symmetric transition to and from
               that intermediate state */
            if (intinint(curr_state -> idxs_to, from, curr_state -> n_tfrom)
                && (curr_state-> type == 1)) {
              /* if (intinint(curr_state -> idxs_to, from, curr_state -> n_tfrom)) { */
              /* printf( "\n\nfrom %d,type %d to %d, type %d, bookmark = %d\n", from, from_type, curr_state -> state_idx, curr_state -> type, bookmark -> type); */
              sym_dat[0][n_add] = (double)curr_state -> state_idx;
              sym_dat[1][n_add] = get_trans(curr_state, from);
              sym_dat[2][n_add] = curr_state -> e_val;
              /* printf( "  sym_dat[0][0] = %le\n", sym_dat[0][0]); */
              /* printf( "  sym_dat = %d %le %le\n", (int)sym_dat[0][n_add], sym_dat[1][n_add], sym_dat[2][n_add]); */
              sym_dat[0][0] = (++n_add)-1;

              /* add the index of this inner state to the matrix */
            }
            next_state = curr_state -> next;
            curr_state = next_state;
          }
          /* did we find any symmetric transitions to the "from" state?  */
          if ((int)sym_dat[0][0] >= 0) {
            /* printf( "  assigning new data1\n" ); */
            /* sleep(1); */

            /* jump back to the from state */
            curr_state = get_state_sil(inode, from);
            idxs_proc[n_proc++] = from;
            /* printf( "  assigning new data\n" ); */

            /* the data we need for all symmetric transitions to the "from" state
               is now stored in sym_dat. update the from state with that data. */
            curr_state = get_state_sil(inode, from);

            curr_state -> idxs_to = appc_d((curr_state -> idxs_to), &sym_dat[0][1] \
                                           , (curr_state -> n_tfrom), (int)sym_dat[0][0]);

            curr_state -> t_moms = appc((curr_state -> t_moms), &sym_dat[1][1] \
                                        , (curr_state -> n_tfrom), (int)sym_dat[0][0]);

            curr_state -> e_vals = appc((curr_state -> e_vals), &sym_dat[2][1]\
                                        , (curr_state -> n_tfrom), (int)sym_dat[0][0]);

            curr_state -> n_tfrom += sym_dat[0][0];

            if (curr_state -> ttypes != NULL) {
              free(curr_state -> ttypes);
            }

            curr_state -> ttypes = malloc((curr_state -> n_tfrom)*sizeof(int));

            /* printf( "sym_dat[0][0] = %d\n", sym_dat[0][0]); */
            /* curr_state -> max_tmom = get_maxl(curr_state ->  t_moms, curr_state -> n_tfrom); */
            /* finally, jump back to the bookmark and continue iterating */
            curr_state = bookmark;

            /* printf( "====back at state = %d\n\n", (curr_state -> state_idx)); */
          }
        }
      }
      /* jump back to the bookmark to keep iterating from the last
         position in the linked list */
      curr_state = bookmark;
    }
    next_state = curr_state -> next;
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
