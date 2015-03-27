#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <signal.h>
#include <string.h>
#include <stdarg.h> /* for var args */
#include "k_meansl.h"
#include "std_char_ops.h"
#include "std_num_ops.h"
#include "std_f.h"
#include "get_numsl.h"
#include "parse_input.h"
#include "input_formats.h"
#include "sci_const.h"
#include "dynarray.h"
#include "rmap_structs.h"
#include "info_ll.h"
#include "e_state_ll.h"

#define BUF_SIZE 256

int n_info_nodes; /* counter for number of  */

info_node info_node_root;

info_node
get_inode (char * id
           ){

  info_node curr_info_node = info_node_root;
  info_node next_info_node;

    /* locate the info_node corresponding to the file name input argument */
  while(strstr((curr_info_node -> str_id),id) == NULL) {

    next_info_node = curr_info_node -> next;
    curr_info_node = next_info_node;

    if (curr_info_node == NULL) {
      fprintf(stderr, "parse_input.c: get_inode, no info node can be found\
 with a str_id == %s\n", id);
      printf( "program terminating due to the previous error.\n");
      exit(1);
    }
  }

  return curr_info_node;
}

mdda_s *
screen_states (char * fn_infile,
             int n_args,
             ...){
  printf( "  -screening states\n\n" );
  info_node inode = get_inode(fn_infile);

  int j,k; /* looping variables */
  int n_states;
  int gs_idx; /* ground state index */
  int is_idx; /* intermediate state index */

  int n_gs = inode -> n_gs;
  int n_is = inode -> n_is;
  int n_fs = inode -> n_fs;

  /* number of screen GS, IS, and FS */
  int n_sgs = 0;
  int n_sis = 0;
  int n_sfs = 0;

  int * tmp_idxs;

  /* arrays used for storing the screened..  */
  mdda_s * igs; /* .. ground states */
  mdda_s * iis; /* .. intermediate states */

  /* some hard-coded threshold values that will get replaced by the call
   values */
  double bw_sum = inode -> bw_sum;

  double thrsh_vals[3] = {0.5, 0.2, 0.2};
  double tmp_thrsh;

  double * tmp_trans;

  n_states = inode -> n_states;

  e_state root_state = inode -> root_e_state;
  e_state curr_state = root_state;
  e_state next_state;
  e_state gs_bm; /* ground state bookmark */
  e_state is_bm; /* intermediate state bookmark */

  va_list argv;
  va_start(argv, n_args);

/* initialize the data structures used to store the screened states */
  igs = mdda_init();
  iis = mdda_init();

  for (j=0; j<n_args; j++) {
    tmp_thrsh = (double)va_arg(argv, double); /* grab the next vararg */
    thrsh_vals[j] = tmp_thrsh;
  }
  /* for a description of the screening algorithm below, look up this function
     in the parse_input.h file */
  /* fprintf(stderr, "\n\n=======Valgrind eject point=======\n\n"); */
  /* exit(1); */
  n_sgs = 0;
  while((next_state = curr_state -> next) != NULL){
    gs_idx = curr_state -> state_idx;
    if ((curr_state -> type) == 1) { /* if we found a ground state */
      printf( "GS: %d %le\n", gs_idx,((curr_state ->  bw)/bw_sum));
      /* perform stage 1 of the screening process */
      /* if (((curr_state ->  bw)/(inode -> bw_sum)) >= thrsh_vals[0]) { */
      if (((curr_state ->  bw)/bw_sum) >= thrsh_vals[0]) {
        /* we found a ground state with high enough boltzmann weight */
        gs_idx = curr_state -> state_idx;
        n_sgs++;
        mdda_set(igs, 0, 0, n_sgs);
        mdda_set(igs, 0, n_sgs, gs_idx );
        /* printf( "found one GS at %d!\n", gs_idx ); */
        printf( "\ngs[%d] = %d %le\n", n_sgs,gs_idx,((curr_state ->  bw)/bw_sum));
        /* store a pointer to the final state matrix in the ground state root node */

        /* sleep(1); */
        /* store a bookmark of where we found the last ground state */
        gs_bm = next_state;

        /* perform stage 2 of the screening process */
        n_sis = 0;
        tmp_trans = curr_state -> t_moms;
        tmp_idxs = curr_state -> idxs_to;
        for (j=0; j<((curr_state -> n_tfrom)-1); j++) {
          /* printf( "IS:j = %d, sidx = %d, trans = %le, reltrans = %le, thrsh = %le\n", j, tmp_idxs[j], tmp_trans[j], (tmp_trans[j]/(inode -> mt_is)), thrsh_vals[1]); */
          /* if IS transition has a transition moment above threshold */
          if ((tmp_trans[j]/(inode -> mt_is)) >= thrsh_vals[1]) {
            is_bm = curr_state;
            is_idx = tmp_idxs[j];
            /* printf( "found one IS at %d!\n", is_idx ); */

           /* locate the candidate intermediate state in the list and varify
               that the state is indeed an IS (has type == 2) */

            /* varify that it has been marked as an intermediate state */
            if (((curr_state = get_state(inode, is_idx)) != NULL)\
                && ((curr_state -> type) == 2)){
              /* from the ground state found above, we also found an intermediate
                 state with high enough relative transition moment */
              n_sis++;
              mdda_set(igs, n_sgs, 0, n_sis);
              mdda_set(igs, n_sgs, n_sis, is_idx);

              /* printf( "  is[%d] = %d %le\n", n_sis, mdda_get(igs, n_sgs,n_sis), (curr_state -> t_moms)[j]/(inode -> mt_is)); */

              n_sfs = 0;

              /* before this is done, check so that the specific is_idx isnt
                 in the iis already, in which case we would find its index in
                 the 0th column of that mdda */
              if (mdda_intinint((iis->root), is_idx) == 0) {

                mdda_set(iis, 0, 0, n_sis);
                mdda_set(iis, 0, n_sis, is_idx);

                /* perform stage 3 of the screening process */
                /* find final state transitions that pass stage 3 screening */
                for (k=0; k < (curr_state -> n_tfrom); k++) {
                  /* printf( "FS:j = %d, sidx = %d, trans = %le, reltrans = %le, thrsh = %le\n", k, ((curr_state -> idxs_to)[k]), ((curr_state -> t_moms)[k]), (((curr_state -> t_moms)[k])/(inode -> mt_fs)), thrsh_vals[2]); */
                  if ((((curr_state -> t_moms)[k])/(inode -> mt_fs)) >= thrsh_vals[2]) {
                    /* from the intermediate state found above, we also found a
                       final state transition with high enough relative
                       transition moment */
                    /* printf( "found one FS!\n" ); */
                    n_sfs++;
                    mdda_set(iis, n_sis, 0, n_sfs);
                    mdda_set(iis, n_sis, n_sfs, (curr_state -> idxs_to)[k]);

                    /* printf("    fs[%d][%d] = %d %le\n", n_sis, n_sfs, mdda_get(iis, n_sis, n_sfs), (((curr_state -> t_moms)[k])/(inode -> mt_fs))); */
                  }
                }
              }
              /* if we looped over all final state transitions for this specific
                 intermediate state and didnt find any FS transitions of high
                 enough intensity, remove that intermediate state from the list */
              /* if (n_sfs == 0) { */
              /*   n_sis--; */
              /*   mdda_set(igs, n_sgs, 0, n_sis); */
              /* } */
            }

            curr_state = is_bm; /* reset the state to the last initial state */
          }
        }
        /* jump back to the last ground state that was found in the list */
        next_state = gs_bm;
      } else {
      /* printf( "unacceptable bw: gs[%d] = %d %le\n", n_sgs,gs_idx,((curr_state ->  bw))); */
      }
    }
    curr_state = next_state;
  }
  /* fprintf(stderr, "\n\n=======Valgrind eject point=======\n\n"); */
  /* exit(1); */
  mdda_set_branch(igs, iis, 0, 0 );


  return (igs->root);
}

int
init_data_branch(double ** pi, /* parsed input */
                 int ns, /* n states */
                 int nt, /* n transitions */
                 char * fs /* input file name string */
                 ){
  int j; /* looping variables */

  /* initialize the linked list structure for the electronic states */
  init_state_ll(fs, ns, nt);

  /* use the extracted input data to set the state of the linked list */
  set_state_ll(pi, ns, nt, fs);

  /* free up parsed_input, since it now lives in the tree */
  for (j=0; j<5; j++) {
    free(pi[j]);
  }

  free(pi);

  return EXIT_SUCCESS;
}

info_node
init_info_node (char * s,
                int ns,
                int nt
                ){

  int j;

  static info_node last_info_node;
  int str_sz = strlen(s);
  char * si = malloc(strlen(s)+1); /* info node id string */
  info_node new_info_node;

  if((new_info_node = malloc(sizeof(struct info_node_s))) == NULL ){
    fprintf(stderr, "parse_input.c:function init_info_node, malloc: failed \
to allocate memory for \"new_info_node\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  for (j=0; j<str_sz; j++) {
    si[j] = s[j];
  }

  new_info_node -> str_id = si;
  new_info_node -> n_states = ns;
  new_info_node -> n_trans = nt;
  new_info_node -> n_gs = 0;
  new_info_node -> n_is = 0;
  new_info_node -> n_fs = 0;

  if (n_info_nodes == 0) { /* there is no root info node defined  */
    n_info_nodes = 1;
    new_info_node -> last = NULL;
    new_info_node -> next = NULL;
    info_node_root = new_info_node;
    last_info_node = info_node_root;
  }
  else { /* the root info node is already defined */
    n_info_nodes++;
    new_info_node -> last = last_info_node;
    last_info_node -> next = new_info_node;
    new_info_node -> next = NULL;
  }
  new_info_node -> idx = n_info_nodes-1;

  return new_info_node;
}

e_state
init_state_ll (char * str_id,
               int n_states,
               int n_trans
               ){

  int j,k,l; /* looping variables */
  int from,from_last,to;

  double tmp_energy;
  double tmp_bw;

  e_state root_state;
  e_state curr_state;
  e_state next_state;
  e_state last_state;

  info_node curr_info_node;

  curr_info_node = init_info_node(str_id, n_states, n_trans); /* this defines the root node and
                                              assigns it a non-NULL index */

  root_state = malloc(sizeof(struct e_state_s));
  root_state -> last = NULL;
  root_state -> list_idx = 0;
  root_state -> info = curr_info_node;

    /* by definition the parsed_input[:][0] should contain information on the
   lowest energy ground state. set the root state node type accordingly. */
  curr_info_node -> root_e_state = root_state;
  curr_state = root_state;

  /* construct the list structure */
  for (j=1; j<(n_states+1); j++) {
    /* note that the idxs_to and t_moms arrays dont get allocated here. this
     is done in the set_state_ll function */
    next_state = malloc(sizeof(struct e_state_s));
    next_state -> last = curr_state;
    next_state -> list_idx = j;
    next_state -> info = curr_info_node;
    curr_state -> next = next_state;
    curr_state -> type = 0;
    curr_state = next_state;
  }
  curr_state -> next = NULL;

  return root_state;
}

int
set_state_node (e_state st,
                int s_idx,
                int * idxs_buf,
                double * evals_buf,
                double * moms_buf,
                int n_trs_from,
                double e_rel,
                double e
                ){

  int j; /* looping variables */
  int * idxs;
  double * tm;
  double * ev;

  if((idxs = malloc(n_trs_from*sizeof(int))) == NULL ){
    fprintf(stderr, "parse_input:function set_state_node, malloc: failed \
to allocate memory for \"idxs_to\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  if((tm = malloc(n_trs_from*sizeof(double))) == NULL ){
    fprintf(stderr, "parse_input:function set_state_node, malloc: failed \
to allocate memory for \"tmoms\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  if((ev = malloc(n_trs_from*sizeof(double))) == NULL ){
    fprintf(stderr, "parse_input:function set_state_node, malloc: failed \
to allocate memory for \"tmoms\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  for (j=0; j<n_trs_from; j++) {
    idxs[j] = idxs_buf[j];
    tm[j] = moms_buf[j];
    ev[j] = evals_buf[j];
  }

  st -> state_idx = s_idx;
  st -> e_val = e;
  st -> n_tfrom = n_trs_from;
  st -> idxs_to = idxs;
  st -> t_moms = tm;
  st -> e_vals = ev;
  st -> max_tmom = get_maxl(tm, n_trs_from);

  return EXIT_SUCCESS;
}

int
set_symtrans (info_node inode
              ){
  int j,k; /* looping variables */

  int from,to, n_add, n_proc;
  int tmp_n_tfrom;
  int n_states = inode -> n_states;

  /* a register of all intermediate states whos symmetric transitions
   have been processed*/
  int * idxs_proc;

  double * tmp_e_vals;
  double * tmp_idxs_to;
  double * tmp_t_moms;

  double ** sym_dat;

  e_state bookmark;
  e_state next_state;
  e_state curr_state;

  if((idxs_proc = malloc((n_states)*sizeof(int))) == NULL ){
    fprintf(stderr, "parse_input:function set_symtrans, malloc: failed \
to allocate memory for \"ixds_proc\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  if((sym_dat = malloc(3*sizeof(double*))) == NULL ){
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
        /* printf( "====new from = %d\n", from); */
        /* first off, check if the state is even in the list.*/
        if (is_state_inlist(inode, from)) {

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

              if (intinint(curr_state -> idxs_to, from, curr_state -> n_tfrom)) {
                /* printf( "  from %d to %d\n", from, curr_state -> state_idx); */
                sym_dat[0][n_add] = curr_state -> state_idx;
                sym_dat[1][n_add] = get_trans(curr_state, from);
                sym_dat[2][n_add] = curr_state -> e_val;
                /* printf( "  sym_dat = %d %le %le\n", (int)sym_dat[0][n_add], sym_dat[1][n_add], sym_dat[2][n_add]); */
                sym_dat[0][0] = ++n_add-1;
                /* add the index of this inner state to the matrix */

              }
              curr_state = next_state;
            }

            curr_state = get_state(inode, from);
            idxs_proc[n_proc++] = from;

            /* the data we need for all symmetric transitions to the "from" state
               is now stored in sym_dat. update the from state with that data. */
            curr_state = get_state(inode, from);

            curr_state -> idxs_to = appc_d((curr_state -> idxs_to), &sym_dat[0][1] \
                                           , (curr_state -> n_tfrom), sym_dat[0][0]);

            curr_state -> t_moms = appc((curr_state -> t_moms), &sym_dat[1][1] \
                                        , (curr_state -> n_tfrom), sym_dat[0][0]);

            curr_state -> e_vals = appc((curr_state -> e_vals), &sym_dat[2][1]\
                                        , (curr_state -> n_tfrom), sym_dat[0][0]);

            curr_state -> n_tfrom += sym_dat[0][0];

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
  }

  printf( "proc_idx=" );
  for (j=0; j<n_proc; j++) {
    printf( "%d, ", idxs_proc[j]);
  }
  printf( "\n" );
  printf( "====done\n" );

  for (j=0; j<3; j++) {
    free(sym_dat[j]);
  }

  free(sym_dat);
  free(idxs_proc);
  e_statelist2s(inode, 1);

  fprintf(stderr, "\n\n=======Valgrind eject point=======\n\n");
  exit(1);
}

int
set_state_ll (double ** parsed_input,
              int n_states,
              int n_trans,
              char * id) {

  int j,k,l,m; /* looping variables */
  int from_last = parsed_input[0][0];
  int from,to,sorted_idx,tmp_idx;
  int tmp_ntfrom, tmp_type;
  int SYM = 1;
  int * g;
  int * ito;
  int ** groups;

  double tmax_is,tmax_fs; /*  */
  double tmp_energy;
  double bw_s; /* sum of boltzmann weight */

  double max_tm = 0;
  info_node curr_info_node;

  e_state end_state; /* pointer to the last state in the llist */
  e_state next_state;
  e_state curr_state;
  e_state tmp_state;

  int * tmp_idxs;
  double * tmp_tmoms;
  double * tmp_evals;
  double * e_vals;


  if((tmp_idxs = malloc(n_states*sizeof(int))) == NULL ){
    fprintf(stderr, "parse_input:function init_state_ll, malloc: failed \
to allocate memory for \"tmp_idxs\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  if((tmp_tmoms = malloc(n_states*sizeof(double))) == NULL ){
    fprintf(stderr, "parse_input:function init_state_ll, malloc: failed \
to allocate memory for \"tmp_tmoms\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  if((tmp_evals = malloc(n_states*sizeof(double))) == NULL ){
    fprintf(stderr, "parse_input:function init_state_ll, malloc: failed \
to allocate memory for \"tmp_evals\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  if((e_vals = malloc(n_states*sizeof(double))) == NULL ){
    fprintf(stderr, "parse_input:function init_state_ll, malloc: failed \
to allocate memory for \"e_vals\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  if((groups = malloc(2*sizeof(int*))) == NULL ){
    fprintf(stderr, "parse_input:function init_state_ll, malloc: failed \
to allocate memory for \"groups\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  for (j=0; j<2; j++) {
    if((groups[j] = malloc((n_states+1)*sizeof(int))) == NULL ){
      fprintf(stderr, "parse_input:function init_state_ll, malloc: failed \
to allocate memory for \"groups\"\n");
      printf( "program terminating due to the previous error.\n");
      exit(1);
    }
  }

  curr_info_node = get_inode(id);
  curr_state = curr_info_node -> root_e_state;

  /* the first node in the llist is the lowest energy ground state */
  curr_state -> type = 1;

  /* grab the transitions and energies */
  for (k=0,j=0,l=0; j<=n_trans; j++) { /* j = read head for the parsed_input matrix */

    if (j < n_trans) {
      from = parsed_input[0][j];
      to = parsed_input[1][j];
    } else {
      from = 0;
    }

    if (from != from_last) { /* we just started reading values for transitions
                                from a new state */
      /* initialize the next state in the ll */
      tmp_energy = parsed_input[2][j-1];

      /* we need these energies to sort the states later */
      e_vals[l] = tmp_energy;
      set_state_node(curr_state, from_last, tmp_idxs, tmp_evals, tmp_tmoms, k,\
                     parsed_input[2][0], tmp_energy);

      /* printf( "\nto=%d, e_val=%le, type=%d\n",(curr_state -> state_idx), \ */
      /*         (curr_state -> e_val), (curr_state -> type)); */
      /* for (m=0; m<k; m++) { */
      /*   printf( "from:%d, ev=%le, mom[%d]=%le\n", (curr_state -> idxs_to)[m],\ */
      /*           (curr_state -> e_vals)[m], m, (curr_state -> t_moms)[m]); */
      /* } */
      /* sleep(1); */
      if (((curr_state -> next) == NULL) && ((l+1) >= n_states)) {
/*       if ((l+1) >= n_states) { */
        fprintf(stderr, "parse_input.c:set_state_ll, the list of electronic \
states ended before all input data was transfered. there are %d out of %d \
states left to process.\n", l, n_states);
        printf( "program terminating due to the previous error.\n");
        exit(1);
      }
      if (l >= (n_states-1)) {
        break;
      }

      next_state = curr_state -> next;
      curr_state = next_state;

      from_last = from;

      k=0;
      j--;
      l++;

    }
    else if(j < n_trans) {
      tmp_idxs[k] = parsed_input[1][j];
      tmp_evals[k] = parsed_input[3][j];
      tmp_tmoms[k] = parsed_input[4][j];
      if (tmp_tmoms[k] > max_tm) {
        max_tm = tmp_tmoms[k];
      }
      k++;  /* increase counter for transitions counted */
    }
  }
  /* printf( "%d, %d, %d\n",n_states, l,n_trans); */
  /* fprintf(stderr, "\n\n=======Valgrind eject point=======\n\n"); */
  /* exit(1); */
  end_state = curr_state->last;

  /* printf( "%d\n", curr_state->state_idx); */
  /* printf( "%d\n", (curr_state-> last) -> state_idx); */
  /* printf( "%d\n", end_state -> last -> state_idx); */
  /* printf( "%d\n", end_state -> state_idx); */
  /* printf( "%d\n", curr_info_node -> root_e_state -> state_idx); */
  /* printf( "%d\n", l); */
  /*     fprintf(stderr, "\n\n=======Valgrind eject point=======\n\n"); */
  /* exit(1); */
  /* remove any unused nodes in the llist iterate forwards and free up\
     the unused nodes. */
  curr_info_node -> n_states = l;
  n_states = l;
  /* if (l < n_states) { */

  /*   for (j=0; j<n_states-l; j++) { */
  /*     printf( "freed!\n" ); */
  /*     next_state = curr_state -> next; */
  /*     free(curr_state); */
  /*     curr_state = NULL; /\* destroy the pointer *\/ */
  /*     curr_state = next_state; */
  /*   } */
  /* } */

  /* the last node in the llist is the highest energy intermediate state */
  end_state -> type = 2;
  end_state -> next = NULL;

  /* use the k-means algorithm to do a preliminary sorting of the states */
  k_meansl(e_vals, groups, l);

  /* based on the fact that the transitions between states inside the same
     group will be too low to get included in the data tree, we can sort
     out any state indices that ended up in the wrong category in the k_means
     sorting above.*/

  tmax_is = tmax_fs = -1;
  bw_s = 0;
  for (j=0; j<2; j++) {
    for (k=1; k<=groups[j][0]; k++) {

      sorted_idx = groups[j][k];
      /* check if the sorted_idx state is even in the llist of electronic states */

      if ((curr_state = get_state(curr_info_node, sorted_idx)) != NULL) {
        curr_state -> type = j+1;
        /* store the sum of boltzmann weights to later on use it for the state
           screening process */

        if (j == 0) {

          if ((curr_state -> max_tmom) > tmax_is) {
            tmax_is = (curr_state -> max_tmom);
          }

          curr_state -> bw = get_rbdist(parsed_input[2][0],curr_state -> e_val);
          if (curr_state -> state_idx != 1) {
            bw_s += curr_state -> bw;
          }

          (curr_info_node -> n_gs) += 1;

        } else {

          if ((curr_state -> max_tmom) > tmax_fs) {
            tmax_fs = (curr_state -> max_tmom);
          }

          (curr_info_node -> n_is) += 1;
        }
        /* printf( "\n\n"); */
        /* /\* count the state types *\/ */
        /* curr_state = curr_info_node -> root_e_state; */
        /* while((next_state = curr_state -> next) != NULL){ */
        /*   printf( "state %d, type = %d\n", curr_state -> state_idx, curr_state -> type); */
        /*   curr_state = next_state; */
        /* } */
      }
    }
  }

  curr_info_node -> bw_sum = bw_s;
  curr_info_node -> mt_is = tmax_is;
  curr_info_node -> mt_fs = tmax_fs;

  /* finally, if we want to take symmetric transitions into account, add the right ground
   states to the corresponding final states */
  if (SYM == 1) {
    set_symtrans(curr_info_node);
  }

  for (j=0; j<2; j++) {
    free(groups[j]);
  }
  free(groups);

  /* trim off any nodes at the end of the list that hasnt gotten allocated. */
  free(tmp_idxs);
  free(tmp_tmoms);
  free(tmp_evals);
  free(e_vals);

  /* e_statelist2s(curr_info_node); */

  return EXIT_SUCCESS;
}

int
parse_input_molcas (char * fn_infile) {

  int j,k,l,m,i,n,k_its; /* control loop variables */
  int mode; /* string matching mode flag */
  int match_start;
  int match_end;
  int last_int;
  int next_int;
  int idx_from;
  int idx_to;
  int n_states,n_trans;

  int * num_idxs1;
  int * num_idxs2;

  double ** parsed_input;
  double * e_eigval;
  double ** trans_idxs;
  double * t_mom;

  FILE * fp_tmpdata;
  FILE * fp_infile;

  const int n_lookup_str = 4; /* number of strings used for searching the input file */
  const int n_maxmatch = n_lookup_str/2;
  int n_matchfound = 0;
  int match_vals[2] = {0,0}; /* place to store the indexes of the lines containing the
                              matches */
  /* place to store the indexes of the lines containing the
     matches. each data block will start after,  and end after the offset values.*/
  int match_ln[n_lookup_str];

  int c; /* temporary char for storing input file characters */
  /* we are looking for four strings in the output file: one at the beginning
     of each data block, and one at the end. */
  char * str_buf = malloc(BUF_SIZE*2);

  const char s1[26] = "        Relative EVac(au)";
  const char s2[51] = " Weights of the five most important spin-orbit-free";
  const char s3[18] = "         To  From";
  const char s4[52] = " ##################################################";
  const char * lookup_str[4] = {s1,s2,s3,s4};

  /* open the input file */
  if((fp_infile = fopen(fn_infile, "r")) == NULL) {
    fprintf(stderr,"parse_input: unable to open the input file %s.\n",fn_infile);
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  if((fp_tmpdata = fopen("/home/kimchi/dev/rmap/tmp/molcas_data.tmp", "w+")) == NULL) {
    fprintf(stderr,"parse_input: unable to open the output file %s.\n",fn_infile);
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  k = 0; /* index for tmp_string */
  l = 0; /* index for lookup string */
  m = 0; /* index for string matches */
  mode = 0; /* start of in string search mode */

  /* read char by char */
  for (j=0; ((c = fgetc(fp_infile)) != EOF); j++, k++) {
    str_buf[k] = (char) c;

    /* check to see if there are any matches left to catch in the input data */
    if (n_matchfound < n_maxmatch) {

      /* keep extracting characters from the input data until an entire line
         has been stored in the temporary str_buf buffer */
      if (str_buf[k] == '\n') {

        /* check every line for a matching substring */
        /* mode = 1 and a line match means that we reached
           the end of this data block */
        if (strstr(str_buf,lookup_str[l]) || (mode == 1)) {

          /* we found the first substring, the data we're looking for is
             inside the coming block of text. switch to mode 1.*/
          if (mode == 0) {
            match_ln[m++] = j;
            l++;
            mode = bin_flip(mode);
          }

          else { /* we're in count the matched lines mode */
            /* if (regexec(&c_regexs[m-1],str_buf,0,NULL,0)) { */
            if (isdigitin(str_buf,k-1) == 1) {

              if (isdashes(str_buf,k) == 0){
                /* spit the line to file */
                for (n=0; n<BUF_SIZE; n++) {
                  if (n<k) {
                    fputc(str_buf[n],fp_tmpdata);
                  } else {
                    fputc(' ',fp_tmpdata);
                  }
                }
                fputc('\n',fp_tmpdata);
              }
              match_vals[n_matchfound]++;
            }

              /* skip empty and dashed lines */
            else if((isempty(str_buf,k) == 0) && (isdashes(str_buf,k) == 0)) {
              /*                        string of the block *\/ */
              l++;
              mode = bin_flip(mode); /* switch back to mode 0 */
              n_matchfound++; /* one data block was found */
            }
          }
        }
        k = 0;
      }
    } else {
      break;
    }
  }

  /* we now know the number of states used in the molcas calculation, as well
     as the number of possible transitions.*/
  n_states = match_vals[0];
  n_trans = match_vals[1];

  /* convert the number of states and transitions read to offsets in bytes
   that can be used with the fseek function */
  match_ln[0] = 0;
  match_ln[1] = (n_states+1)*256;
  match_ln[2] = (n_states)*256;
  match_ln[3] = match_ln[2] + ((n_trans+1)*256);
  /* for (j=0; j<4; j++) { */
  /*   printf( "match = %d\n", match_ln[j]); */
  /* } */

  if((parsed_input = malloc(5*sizeof(double*))) == NULL ){
    fprintf(stderr, "parse_input_molcas, malloc: failed to allocate memory for\
 \"input_data\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  /* allocate space for the "parsed input matrix" that will be filled with data
     in the remaining sections of this function */
  for (j=0; j<5; j++) {
      if((parsed_input[j] = malloc(n_trans*sizeof(double))) == NULL ){
        fprintf(stderr, "parse_input_molcas, malloc: failed to allocate memory \
for pointers in \"input_data\"\n");
        printf( "program terminating due to the previous error.\n");
        exit(1);
    }
  }

  /* storage for the energy eigenvalues */
  if((e_eigval = malloc(n_states*sizeof(double))) == NULL ){
    fprintf(stderr, "parse_input_molcas, malloc: failed to allocate memory for\
 \"e_eigval\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  /* storage for the transition moments */
  if((t_mom = malloc(n_trans*sizeof(double))) == NULL ){
    fprintf(stderr, "parse_input_molcas, malloc: failed to allocate memory for\
 \"e_eigval\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }
    /* storage for the transition indexes, column 1 is from a state
     column 2 is to state index */
  if((trans_idxs = malloc(2*sizeof(double*))) == NULL ){
    fprintf(stderr, "parse_input_molcas, malloc: failed to allocate memory for\
 \"e_eigval\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  for (j=0; j<2; j++) {
    if((trans_idxs[j] = malloc(n_trans*sizeof(double))) == NULL ){
        fprintf(stderr, "parse_input_molcas, malloc: failed to allocate memory for\
 \"e_eigval\"\n");
        printf( "program terminating due to the previous error.\n");
        exit(1);
      }
  }

  if((num_idxs1 = malloc(sizeof(int))) == NULL ){
    fprintf(stderr, "parse_input_molcas, malloc: failed to allocate memory for\
 \"e_eigval\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  if((num_idxs2 = malloc(3*sizeof(int))) == NULL ){
    fprintf(stderr, "parse_input_molcas, malloc: failed to allocate memory for\
 \"e_eigval\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  /* int num_idxs1[1] = {2}; */
  num_idxs1[0] = 1;
  int n_idxs1 = 1;

  num_idxs2[0] = 0;
  num_idxs2[1] = 1;
  num_idxs2[2] = 6;
  int n_idxs2 = 3;

  for (j=0; j<n_lookup_str; j+=2) {

    l = 0; /* index for string matches */
    m = 0;

    match_start = match_ln[j];
    match_end = match_ln[j+1];
    /* printf( "\nstart = %d , end = %d\n", match_start, match_end); */
    fseek(fp_tmpdata, match_start, 0);

    k_its = match_end-match_start;

    for (k=0; k<k_its; ) {

      if ((c = fgetc(fp_tmpdata)) == EOF) {
        break;
      }

      str_buf[l] = (char)c;

      if (j == 0) {
        k++;
      }

      if ((str_buf[l] == '\n') && (l > 0)) { /* dont send blank lines */

        if ((j == 0) && (isempty(str_buf,l) != 1)) {

          /* extract energy eigenvalues and state indexes */
          get_numsl(str_buf,num_idxs1,l,n_idxs1,&e_eigval[m]);
          /* printf( "e_eigval[%d] = %le\n", m, e_eigval[m]); */

          m++;
        }

        if ((j == 2) && (isempty(str_buf,l) != 1)) {
           /* extract transition moments and transition indexes */
          get_numsl(str_buf,num_idxs2,l,n_idxs2,&trans_idxs[0][m],\
                    &trans_idxs[1][m],&t_mom[m]);
          /* printf( "to %le from %le, %le \n", trans_idxs[0][m], trans_idxs[1][m], t_mom[m]); */
          /* sleep(1); */
          m++;
        }
        l=0; /* reset the buffer write head to start reading a the next line */
      }
      l++;
    }
  }

  /* finally, store the data in the parsed_input matrix for the parse_input function  */
  for (j=0; j<n_trans; j++) {
    idx_from = trans_idxs[0][j];
    idx_to = trans_idxs[1][j];
    parsed_input[0][j] = idx_from;
    parsed_input[1][j] = idx_to;
    parsed_input[2][j] = e_eigval[idx_from-1];
    parsed_input[3][j] = e_eigval[idx_to-1];
    parsed_input[4][j] = t_mom[j];
  }

  /* for (k=0; k<n_trans; k++) { */
  /*   printf( "%le   %le   %le   %le   %le\n", parsed_input[0][k], parsed_input[1][k], parsed_input[2][k], parsed_input[3][k], parsed_input[4][k]); */
  /* } */

  fclose(fp_infile);
  fclose(fp_tmpdata);

  free(str_buf);
  free(num_idxs1);
  free(num_idxs2);

  free(e_eigval);

  free(t_mom);
  free(trans_idxs);


  return init_data_branch(parsed_input,n_states,n_trans,fn_infile);
}

int
parse_input (char * fn_infile, /* name of input file */
             int len_fn){

  int j, k, l; /* looping variables */
  int n_states,n_trans;

  char format[BUF_SIZE] = {0};
  e_state root_e_state;

  /* after having been used in a reference call to a parsing function the
     parsed_input array should contain an n_trans * 5 matrix loaded with
     data for each electronic transition that was calculated:
     [idx_f,idx_t,e_f,e_t,osc] where...
     idx_f = the index the electronic state FROM which the transition took place
     idx_t = the index the electronic state where the transition went TO
     e_f = energy corresponding to idx_f
     e_t = energy corresponding to idx_t
     mom = the transition moment for the transition
  */

  /* loop over the input file name and extract the ending */
  for (j=len_fn; j>0; j--) {
    format[j] = fn_infile[j-3];
    if (format[j] = '.') {

      if (strcmp(format,MOLCAS_FORMAT)){
        parse_input_molcas(fn_infile);
        break; /* for */
      }

      else {
        fprintf(stderr, "format not found in input_formats.h. are you sure the\
 file extension %s is correct?\n",format);
        printf( "program terminating due to the previous error.\n");
        exit(1);
      }
    }
  }

  return EXIT_SUCCESS;
}
