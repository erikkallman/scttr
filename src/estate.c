#include <stdlib.h>
#include <stdio.h>
#include "k_meansl.h"
#include "estate.h"
#include "structs.h"
#include "info_node.h"
#include "appc.h"
#include "sci_const.h"
#include "std_num_ops.h"

void
sort_states (double * state_er,
             double * e_vals,
             int ** groups,
             int n_states){
  int j,k,l; /* looping variables */

  k = l = 1;
  for (j=0; j<n_states; j++) {
    if ((state_er[1] < e_vals[j]) && (state_er[2] > e_vals[j])) {
      groups[0][k] = j;
      groups[0][0] = k++;
    }
    else if ((state_er[3] < e_vals[j]) && (state_er[4] > e_vals[j])){
      groups[1][l] = j;
      groups[1][0] = l++;
    }
  }
}

int
set_estate_list (double * state_er,
                 double ** parsed_input,
                 int n_states,
                 int n_trans,
                 char * id) {

  int j,k,l,m; /* looping variables */
  int from_last = parsed_input[0][0];
  int from,to,sorted_idx,tmp_idx,state_bookmark;
  int tmp_ntfrom, tmp_type;
  int SYM = 1;
  int * g;
  int * ito;
  int ** groups;

  double tmp_energy;
  double bw_s; /* sum of boltzmann weight */

  double max_tm = 0;
  info_node curr_inode;

  estate end_state; /* pointer to the last state in the llist */
  estate next_state;
  estate curr_state;
  estate tmp_state;

  int * tmp_idxs;
  double * tmp_tmoms;
  double * tmp_evals;
  double * e_vals;


  if((tmp_idxs = malloc(n_states*sizeof(int))) == NULL ){
    fprintf(stderr, "parse_input:function init_estate_list, malloc: failed \
to allocate memory for \"tmp_idxs\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  if((tmp_tmoms = malloc(n_states*sizeof(double))) == NULL ){
    fprintf(stderr, "parse_input:function init_estate_list, malloc: failed \
to allocate memory for \"tmp_tmoms\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  if((tmp_evals = malloc(n_states*sizeof(double))) == NULL ){
    fprintf(stderr, "parse_input:function init_estate_list, malloc: failed \
to allocate memory for \"tmp_evals\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  if((e_vals = malloc(n_states*sizeof(double))) == NULL ){
    fprintf(stderr, "parse_input:function init_estate_list, malloc: failed \
to allocate memory for \"e_vals\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  if((groups = malloc(2*sizeof(int*))) == NULL ){
    fprintf(stderr, "parse_input:function init_estate_list, malloc: failed \
to allocate memory for \"groups\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  for (j=0; j<2; j++) {
    if((groups[j] = malloc((n_states+1)*sizeof(int))) == NULL ){
      fprintf(stderr, "parse_input:function init_estate_list, malloc: failed \
to allocate memory for \"groups\"\n");
      printf( "program terminating due to the previous error.\n");
      exit(1);
    }
  }

  curr_inode = get_inode(id);
  curr_state = curr_inode -> root_e_state;

  /* the first node in the llist is the lowest energy ground state */
  curr_state -> type = 1;

  /* grab the transitions and energies */
  for (k=0,j=0,l=0; j<=n_trans; j++) { /* j = read head for the parsed_input matrix */

    if (j < n_trans) {
      from = parsed_input[0][j];
      to = parsed_input[1][j];
    } else {
      /* store the data of the final state in the transition list */
      from = 0;
    }

    if (from != from_last) { /* we just started reading values for transitions
                                from a new state */
      /* initialize the next state in the ll */
      tmp_energy = parsed_input[2][j-1];

      /* we need these energies to sort the states later */
      e_vals[l] = tmp_energy;
      set_estate(curr_state, from_last, tmp_idxs, tmp_evals, tmp_tmoms, k,\
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
        fprintf(stderr, "parse_input.c:set_estate_list, the list of electronic \
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
  /* the last state that had transitions found in the output */
  end_state = curr_state -> last;


  /* printf( "%d, %d, %d\n",n_states, l,n_trans); */

  /* printf( "%d\n", curr_state->state_idx); */
  /* printf( "%d\n", (curr_state-> last) -> state_idx); */
  /* printf( "%d\n", end_state -> last -> state_idx); */
  /* printf( "%d\n", end_state -> state_idx); */
  /* printf( "%d\n", curr_inode -> root_e_state -> state_idx); */

  if (SYM == 1) {
    /* loop over the state list and check so that every "from" index is in the
       the list. if not, add it to the end of the list */
    state_bookmark = end_state -> list_idx;
    curr_state = curr_inode -> root_e_state;
    next_state = curr_state;

    while((curr_state = next_state) != NULL && (curr_state -> list_idx) != state_bookmark+1){

      for (j=0; j<curr_state -> n_tfrom; j++) {

        tmp_idx = (curr_state -> idxs_to)[j];
        /* printf( "idx %d\n", tmp_idx); */
        if (get_state_si(curr_inode,tmp_idx) == NULL) {

          /* printf( "found state %d not in list\n", tmp_idx); */
          /* sleep(1); */
          set_estate(end_state -> next, tmp_idx, NULL, NULL, NULL, 0, \
                         parsed_input[2][0],(curr_state -> e_vals)[j]);
          e_vals[l++] = (curr_state -> e_vals)[j];
          end_state = end_state -> next;
        }
      }
      next_state = curr_state -> next;
      /* printf( "next state %d\n", next_state -> list_idx); */
    }

    curr_state = curr_state-> last;

  } else {




    curr_state = end_state;
    /* remove any unused nodes in the llist iterate forwards and free up \
       the unused nodes. */
    /* for (j=0; j<n_states-l; j++) { */
    /*   printf( "freed!\n" ); */
    /*   next_state = curr_state -> next; */
    /*   free(curr_state); */
    /*   curr_state = NULL; /\* destroy the pointer *\/ */
    /*   curr_state = next_state; */
    /* } */
  }
  end_state -> next = NULL;
  n_states = l;
  curr_inode -> n_states = l;

  /* the last node in the llist is the highest energy intermediate state */

  /* for (j=0; j<n_states; j++) { */
  /*   printf( "%le\n", e_vals[j]); */
  /* } */
  if (state_er[0] = 0) {
    /* use the k-means algorithm to do a preliminary sorting of the states */
    k_meansl(e_vals, groups, n_states);
  } else {
    sort_states(state_er, e_vals, groups, n_states);
  }

  /* check so that all states were classified */
  if ((groups[0][0] + groups[1][0]) != n_states) {
    fprintf(stderr, "\n\nError: parse_input.c:set_estate_list, some states\
 from the list of eigenvalues could not be classified as either initial, \
final, or intermediate.\n", l, n_states);
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }
  for (j=0; j<2; j++) {
    for (k=1; k<=groups[j][0]; k++) {
      printf( "groups[%d][%d] = %d\n", j, k, groups[j][k]);
    }
  }
  fprintf(stderr, "\n\n=======Valgrind eject point=======\n\n");
  exit(1);

  /* based on the fact that the transitions between states inside the same
     group will be too low to get included in the data tree, we can sort
     out any state indices that ended up in the wrong category in the k_means
     sorting above.*/
  bw_s = 0;
  for (j=0; j<2; j++) {
    for (k=1; k<=groups[j][0]; k++) {

      sorted_idx = groups[j][k];

      /* check if the sorted_idx state is even in the llist of electronic states */
      if ((curr_state = get_state_li(curr_inode, sorted_idx-1)) != NULL) {
        curr_state -> type = j+1;
        /* store the sum of boltzmann weights to later on use it for the state
           screening process */

        if (j == 0) {

          curr_state -> bw = get_rbdist(parsed_input[2][0],curr_state -> e_val);
          if (curr_state -> state_idx != 1) {
            bw_s += curr_state -> bw;
          }

          (curr_inode -> n_gfs) += 1;

        } else {

          (curr_inode -> n_is) += 1;
        }
        /* printf( "\n\n"); */
        /* /\* count the state types *\/ */
        /* curr_state = curr_inode -> root_e_state; */
        /* while((next_state = curr_state -> next) != NULL){ */
        /*   printf( "state %d, type = %d\n", curr_state -> state_idx, curr_state -> type); */
        /*   curr_state = next_state; */
        /* } */
      }
    }
  }

  if (SYM == 1) {
    set_symtrans(curr_inode);
  }

  curr_inode -> bw_sum = bw_s;
  reset_info_maxvals(curr_inode);

  /* finally, if we want to take symmetric transitions into account, add the right ground
   states to the corresponding final states */

  /* e_statelist2s(curr_inode,1); */

  for (j=0; j<2; j++) {
    free(groups[j]);
  }
  free(groups);

  /* trim off any nodes at the end of the list that hasnt gotten allocated. */
  free(tmp_idxs);
  free(tmp_tmoms);
  free(tmp_evals);
  free(e_vals);

  return EXIT_SUCCESS;
}

int
set_estate (estate st,
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
    fprintf(stderr, "parse_input:function set_estate, malloc: failed \
to allocate memory for \"idxs_to\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  if((tm = malloc(n_trs_from*sizeof(double))) == NULL ){
    fprintf(stderr, "parse_input:function set_estate, malloc: failed \
to allocate memory for \"tmoms\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  if((ev = malloc(n_trs_from*sizeof(double))) == NULL ){
    fprintf(stderr, "parse_input:function set_estate, malloc: failed \
to allocate memory for \"tmoms\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  st -> state_idx = s_idx;
  st -> e_val = e;
  st -> n_tfrom = n_trs_from;

  if (idxs_buf != NULL) {
    for (j=0; j<n_trs_from; j++) {
      idxs[j] = idxs_buf[j];
      tm[j] = moms_buf[j];
      ev[j] = evals_buf[j];
    }
    st -> idxs_to = idxs;
    st -> t_moms = tm;
    st -> e_vals = ev;
    st -> max_tmom = get_maxl(tm, n_trs_from);
  } else {
    st -> idxs_to = NULL;
    st -> t_moms = NULL;
    st -> e_vals = NULL;
    st -> max_tmom = 0;
  }

  return EXIT_SUCCESS;
}

estate
init_estate_list (char * str_id,
               int n_states,
               int n_trans
               ){

  int j,k,l; /* looping variables */
  int from,from_last,to;

  double tmp_energy;
  double tmp_bw;

  estate root_state;
  estate curr_state;
  estate next_state;
  estate last_state;

  info_node curr_inode;

  curr_inode = init_inode(str_id, n_states, n_trans); /* this defines the root node and
                                              assigns it a non-NULL index */

  root_state = malloc(sizeof(struct e_state_s));
  root_state -> last = NULL;
  root_state -> list_idx = 0;
  root_state -> info = curr_inode;

    /* by definition the parsed_input[:][0] should contain information on the
   lowest energy ground state. set the root state node type accordingly. */
  curr_inode -> root_e_state = root_state;
  curr_state = root_state;

  /* construct the list structure */
  for (j=1; j<(n_states+1); j++) {
    /* note that the idxs_to and t_moms arrays dont get allocated here. this
     is done in the set_estate_list function */
    next_state = malloc(sizeof(struct e_state_s));
    next_state -> last = curr_state;
    next_state -> list_idx = j;
    next_state -> info = curr_inode;
    curr_state -> next = next_state;
    curr_state -> type = 0;
    curr_state = next_state;
  }
  curr_state -> next = NULL;

  return root_state;
}

int
is_state_inlist (info_node inode,
                 int idx
                 ) {

  estate next_state = inode -> root_e_state;
  estate curr_state;

  while(next_state != NULL){
    curr_state = next_state;
    if (curr_state -> state_idx == idx) {
      return 1;
    }

    next_state = curr_state -> next;
  }
  return 0;
}

void
swapd_estate(estate e1,
             estate e2){

  e2 -> next = e1 -> next;
  e2 -> last = e1 -> last;
  e2 -> info = e1 -> info;
  e2 -> list_idx = e1 -> list_idx;

  dstruct_estate(e1);
}

void
dstruct_estate(estate e){

  free(e -> idxs_to);
  free(e -> t_moms);
  free(e -> e_vals);
  free(e);
  e = NULL;
}

double
get_trans (estate es, /* root of the electronic state llist */
           int idx_to /* index of the state transitioning to */
           ) {
  /*                if (a1 == NULL) { */
  /*   fprintf(stderr, "Found NULL\n"); */
  /*   printf( "program terminating due to the previous error.\n"); */
  /*   exit(1); */
  /* } */
  int j; /* looping variables */
  int trs_max = es -> n_tfrom;

  int * ti = es -> idxs_to;
  /* printf( "looping\n" ); */
  for (j=0; j<trs_max; j++) {
    /* printf( "in state %d: to = %d, from = %d, %d\n", es -> state_idx,idx_to, ti[j], trs_max); */
    if (ti[j] == idx_to) {
      return (es -> t_moms)[j];
    }
  }
  /* we looped over the entire list of electronic states but didnt find the
   one that was requested */
  fprintf(stderr, "estate.c, function get_trans: unable to locate state of \
index %d in the list of transitions for state %d.\n", idx_to, es -> state_idx);
  e_state2s(es, 1);
  printf( "program terminating due to the previous error.\n");
  exit(1);
}

double
get_eval (estate es, /* root of the electronic state llist */
           int idx_to /* index of the state transitioning to */
           ) {
  int j; /* looping variables */
  int trs_max = es -> n_tfrom;

  int * ti = es -> idxs_to;
  /* printf( "looping\n" ); */
  for (j=0; j<trs_max; j++) {
    /* printf( "to = %d, from = %d, %d\n",idx_to, ti[j], trs_max); */
    if (ti[j] == idx_to) {
      return (es -> e_vals)[j];
    }
  }
  /* we looped over the entire list of electronic states but didnt find the
   one that was requested */
  fprintf(stderr, "estate.c, function get_trans: unable to locate state of \
index %d in the list of transitions for state %d.\n", idx_to, es -> state_idx);
  printf( "program terminating due to the previous error.\n");
  exit(1);
}

estate
get_state_si (info_node inode, /* the info node at root of the state ll */
           int s_idx /* index of the state to get */
           ){

  int n_s = inode -> n_states;
  int n_t = inode -> n_trans;

  estate curr_st = (inode -> root_e_state);
  estate next_st = curr_st;

  while((curr_st = next_st) != NULL){

    if ((curr_st -> state_idx) == s_idx) {
      return curr_st;
    }
    next_st = curr_st -> next;
  }

  /* unable to locate state of in the list of electronic states */
  return NULL;
}

estate
get_state_li (info_node inode, /* the info node at root of the state ll */
              int l_idx /* list index of the state to get */
              ){

  int n_s = inode -> n_states;
  int n_t = inode -> n_trans;

  estate curr_st = (inode -> root_e_state);
  estate next_st = curr_st;

  while((curr_st = next_st) != NULL){

    if ((curr_st -> list_idx) == l_idx) {
      return curr_st;
    }
    next_st = curr_st -> next;
  }

  /* unable to locate state of in the list of electronic states */
  return NULL;
}

double
get_ediff (info_node inode, /* root of the electronic state llist */
           int idx_es1,
           int idx_es2
           ) {
  int j;
  estate es = get_state_si(inode, idx_es1);
  int trs_max = es -> n_tfrom;
  /* printf( "\neval = %le\n", (es -> e_val)); */
  int * ti = es -> idxs_to;
  /* printf( "looping\n" ); */
  for (j=0; j<=trs_max; j++) {
    /* printf( "%le\n", (es -> e_vals)[j]); */
    if (ti[j] == idx_es2) {
      return ((es -> e_vals)[j] - (es -> e_val));
    }
  }
  fprintf(stderr, "estate.c, get_ediff: unable to locate state of \
index %d in the list of transitions from state %d.\n", idx_es2, idx_es1);
  printf( "program terminating due to the previous error.\n");
  exit(1);
}

void
e_statelist2s(info_node inode,
              int flag){

  int j;

  char * fn_in = inode -> str_id;
  int n_s = inode -> n_states;
  int n_t = inode -> n_trans;
  double bw_sum = inode -> bw_sum;
  double gs_eval =  (inode -> root_e_state) -> e_val;
  estate curr_st = (inode -> root_e_state);
  estate next_st = curr_st;
  printf( "\n  -printing the content of the electronic state list:");

  for (j=0; j<n_s; j++) {

    curr_st = next_st;
    if (flag == 1) {
      e_state2s(curr_st,1);
    } else {

      printf( "\n   state[%d/%d], type %d\n", curr_st->list_idx + 1, n_s, curr_st->type);
      printf( "     state_idx = %d\n", curr_st->state_idx);

      if ((curr_st->type) != 2) {
        printf( "     boltzmann weight = %le\n", curr_st->bw);
      }

      printf( "     n_tfrom = %d\n", curr_st->n_tfrom);
      printf( "     e_val(au) = %le\n", curr_st->e_val);
      printf( "     max_tmom = %le\n", curr_st->max_tmom);
      printf( "     delta e(ev) = %le\n", (gs_eval - (curr_st->e_val))*(double)AUTOEV);

    }
    next_st = curr_st -> next;
  }
}

void
  e_state2s(estate es,
            int flag){
  int j;

  printf( "\n   state_idx = %d\n", es->state_idx);
  printf( "     list_idx = %d\n", es->list_idx);
  printf( "     state type %d\n", es->type);

  if ((es->type) != 2) {
    printf( "     boltzmann weight = %le\n", es->bw);
  }

  printf( "     n_tfrom = %d\n", es->n_tfrom);
  printf( "     e_val(au) = %le\n", es->e_val);
  printf( "     max_tmom = %le\n", es->max_tmom);
  if (flag == 1) {
    for (j=0; j<es->n_tfrom; j++) {
      printf( "     -> %d [%d/%d] e = %le, tmom = %le\n", es->idxs_to[j], j+1, j<es->n_tfrom, es->e_vals[j], es->t_moms[j]);
    }
  }
}

void
reset_info_maxvals (info_node inode) {

  double maxt_fs, maxt_is;
  double tmp_max;

  double * trans;

  estate curr_st = (inode -> root_e_state);
  estate next_st = curr_st;
  inode -> mt_is = 0;
  inode -> mt_fs = 0;

  maxt_fs = -0.1;
  maxt_is = -0.1;

  while((curr_st = next_st) != NULL){
    tmp_max = get_maxl(curr_st -> t_moms, curr_st -> n_tfrom);
    curr_st -> max_tmom = tmp_max;
    if ((curr_st -> state_idx) != 2) {
      if (tmp_max >= maxt_fs) {
        maxt_fs = tmp_max;
      }
    }
    else {
      if (tmp_max >= maxt_is) {
        maxt_is = tmp_max;
      }
    }
    next_st = curr_st -> next;
  }

  inode -> mt_is = maxt_is;
  inode -> mt_fs = maxt_fs;

}
