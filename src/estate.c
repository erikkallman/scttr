#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include "k_meansl.h"
#include "estate.h"
#include "structs.h"
#include "info_node.h"
#include "appc.h"
#include "sci_const.h"
#include "std_num_ops.h"
#include "std_char_ops.h"
#define AR_FLAG 0

void
set_ttypes (double * state_er,
            int * mom,
            info_node inode) {

  int n_t;
  double e_val;
  double e_ref;

  /* indices of states stored in.. */
  int ie_idx; /* the info node */
  int ee_idx; /* the estates */

  /* indices stored in the information node */

  estate curr_st = (inode -> root_e_state);
  estate next_st = curr_st;
  estate trs_st;

  while((curr_st = next_st) != NULL){

    /* e_val = ((curr_st -> e_val) - ((inode -> e_states)[(inode->ev_idxs)[1]] -> e_val))*AUTOEV; */

    /* if ((state_er[1] <= e_val) && (state_er[2] >= e_val)) { */
    /*   curr_st -> ttype = 0; */
    /* } */
    /* else if (((state_er[3] <= e_val) && (state_er[4] >= e_val)) && */
    /*          !(e_val >= (state_er[4]+1000))){ */
    /*   curr_st -> ttype = mom[2]; */
    /* } */
    if (curr_st -> type == 1) {
      curr_st -> ttype = 0;
    }
    else {
      curr_st -> ttype = mom[curr_st -> type];
    }

    next_st = curr_st -> next;
  }

  curr_st = (inode -> root_e_state);
  next_st = curr_st;
  e_ref = ((inode -> e_states)[(inode->ev_idxs)[0]] -> e_val);

  while((curr_st = next_st) != NULL){

    /* e_state2s(curr_st,1);  */
    /* sleep(2); */
    n_t = curr_st -> n_tfrom;
    /* printf( "n_t = %d\n", n_t); */
    /* printf( "type = %d\n", curr_st -> type); */

    while(--n_t != -1){
      /* printf( "n_t=%d n_from=%d , idx_to=%d\n", n_t,curr_st -> n_tfrom,(curr_st -> idxs_to)[n_t]); */
      trs_st = get_state_sil(inode,(curr_st -> idxs_to)[n_t]);
      /* printf( "got state state_idx=%d ttype=%d n_tfrom=%d!\n", trs_st -> state_idx,trs_st -> ttype, trs_st -> n_tfrom ); */
      /* e_val = ((curr_st -> e_val) - e_ref)*AUTOEV; */
      /* printf( "got eval\n" ); */
      if (curr_st -> type == 1) {
        /* transition from ground to intermediate state */
        /* printf( "setting type1\n"); */
        (curr_st -> ttypes[n_t]) = mom[trs_st -> type];
        /* printf( "set type1\n" ); */
      }
      else {
        /* transition from intermediate to ground state */
        /* printf( "setting type2\n" ); */
        (curr_st -> ttypes)[n_t] = mom[curr_st -> type];
        /* printf( "set type2\n" ); */
      }
      /* (curr_st -> ttypes[n_t]) = mom[trs_st -> type]; */
    }

    next_st = curr_st -> next;
  }
}
void
sort_states (double * state_er,
             double * e_vals,
             int ** groups,
             int n_states){

  int j,k,l,m; /* looping variables */
  int x,y;
  int lr,hr,tot,last_tot;
  int r_val; /* index to the range value that needs replacing */
  int group_idx;

  double e_val;
  /* message sent to user if values outside of the ranges specified are fund */
  /* printf( "current ranges:\n" ); */
  /* for (j=1; j<=state_er[0]; j++) { */
  /*   printf( "%le,", state_er[j]); */
  /* } */
  /* printf( "\n" ); */

  /* the lowest energy ground state is already added to group 1 */
  /* if (groups[0][0] == 0) { */
  /* groups[0][0] = 2; */
  /* groups[0][1] = 1; */
  /* } */

  /* for (j=0; j<n_states; j++) { */
  /*   printf( "%le\n", e_vals[j]); */
  /*   sleep(1); */
  /* } */
  /* fprintf(stderr, "\n\n=======Valgrind eject point=======\n\n"); */
  /* exit(1); */
  /* for (j=0; j<(int)(state_er[0])/2; j++) { */
  /*   /\* if (groups[j][0] == 0) { *\/ */
  /*   groups[j][0] = 1; */
  /*   /\* } *\/ */
  /* } */

  while(1){
  restart:
    tot = 0;
    /* } */

    for (j=0; j<(int)(state_er[0])/2; j++) {
      /* if (groups[j][0] == 0) { */
      groups[j][0] = 1;
      /* } */
    }

    /* printf( "\n\n\nstarting!\n\n\n"); */
    for (j=0; j<n_states; j++) {
      /* e_val = (e_vals[j] - e_vals[0])*AUTOEV; */
      e_val = e_vals[j];
      last_tot = tot;
      for (k=1; k<(int)state_er[0]; k=k+2) {
        /* l = (k/2) + k % 2; /\* group index from state_er index *\/ */
        l = (k/2); /* group index from state_er index */
        /* printf( "l=%d k=%d\n", l, k); */
        /* printf( "e_val = %le, \n", e_val); */
        /* printf( "\n\n"); */
        /* sleep(1); */
        /* is the energy value within the (k/2)th range? */
        if ((state_er[k] <= e_val) && (state_er[k+1] >= e_val)) {
          /* if (l == 2) { */
          /*   printf( "tot = %d, group index=%d, state_er[%d] = %le to %le\n",tot, l, k, state_er[k], state_er[k+1]); */
          /*   printf( "e_val = %le, \n", e_val); */
          /*   printf( "\n\n"); */
          /*   sleep(1); */
          /* } */
          tot += 1;
          groups[l][groups[l][0]] = tot;
          groups[l][0] += 1;
        }
      }

      /* we didnt find a suiting interval for e_val. can the existing ones be
         extended to include it? */
      if (last_tot == tot) {
        fprintf(stderr, "\n\n=======Automatic range extension engaged=======\n\n");
        printf( "%le\n",e_val );
        fprintf(stderr, "\n\n=======Valgrind eject point=======\n\n");
        exit(1);
        printf( "\n\n\nnot same!\n\n\n" );
        if (e_val < state_er[1]){
          r_val = 1;
          printf( "smaller than small %d\n",r_val );
          /* sleep(1); */
        }
        else if (e_val > state_er[(int)state_er[0]]){
          r_val = state_er[0];
          printf( "larger than large %d\n",r_val );
          /* sleep(1); */
        }
        else if (AR_FLAG == 0){
          printf( "automatique!\n" );
          /* sleep(1); */
          for (k=1; k<=(int)state_er[0]; k=k+2) {
            if (pyth_distl(e_val,state_er[k]) < \
                pyth_distl(e_val,state_er[k+1])) {
              r_val = k;
              /* state_er[k] = e_val; */
            } else {
              r_val = k+1;
              /* state_er[k+1] = e_val; */
            }
          }
          printf( "automatique done! k = %d, r_val = %d\n",k, r_val );
        }
        else {
          printf( "asque!\n" );
          /* manually edit the ranges from user input */
          r_val = send_range_qmsg(state_er,e_val);
        }
        state_er[r_val] = e_val;

        /* if ((r_val == 1) || (r_val == 2)) { */
        /*   group_idx = 0; */
        /* } else { */
        /* p  group_idx = r_val/2-1; */
        /* } */

        /* m = r_val/2; */
        /* /\* printf( "r_val = %d, m = %d\n", r_val, m); *\/ */
        /* /\* add the value to its corresponding group *\/ */
        /* tot += 1; */
        /* groups[group_idx][groups[group_idx][0]] = tot; */
        /* groups[group_idx][0] += 1; */
        /* for (j=0; j<state_er[0]/2; j++) { */
        /*   printf( "total in group %d = %d\n", j, groups[j][0]); */
        /*   for (k=1; k<groups[j][0]; k++) { */
        /*     printf( "groups[%d][%d] = %d\n", j, k, groups[j][k]); */
        /*     groups[j][k] = 0; */
        /*   } */
        /* } */
        /* printf( "\n\n\n" ); */
        /* sleep(1); */
        goto restart;
      }

      /* printf( "\n\n==========\n\n"); */
      /* for (x=0; x<state_er[0]/2; x++) { */
      /*   printf( "\n\ntotal in group %d = %d\n", x, groups[x][0]); */
      /*   for (y=0; y<groups[x][0]; y++) { */
      /*     printf( "  groups[%d][%d] = %d\n", x, y, groups[x][y]); */
      /*   } */
      /* } */
      /* printf( "\n\n==========\n\n"); */
      /* usleep(50000); */
      if ((j+1) == n_states) {
        goto convergence;
      }
    }
  }
 convergence:

  /* for (j=0; j<state_er[0]/2; j++) { */
  /*   printf( "total in group %d = %d\n", j, groups[j][0]); */
  /*   for (k=0; k<groups[j][0]; k++) { */
  /*     printf( "groups[%d][%d] = %d\n", j, k, groups[j][k]); */
  /*   } */
  /* } */
  for (j=0; j<=state_er[0]; j++) {
    printf( "state_er[%d] = %le\n", j, state_er[j]);
  }
}

int
set_estate_list (double * state_er,
                 double ** parsed_input,
                 int * mom,
                 int n_states,
                 int n_trans,
                 char * id) {

  int j,k,l,m; /* looping variables */
  int from_last = parsed_input[0][0];
  int from,to,sorted_idx,tmp_idx,state_bookmark;
  int tmp_ntfrom, tmp_type;
  int n_added = 0;
  int SYM = 0;
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
  int * added_idxs;
  int * idxs_e_vals;
  double * tmp_tmoms;
  double * tmp_evals;
  double * e_vals;


  if((tmp_idxs = malloc(n_states*sizeof(int))) == NULL ){
    fprintf(stderr, "parse_input:function set_estate_list, malloc: failed \
to allocate memory for \"tmp_idxs\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  if((tmp_tmoms = malloc(n_states*sizeof(double))) == NULL ){
    fprintf(stderr, "parse_input:function set_estate_list, malloc: failed \
to allocate memory for \"tmp_tmoms\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  if((tmp_evals = malloc(n_states*sizeof(double))) == NULL ){
    fprintf(stderr, "parse_input:function set_estate_list, malloc: failed \
to allocate memory for \"tmp_evals\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  if((e_vals = malloc(n_states*2*sizeof(double))) == NULL ){
    fprintf(stderr, "parse_input:function set_estate_list, malloc: failed \
to allocate memory for \"e_vals\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }
  if((idxs_e_vals = malloc(n_states*2*sizeof(int))) == NULL ){
    fprintf(stderr, "parse_input:function set_estate_list, malloc: failed \
to allocate memory for \"e_vals\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  if((added_idxs = malloc(n_states*2*sizeof(double))) == NULL ){
    fprintf(stderr, "parse_input:function set_estate_list, malloc: failed \
to allocate memory for \"e_vals\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  if((groups = malloc(mom[0]*sizeof(int*))) == NULL ){
    fprintf(stderr, "parse_input:function set_estate_list, malloc: failed \
to allocate memory for \"groups\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  for (j=0; j<mom[0]; j++) {
    if((groups[j] = malloc((n_states+1)*sizeof(int))) == NULL ){
      fprintf(stderr, "parse_input:function set_estate_list, malloc: failed \
to allocate memory for \"groups\"\n");
      printf( "program terminating due to the previous error.\n");
      exit(1);
    }
    groups[j][0] = 0; /* zero the group counters (see sort_states) */
  }

  curr_inode = get_inode(id);
  curr_state = curr_inode -> root_e_state;

  /* the first node in the llist is the lowest energy ground state */
  curr_state -> type = 1;
  printf( "\n\n=========Constructing the estate list=========\n\n" );
  /* e_vals[0] = parsed_input[2][0]; */
  /* state the state of the previousl initialized estate list */
  for (k=0,j=0,l=1; j<=n_trans; j++) { /* j = read head for the parsed_input matrix */

    if (j < n_trans) {
      from = parsed_input[0][j];
      to = parsed_input[1][j];
    }/*  else { */
    /*   /\* store the data of the final state in the transition list *\/ */
    /*   from = 0; */
    /* } */
    printf( "from = %d to=%d, from_last = %d, %d %d\n", from, to, from_last, j , n_trans);
    if((j == n_trans) || (from != from_last)) { /* we just started reading values for transitions
                                from a new state */

      /* check so that the state is not already in the list */
      if ((tmp_state = get_state_sil(curr_inode,from)) == NULL) {

        /* initialize the next state in the ll */
        tmp_energy = parsed_input[2][j-1];

        /* we need these energies to sort the states later */
        /* e_vals[l] = tmp_energy; */

        set_estate(curr_state, from_last, tmp_idxs, tmp_evals, tmp_tmoms, k,\
                   parsed_input[2][0], tmp_energy);
        next_state = curr_state -> next;
        curr_state = next_state;

        from_last = from;

        /* printf( "\nto=%d, e_val=%le, type=%d\n",(curr_state -> state_idx), \ */
        /*         (curr_state -> e_val), (curr_state -> type)); */
        /* for (m=0; m<k; m++) { */
        /*   printf( "from:%d, ev=%le, mom[%d]=%le\n", (curr_state -> idxs_to)[m],\ */
        /*           (curr_state -> e_vals)[m], m, (curr_state -> t_moms)[m]); */
        /* } */
        /* sleep(1); */
        if (j == n_trans) {
          break;
        }
      } else {
        /* it is in the list */

        tmp_state -> idxs_to = appc_dd((tmp_state -> idxs_to), &tmp_idxs[0] \
                                       , (tmp_state -> n_tfrom), k);

        tmp_state -> t_moms = appc((tmp_state -> t_moms), &tmp_tmoms[0] \
                                    , (tmp_state -> n_tfrom), k);

        tmp_state -> e_vals = appc((tmp_state -> e_vals), &tmp_evals[1]\
                                    , (tmp_state -> n_tfrom), k);

        tmp_state -> n_tfrom += k;

        free(tmp_state -> ttypes);

        tmp_state -> ttypes = malloc((tmp_state -> n_tfrom)*sizeof(int));

      }

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


      k=0;
      j--;
      l++;

    }

    else if((j < n_trans) && (fabs(parsed_input[3][j] - e_vals[0])*AUTOEV < \
                              (state_er[(int)(state_er[0])] + 100 ))){
      /* printf( "%le,%le, %le vs %le\n",parsed_input[3][j],e_vals[0],fabs(parsed_input[3][j] - e_vals[0])*(double)AUTOEV, (state_er[(int)(state_er[0])] + 100 )); */
      /* sleep(1); */


      tmp_idxs[k] = parsed_input[1][j];
      tmp_evals[k] = parsed_input[3][j];
      tmp_tmoms[k] = parsed_input[4][j];
      if (tmp_tmoms[k] > max_tm) {
        max_tm = tmp_tmoms[k];
      }
      if (from == 1) {
        printf( "%d %le %le\n",tmp_idxs[k],tmp_evals[k],tmp_tmoms[k]);
      }
      k++;  /* increase counter for transitions counted */
    }
    /* else{ */
    /*   printf( "%le,%le, %le vs %le\n",parsed_input[3][j],e_vals[0],fabs(parsed_input[3][j] - e_vals[0])*(double)AUTOEV, (state_er[(int)(state_er[0])] + 100 )); */
    /*   sleep(1); */
    /* } */
  }

  /* the last state that had transitions found in the output */
  /* end_state = curr_state -> last; */
  end_state = curr_state; /* -> last; */

  k = 0;
  curr_state = curr_inode -> root_e_state;
  while(curr_state ->state_idx != -1){
    added_idxs[n_added++] = curr_state -> state_idx;
    next_state = curr_state -> next;
    curr_state = next_state;
  }
  /* for (j=0; j<n_added; j++) { */
  /*   printf( "added_idxs[%d/%d] = %d\n",j,n_added,added_idxs[j]); */
  /* } */

  /* e_state2s(get_state_sil(curr_inode,1),1); */

  /* for (j=0; j<l; j++) { */
  /*   printf( "%le\n", e_vals[j]); */
  /* } */

  /* printf( "%d, %d, %d\n",n_states, l,n_trans); */

  /* printf( "%d\n", curr_state->state_idx); */
  /* printf( "%d\n", curr_state->list_idx); */
  /* printf( "%d\n", (curr_state-> last) -> state_idx); */
  /* printf( "%d\n", (curr_state-> next) -> state_idx); */
  /* printf( "%d\n", end_state -> last -> state_idx); */
  /* printf( "%d\n", end_state -> state_idx); */
  /* printf( "%d\n", curr_inode -> root_e_state -> state_idx); */
  /* printf( "%d\n", curr_inode -> root_e_state -> list_idx); */

  k = 0;

  if (SYM == 0) {
    printf( "\n\n=========Adding from states=========\n\n" );
    /* loop over the state list and check so that every "from" index is in the
       the list. if not, add it to the end of the list */
    state_bookmark = end_state -> list_idx;
    curr_state = curr_inode -> root_e_state;
    next_state = curr_state;

    while((curr_state != NULL) && ((curr_state -> list_idx) != state_bookmark)){
      /* printf( "curr_idx = %d != state_bookmark = %d\n", curr_state -> list_idx,state_bookmark); */
      for (j=0; j<curr_state -> n_tfrom; j++) {

        tmp_idx = (curr_state -> idxs_to)[j];
        /* printf( "idx %d\n", tmp_idx);  */
        if ((get_state_sil(curr_inode,tmp_idx) == NULL) &&\
            !intinint(added_idxs,tmp_idx,n_added) &&\
            (fabs((curr_state -> e_vals)[j] - e_vals[0])*AUTOEV < \
             (state_er[(int)(state_er[0])] + 100 ))) {

          /* printf( "found state %d not in list, n_states=%d, %le\n", tmp_idx, curr_inode -> n_states,((curr_state -> e_vals)[j] - e_vals[0])*AUTOEV); */
          /* tmp_idx is not in the estate list of curr_inode. find the
             first state that hasnt had its state_idx set */
          tmp_state = end_state;
          while((tmp_state -> next != NULL) && ((tmp_state -> state_idx) != -1)){
            /* printf( "idx = %d\n", tmp_idx); */
            tmp_state = tmp_state -> next;
          };
          /* e_statelist2s(curr_inode, 0); */
          end_state = tmp_state;
          /* check the outcome of the above while loop */
          if ((end_state -> state_idx) != -1) {
            /* end_state = end_state -> last; */
            /* printf( "\n\n========= NULL =========\n\n" ); */
            /* printf( "initializing the state!\n" ); */
            inita_estate(end_state,id);
            /* printf( "initialized the state!\n" ); */
          } else {
            /* printf( "\n\n========= is -1 =========\n\n" ); */
          }
          /* printf( "setting the state!\n" ); */
          set_estate(end_state, tmp_idx, NULL, NULL, NULL, 0, \
                     parsed_input[2][0],(curr_state -> e_vals)[j]);
          /* printf( "set the state!\n" ); */
          if (get_state_sil(curr_inode,tmp_idx) == NULL) {
            fprintf(stderr, "estate.c, function set_estate_list: state %d at %d\
 not added properly to the estate list.\n",tmp_idx, end_state -> list_idx);
            printf( "program terminating due to the previous error.\n");
            exit(1);
          } else {
            /* e_vals[l++] = (curr_state -> e_vals)[j]; */
            added_idxs[n_added++] = tmp_idx;
            /* e_statelist2s(curr_inode, 0); */

          }
        }
      }

      next_state = curr_state -> next;
      curr_state = next_state;
      /* printf( "next state %d\n", next_state -> list_idx); */
    }
    curr_state = curr_state-> last;
    printf( "\n\n=========Done adding from states=========\n\n" );
    /* e_statelist2s(curr_inode, 0); */

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
  n_states = curr_inode -> n_states;
  /* e_state2s(get_state_sil(curr_inode,1),1); */
  /* e_state2s(get_state_sil(curr_inode,619),0); */
  /* e_state2s(get_state_sil(curr_inode,620),0); */

  k = 0;
  curr_state = curr_inode -> root_e_state;
  while(curr_state != NULL){
    idxs_e_vals[k] = curr_state -> state_idx;
    e_vals[k] = ((curr_state -> e_val) - parsed_input[2][0])*AUTOEV;
    /* printf( "%le\n", ((curr_state -> e_val) - parsed_input[2][0])*AUTOEV); */
    /* printf( "%le\n\n", ((curr_state -> e_val) - parsed_input[2][0])*AUTOEV); */
    /* sleep(1); */
    k++;
    next_state = curr_state -> next;
    curr_state = next_state;
  }

  /* now that all states are in the elist, set the ev_idxs array */
  printf( "\n\n=========Setting the ev_idxs array=========\n\n" );
  for (j=0; j<n_states; j++) {
    /* printf( "j=%d si=%d\n", j, (curr_inode -> e_states)[j] -> state_idx); */
    (curr_inode -> ev_idxs)[(curr_inode -> e_states)[j] -> state_idx] = j;
  }
  printf( "\n\n=========Done setting the ev_idxs array=========\n\n" );

  /* the last node in the llist is the highest energy intermediate state */

  if ((int)state_er[0] == 0) {
    /* use the k-means algorithm to do a preliminary sorting of the states */
    k_meansl(e_vals, groups, n_states);
  } else {
    printf( "\n\n=========Sorting states=========\n\n" );
    sort_states(state_er, e_vals, groups, n_states);
    printf( "\n\n=========Done sorting states=========\n\n" );
  }

  /* for (j=0; j<mom[0]; j++) { */
  /*   printf( "total in group %d = %d\n", j, groups[j][0]); */
  /*   for (k=0; k<groups[j][0]; k++) { */
  /*     printf( "groups[%d][%d] = %d = %d\n", j, k, groups[j][k],idxs_e_vals[groups[j][k]-1]); */
  /*   } */
  /* } */
  /* printf( "\n\n" ); */
  /* printf( "current ranges:\n" ); */
  /* for (j=1; j<=state_er[0]; j++) { */
  /*   printf( "%le,", state_er[j]); */
  /* } */
  /* printf( "\n" ); */

  /* fprintf(stderr, "\n\n=======Valgrind eject point=======\n\n"); */
  /* exit(1); */
  /* check so that all states were classified */

  for (k=0,j=0; j<mom[0]; j++) {
    k += groups[j][0]-1;
  }

  if (k != n_states) {
    fprintf(stderr, "\n\nError: parse_input.c:set_estate_list, some states\
 from the list of eigenvalues could not be classified as either initial, \
final, or intermediate. counted %d, expected %d.\n", k, n_states);
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  /* based on the fact that the transitions between states inside the same
     group will be too low to get included in the data tree, we can sort
     out any state indices that ended up in the wrong category in the k_means
     sorting above.*/
  bw_s = 0;
  for (j=0; j<mom[0]; j++) {
    for (k=1; k<groups[j][0]; k++) {
      /* printf( "groups[%d][%d] = %d = %d\n", j, k, groups[j][k],idxs_e_vals[groups[j][k]-1]); */
      sorted_idx = idxs_e_vals[groups[j][k]-1];

      /* check if the sorted_idx state is even in the llist of electronic states */
      if ((curr_state = get_state_sil(curr_inode, sorted_idx)) != NULL) {
        curr_state -> type = j+1;
        /* store the sum of boltzmann weights to later on use it for the state
           screening process */

        if (j == 0) {

          curr_state -> bw = get_rbdist(parsed_input[2][0],curr_state -> e_val);
          /* curr_state -> bw = get_bdist(curr_state -> e_val); */
          /* if (curr_state -> state_idx != 1) { */
          bw_s += curr_state -> bw;
          /* } */

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
  /* e_state2s(get_state_sil(curr_inode,609),0); */
  /* e_state2s(get_state_sil(curr_inode,610),0); */
  /* e_state2s(get_state_sil(curr_inode,611),0); */
  /* e_state2s(get_state_sil(curr_inode,40),0); */
  /* e_state2s(get_state_sil(curr_inode,41),0); */
  /* fprintf(stderr, "\n\n=======Valgrind eject point=======\n\n"); */
  /* exit(1); */
    /* e_state2s(get_state_si(curr_inode,5),1); */
  /* e_statelist2s(curr_inode,1); */
  /* if we want to take symmetric transitions into account, add the right ground
   states to the corresponding final states */
  printf( "\n\n=========Setting symtrans=========\n\n" );
  if (SYM == 1) {
    set_symtrans(curr_inode);
  }

  printf( "\n\n=========Done setting symtrans=========\n\n" );
  curr_inode -> bw_sum = bw_s;

  /* if more than one range was provided in the CLI, this means the user
   wants to read more than one range of intermediate states (dipole transitions for example) */


  /* e_state2s(get_state_si(curr_inode,5),1); */

  set_ttypes(state_er,mom,curr_inode);

  reset_info_maxvals(curr_inode);


  /* e_state2s(get_state_sil(curr_inode,609),0); */
  /* e_state2s(get_state_sil(curr_inode,610),0); */
  /* e_state2s(get_state_sil(curr_inode,611),0); */
  /* e_state2s(get_state_sil(curr_inode,40),0); */
  /* e_state2s(get_state_sil(curr_inode,41),0); */
  /* e_state2s(get_state_sil(curr_inode,1),1); */
  /* fprintf(stderr, "\n\n=======Valgrind eject point=======\n\n"); */
  /* exit(1); */

  /* e_statelist2s(curr_inode,1); */

  /* e_state2s(get_state_si(curr_inode,620),1); */

  for (j=0; j<mom[0]; j++) {
    free(groups[j]);
  }
  free(groups);

  /* trim off any nodes at the end of the list that hasnt gotten allocated. */

  free(tmp_idxs);

  /* free(tmp_tmoms); */

  /* free(tmp_evals); */

  free(e_vals);
  free(idxs_e_vals);

  return EXIT_SUCCESS;
}

int
set_estate(estate st,
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
  int * tt;
  double * tm;
  double * ev;
  if (idxs_buf != NULL) {
    if((idxs = malloc(n_trs_from*sizeof(int))) == NULL ){
      fprintf(stderr, "parse_input:function set_estate, malloc: failed \
to allocate memory for \"idxs_to\"\n");
      printf( "program terminating due to the previous error.\n");
      exit(1);
    }

    if((tt = malloc(n_trs_from*sizeof(int))) == NULL ){
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
    st -> ttypes = tt;
    st -> idxs_to = idxs;
    st -> t_moms = tm;
    st -> e_vals = ev;
    st -> max_tmom = get_maxl(tm, n_trs_from);
  } else {
    st -> ttypes = NULL;
    st -> idxs_to = NULL;
    st -> t_moms = NULL;
    st -> e_vals = NULL;
    st -> max_tmom = 0;
  }

  return EXIT_SUCCESS;
}

estate
inita_estate (estate curr_state,
              char * id
              ){

  int j,k,l; /* looping variables */
  int from,from_last,to;
  int n_states;

  estate next_state;

  info_node curr_inode;
  /* printf( "0\n" ); */
  curr_inode = get_inode(id);
  /* printf( "0.1\n" ); */
  n_states = curr_inode -> n_states;
  /* printf( "0.2\n" ); */
  (curr_inode -> e_states)[n_states] = curr_state;
  /* printf( "0.5\n" ); */
  /* note that the idxs_to and t_moms arrays dont get allocated here. this
     is done in the set_estate function */
  next_state = malloc(sizeof(struct e_state_s));
    /* printf( "1\n" ); */
  next_state -> max_tmom = -1;
  next_state -> last = curr_state;
  next_state -> next = NULL;
  next_state -> list_idx = n_states;
  next_state -> info = curr_inode;
    /* printf( "2\n" ); */
  curr_state -> next = next_state;
  curr_state -> type = 0;
  curr_state = next_state;
  /* printf( "allocated!\n" ); */
  curr_inode -> n_states += 1;
  return curr_state;
}

estate
init_estate_list (char * str_id,
                  int * mom,
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

  curr_inode = init_inode(str_id, mom, n_states, n_trans); /* this defines the root node and
                                              assigns it a non-NULL index */

  root_state = malloc(sizeof(struct e_state_s));
  root_state -> last = NULL;
  root_state -> list_idx = 0;
  root_state -> max_tmom = -1;
  root_state -> info = curr_inode;

    /* by definition the parsed_input[:][0] should contain information on the
   lowest energy ground state. set the root state node type accordingly. */
  curr_inode -> root_e_state = root_state;
  curr_state = root_state;
  for (j=0; j<n_states; j++) {
    (curr_inode -> ev_idxs)[j] = 0;
  }
  /* construct the list structure */
  for (j=1; j<(n_states+1); j++) {
    (curr_inode -> e_states)[j-1] = curr_state;

    /* note that the idxs_to and t_moms arrays dont get allocated here. this
     is done in the set_estate_list function */
    next_state = malloc(sizeof(struct e_state_s));
    next_state -> last = curr_state;
    next_state -> max_tmom = -1;
    next_state -> list_idx = j;
    next_state -> state_idx = -1; /* sign of state idx not having been set */
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
  for (j=0; j<=trs_max; j++) {
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
  /* int j; */
  /* int n_s = inode -> n_states; */
  /* estate * sl = inode -> e_states; */

  /* for (j=0; j<n_s; j++) { */
  /*   if ((sl[j] -> state_idx) == s_idx) { */
  /*     return sl[j]; */
  /*   } */
  /* } */

  /* int n_s = inode -> n_states; */
  /* estate * sl = inode -> e_states; */

  /* while(--n_s != -1){ */
  /*   if ((sl[n_s] -> state_idx) == s_idx) { */
  /*     return sl[n_s]; */
  /*   } */
  /* } */
  /* printf( "si1:%d evi:%d si2:%d\n", s_idx, inode -> ev_idxs[s_idx], (inode->e_states)[inode -> ev_idxs[s_idx]] -> state_idx); */
  return (inode->e_states)[inode -> ev_idxs[s_idx]];

  /* estate curr_st = (inode -> root_e_state); */
  /* estate next_st = curr_st; */

  /* while((curr_st = next_st) != NULL){ */

  /*   if ((curr_st -> state_idx) == s_idx) { */
  /*     return curr_st; */
  /*   } */
  /*   next_st = curr_st -> next; */
  /* } */

  /* /\* unable to locate state of in the list of electronic states *\/ */
  /* return NULL; */
}

estate
get_state_sil (info_node inode, /* the info node at root of the state ll */
              int s_idx /* index of the state to get */
              ){
  /* int j; */
  /* int n_s = inode -> n_states; */
  /* estate * sl = inode -> e_states; */

  /* for (j=0; j<n_s; j++) { */
  /*   if ((sl[j] -> state_idx) == s_idx) { */
  /*     return sl[j]; */
  /*   } */
  /* } */

  /* int n_s = inode -> n_states; */
  /* estate * sl = inode -> e_states; */

  /* while(--n_s != -1){ */
  /*   if ((sl[n_s] -> state_idx) == s_idx) { */
  /*     return sl[n_s]; */
  /*   } */
  /* } */
  /* printf( "si1:%d evi:%d si2:%d\n", s_idx, inode -> ev_idxs[s_idx], (inode->e_states)[inode -> ev_idxs[s_idx]] -> state_idx); */
  /* return (inode->e_states)[inode -> ev_idxs[s_idx]]; */

  estate curr_st = (inode -> root_e_state);
  estate next_st = curr_st;

  while((curr_st = next_st) != NULL){

    if ((curr_st -> state_idx) == s_idx) {
      return curr_st;
    }
    next_st = curr_st -> next;
  }

  /* /\* unable to locate state of in the list of electronic states *\/ */
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

/* double */
/* get_ediff (info_node inode, /\* root of the electronic state llist *\/ */
/*            int idx_es1, */
/*            int idx_es2 */
/*            ) { */
/*   int j; */
/*   estate es = get_state_si(inode, idx_es1); */
/*   int trs_max = es -> n_tfrom; */
/*   /\* printf( "\neval = %le\n", (es -> e_val)); *\/ */
/*   int * ti = es -> idxs_to; */
/*   /\* printf( "looping\n" ); *\/ */
/*   for (j=0; j<=trs_max; j++) { */
/*     /\* printf( "%le\n", (es -> e_vals)[j]); *\/ */
/*     if (ti[j] == idx_es2) { */
/*       return ((es -> e_vals)[j] - (es -> e_val)); */
/*     } */
/*   } */
/*   fprintf(stderr, "estate.c, get_ediff: unable to locate state of \ */
/* index %d in the list of transitions from state %d.\n", idx_es2, idx_es1); */
/*   printf( "program terminating due to the previous error.\n"); */
/*   exit(1); */
/* } */

double
get_ediff (estate es1,
           estate es2
           ){
  return (es2->e_val)-(es1->e_val);
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
      printf( "     list_idx = [%d/%d]\n", curr_st->list_idx, inode -> n_states);
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
  if (es == NULL) {
    fprintf(stderr, "estate.c, function e_state2s: input state not found in the list of electronic states.\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  int j;
  info_node inode = es->info;

  printf( "\n   state_idx = %d\n", es->state_idx);
  printf( "     list_idx = %d\n", es->list_idx);
  printf( "     state type %d\n", es->type);
  printf( "     state ttype %d\n", es->ttype);
  if ((es->type) != 2) {
    printf( "     boltzmann weight = %le\n", es->bw);
  }

  printf( "     n_tfrom = %d\n", es->n_tfrom);
  printf( "     e_val(au) = %le\n", es->e_val);
  printf( "     max_tmom = %le\n", es->max_tmom);
  if (flag == 1) {
    for (j=0; j<es->n_tfrom; j++) {
      printf( "     -> %d [%d/%d] e = %le, tmom = %le, type=%d, ttype=%d\n", es->idxs_to[j], j+1, j<es->n_tfrom, es->e_vals[j], es->t_moms[j],
              /* inode->e_states[inode->ev_idxs[es->idxs_to[j]]]->type */1, es->ttypes[j]);
    }
  }
}

void
reset_info_maxvals (info_node inode) {

  int j,k,l;
  int from_type;

  double * trans;
  double * tmp_max;

  estate curr_st = (inode -> root_e_state);
  estate next_st = curr_st;
  /* for (l=0; l<(curr_st -> n_tfrom); l++) { */
  /*   printf( "%d %d %d\n",(curr_st -> type),(curr_st -> ttypes)[l],(curr_st -> idxs_to)[l] ); */
  /* } */

  while((curr_st = next_st) != NULL){
    /* printf( "%d\n",(curr_st -> type)); */
    curr_st -> max_tmom = -1;
    for (j=0; j<(curr_st -> n_tfrom); j++) {

      /* printf( "storing temp\n" ); */
      from_type = (inode -> e_states[inode -> ev_idxs[curr_st -> idxs_to[j]]]) -> type;
      /* from_type = curr_st -> type; */
      tmp_max = &(inode -> mt)[(curr_st -> ttypes[j])-1][from_type-1];
      /* printf( "%d %d\n",from_type-1, (curr_st -> ttypes[j])-1); */
      /* printf( "stored temp\n" ); */
      /* printf( "from %d to %d, type from =%d, type to =%d, to ttype=%d\n",curr_st -> state_idx, curr_st -> idxs_to[j], curr_st ->type, from_type, (curr_st -> ttypes)[j] ); */
      if ((curr_st-> t_moms)[j] > *tmp_max) {
        /* set info node specific transition max */
        *tmp_max = (curr_st-> t_moms)[j];
        /* printf( "new max1! %le\n", *tmp_max); */
        /* printf( "new max2!%d,%d %le\n",from_type-1,(curr_st -> ttypes[j])-1, (inode -> mt)[from_type-1][(curr_st -> ttypes[j])-1]); */
      }
      if ((curr_st-> t_moms)[j] > curr_st -> max_tmom) {
        /* set state-specific transition max */
        curr_st -> max_tmom = (curr_st-> t_moms)[j];
        /* printf( "new local max! %le\n", curr_st -> max_tmom); */
      }
      /* printf( "\n\n" ); */
    }
    /* printf( "done analyzing\n" ); */
    /* fprintf(stderr, "\n\n=======Valgrind eject point=======\n\n"); */
    /* exit(1); */
    /* printf( "\n\n" ); */
    /* for (j=0; j<(inode -> mom_types[0])-1; j++) { */
    /*   printf( "%le %le\n", (inode -> mt)[0][j],(inode -> mt)[1][j]); */
    /* } */
    /* printf( "\n\n" ); */
    next_st = curr_st -> next;
  }
  (inode -> mt)[0][2] = 5.993278e+19;
  for (j=0; j<(inode -> mom_types[0]); j++) {

    printf( "%le %le\n", (inode -> mt)[0][j],(inode -> mt)[1][j]);
  }
}
