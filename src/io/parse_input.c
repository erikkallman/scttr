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
#include "get_numsl.h"
#include "parse_input.h"
#include "input_formats.h"
#include "sci_const.h"
#include "dynarray.h"

#define BUF_SIZE 256

int n_info_nodes; /* counter for number of  */
double * input_data[4];
int ** state_indices;

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

int **
screen_states (char * fn_infile,
             int n_args,
             ...){
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
  int ** igs; /* .. ground states */
  int ** tmp_iis; /* .. intermediate states */
  int * tmp_ifs; /* .. final states */

  /* some hard-coded threshold values that will get replaced by the call
   values */
  float thrsh_vals[3] = {0.5, 0.2, 0.2};
  float tmp_thrsh;

  double * tmp_trans;

  n_states = inode -> n_states;

  e_state root_state = inode -> root_e_state;
  e_state curr_state = root_state;
  e_state next_state;
  e_state gs_bm; /* ground state bookmark */
  e_state is_bm; /* intermediate state bookmark */

  va_list argv;
  va_start(argv, n_args);

  mdda_s igs1;
  mdda_init(&igs1);
  mdda_set(&igs1, 42, 42, 42);

  printf( "got:%d\n", mdda_get(&igs1, 42, 42));
  mdda_free(&igs1);

  fprintf(stderr, "\n\n=======Valgrind eject point=======\n\n");
  exit(1);

  if((igs = malloc((n_gs+1)*sizeof(int*))) == NULL ){
    fprintf(stderr, "parse_input:function screen_states, malloc: failed \
to allocate memory for \"igs\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  if((igs[0] = malloc((n_gs+1)*sizeof(int))) == NULL ){
    fprintf(stderr, "parse_input:function screen_states, malloc: failed \
to allocate memory for \"igs[0]\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  for (j=0; j<n_args; j++) {
    tmp_thrsh = (float)va_arg(argv, double); /* grab the next vararg */
    thrsh_vals[j] = tmp_thrsh;
  }
  /* for a description of the screening algorithm below, look up this function
     in the parse_input.h file */

  n_sgs = 0;
  while((next_state = curr_state -> next) != NULL){

    if ((curr_state -> type) == 1) { /* if we found a ground state */

      /* perform stage 1 of the screening process */
      if (((curr_state ->  bw)/(inode -> bw_sum)) >= thrsh_vals[0]) {
        /* we found a ground state with high enough boltzmann weight */
        gs_idx = curr_state -> state_idx;
        n_sgs++;
        igs[0][0] = n_sgs;
        igs[0][n_sgs] = gs_idx;

        printf( "gs[%d] = %d %le\n", n_sgs,gs_idx,((curr_state ->  bw)));
        sleep(1);
        /* store a bookmark of where we found the last ground state */
        gs_bm = next_state;

        /* perform stage 2 of the screening process */
        if((tmp_iis = malloc((n_is+1)*sizeof(int*))) == NULL ){
          fprintf(stderr, "parse_input:function screen_states, malloc: failed \
to allocate memory for \"iis\"\n");
          printf( "program terminating due to the previous error.\n");
          exit(1);
        }

        if((tmp_iis[0] = malloc((n_is+1)*sizeof(int))) == NULL ){
          fprintf(stderr, "parse_input:function screen_states, malloc: failed \
to allocate memory for \"iis[0]\"\n");
          printf( "program terminating due to the previous error.\n");
          exit(1);
        }

        igs[n_sgs] = *tmp_iis; /* store a pointer to the table of IS idices */
        n_sis = 0;
        tmp_trans = curr_state -> t_moms;
        tmp_idxs = curr_state -> idxs_to;
        for (j=0; j<((curr_state -> n_tfrom)-1); j++) {
          printf( "IS:j = %d, sidx = %d, n_istrans = %d\n", j, tmp_idxs[j], (curr_state -> n_tfrom));
          /* if IS transition has a transition moment above threshold */
          if ((tmp_trans[j]/(inode -> max_tmom)) >= thrsh_vals[1]) {
            is_bm = curr_state;
            is_idx = tmp_idxs[j];

            /* locate the candidate intermediate state in the list and varify
               that the state is indeed an IS (has type == 2) */
            curr_state = root_state;
            while((next_state = curr_state -> next) != NULL){
              if ((curr_state -> state_idx) == is_idx) { /* find the state in the llist */

                /* varify that it has been marked as an intermediate state */
                if ((curr_state -> type) != 2){
                  printf( "wrong type..\n" );
                  break; /* it wasnt. move on to the next transition */
                }
                printf( "found is idxs=%d\n", is_idx);

                /* from the ground state found above, we also found an intermediate
                   state with high enough relative transition moment */
                n_sis++;
                tmp_iis[0][0] = n_sis;
                tmp_iis[0][n_sis] = is_idx;

                printf( "  is[%d] = %d %le\n", n_sis, tmp_iis[0][n_sis], (curr_state -> t_moms)[j]/(inode -> max_tmom));

                if((tmp_ifs = malloc((n_fs+1)*sizeof(int))) == NULL ){
                  fprintf(stderr, "parse_input:function screen_states, malloc: failed \
to allocate memory for \"tmp_ifs\"\n");
                  printf( "program terminating due to the previous error.\n");
                  exit(1);
                }
                printf( "allocated\n" );
                /* array for storing all final states reachable from the screened
                   intermediate state */
                tmp_iis[n_sis] = tmp_ifs;
                n_sfs = 0;

                /* perform stage 3 of the screening process */
                /* find final state transitions that pass stage 3 screening */
                for (k=0; k < (curr_state -> n_tfrom); k++) {
                  if ((((curr_state -> t_moms)[k])/(inode -> max_tmom)) >= thrsh_vals[2]) {
                    /* from the intermediate state found above, we also found a
                       final state transition with high enough relative
                       transition moment */
                    n_sfs++;
                    tmp_ifs[0] = n_sfs;
                    tmp_ifs[n_sfs] = (curr_state -> idxs_to)[k];
                    printf( "    fs[%d] = %d %le\n", n_sfs, tmp_ifs[n_sfs], (((curr_state -> t_moms)[k])/(inode -> max_tmom)));
                  }
                }
                printf( "scanned all final. breaking!\n" );
                break; /* the IS has been processed , but is this break needed?*/
              }
              curr_state = next_state;
            }
            curr_state = is_bm; /* reset the state to the last initial state */
          }
        }
        /* jump back to the last ground state that was found in the list */
        printf( "broke out!\n");
        next_state = gs_bm;
      } else {
      printf( "fail! gs[%d] = %d %le\n", n_sgs,gs_idx,((curr_state ->  bw)));
      }
    }
    curr_state = next_state;
  }

 return igs;
}

double **
reduce_input (int ** sidxs /* state indices */
              ){
  int j; /* looping variables */

  /* placeholder allocation to test out the function structure */
  double ** ret = malloc(1*sizeof(double*));

  ret[0] = malloc(sizeof(double));

  /* for (j=0; j<3; j++) { */
  /*   free(sidxs[j]); */
  /* } */
  /* free(sidxs); */

  return ret;
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

  if((new_info_node = (info_node)malloc(sizeof(struct info_node_s))) == NULL ){
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

  root_state = (e_state)malloc(sizeof(struct e_state_s));
  root_state -> last = NULL;
  root_state -> next = NULL;
  root_state -> list_idx = 0;
  root_state -> info = curr_info_node;

    /* by definition the parsed_input[:][0] should contain information on the
   lowest energy ground state. set the root state node type accordingly. */
  curr_info_node -> root_e_state = root_state;
  curr_state = root_state;

  /* construct the list structure */
  for (j=1; j<n_states; j++) {
    /* note that the idxs_to and t_moms arrays dont get allocated here. this
     is done in the set_state_ll function */
    next_state = (e_state)malloc(sizeof(struct e_state_s));
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
  st -> bw = exp(-(e - e_rel)*AUTOEV/(TEXP*TTOEV))/2;
  st -> e_val = e;
  st -> n_tfrom = n_trs_from;
  st -> idxs_to = idxs;
  st -> t_moms = tm;
  st -> e_vals = ev;

  return EXIT_SUCCESS;
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
  int * g;
  int * ito;
  int ** groups;

  double tmp_energy;
  double tmp_sum; /* sum of all boltzmann weights for the electronic states
                   in the system */
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

  if((groups = malloc(3*sizeof(int*))) == NULL ){
    fprintf(stderr, "parse_input:function init_state_ll, malloc: failed \
to allocate memory for \"groups\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  for (j=0; j<3; j++) {
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

    from = parsed_input[0][j];
    to = parsed_input[1][j];

    if (from != from_last) { /* we just started reading values for transitions
                                from a new state */

      /* initialize the next state in the ll */
      tmp_energy = parsed_input[2][j-1];

      /* we need these energies to sort the states later */
      e_vals[l++] = tmp_energy;

      set_state_node(curr_state, from_last, tmp_idxs, tmp_evals, tmp_tmoms, k,\
                     parsed_input[2][0], tmp_energy);
      tmp_sum += curr_state -> bw;
      /* printf( "\nto=%d, e_val=%le, type=%d\n",(curr_state -> state_idx), \ */
      /*         (curr_state -> e_val), (curr_state -> type)); */
      /* for (m=0; m<k; m++) { */
      /*   printf( "from:%d, ev=%le, mom[%d]=%le\n", (curr_state -> idxs_to)[m],\ */
      /*           (curr_state -> e_vals)[m], m, (curr_state -> t_moms)[m]); */
      /* } */
      /* sleep(1); */
      if (((curr_state -> next) == NULL) && ((l+1) <= n_states)) {
        fprintf(stderr, "parse_input.c:set_state_ll, the list of electronic \
states ended before all input data was transfered. there are %d states \
left to process.\n",l);
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
    }
    else {
      tmp_idxs[k] = parsed_input[1][j];
      tmp_evals[k] = parsed_input[3][j];
      tmp_tmoms[k] = parsed_input[4][j];
      if (tmp_tmoms[k] > max_tm) {
        max_tm = tmp_tmoms[k];
      }
      k++;  /* increase counter for transitions counted */
    }
  }
  curr_info_node -> max_tmom = max_tm;
  end_state = curr_state;
  /* remove any unused nodes in the llist iterate forwards and free up\
     the unused nodes. */

  for (j=0; j<l-n_states; j++) {
    printf( "freed!\n" );
    next_state = curr_state -> next;
    free(curr_state);
    curr_state = next_state;
  }

  end_state -> next = NULL;

  /* the last node in the llist is the highest energy intermediate state */
  end_state -> type = 3;

  curr_info_node -> bw_sum = tmp_sum;
  curr_info_node -> n_states = n_states - l;

  /* use the k-means algorithm to do a preliminary sorting of the states */
  k_meansl(e_vals, groups, l);

  end_state = curr_state -> last;

  /* based on the fact that the transitions between states inside the same
     group will be too low to get included in the data tree, we can sort
     out any state indices that ended up in the wrong category in the k_means
     sorting above.*/
  for (j=0; j<3; j++) {
    for (k=1; k<=groups[j][0]; k++) {

      sorted_idx = groups[j][k];
      g = &groups[j][1];
      curr_state = curr_info_node -> root_e_state;

      while((next_state = curr_state -> next) != NULL){

        to = (curr_state -> state_idx);

        if (to == sorted_idx) {
          /* store a pointer to the "to" state idxs_to array and check so that every foundstate is not inside of it as well */
          ito = (curr_state -> idxs_to);
          tmp_ntfrom = curr_state -> n_tfrom;
          tmp_type = curr_state -> type;
          tmp_state = curr_state;

          /* loop over the llist once more and find all states that do not have
             any transitions to the selected state */
          curr_state = curr_info_node -> root_e_state;
          while((next_state = curr_state -> next) != NULL){
            tmp_idx = curr_state -> state_idx;

            /* check if the state_idxs is even in the current group */
            if (intinint(g,tmp_idx, groups[j][0]) == 1){

              if ((tmp_idx != to) && (tmp_idx != 1)){
                /* finally check if there are any transitions from the current
                   state to the selected "to" index, and vice versa */
                if ((intinint(curr_state -> idxs_to, to, curr_state -> n_tfrom) \
                     == 0) && (intinint(ito, tmp_idx, curr_state -> n_tfrom) \
                               == 0)){

                  if (tmp_type == 0) {
                    curr_state -> type = 2;
                    if (j == 2) {
                      curr_state -> type = 3;
                      tmp_state -> type = 3;
                    }
                  } else if ((j == 0) && (tmp_type == 1)) {
                    curr_state -> type = 1;
                  } else if ((tmp_type == 2) && (curr_state -> type == 0)) {
                    curr_state -> type = 2;
                  }
                }
              }
            }
            curr_state = next_state;
          }
          break; /* the state of index "to" has been sorted. */
        }
        if ((curr_state -> next) == NULL) {
          fprintf(stderr, "parse_input.c, function set_state_ll: the state of \
index %d could not be found in the linked list of states. error in state type\
 is expected as a consequence of this.\n",to);
          printf( "program terminating due to the previous error.\n");
          exit(EXIT_FAILURE);
        }
        curr_state = next_state;
      }
    }
    curr_state = curr_info_node -> root_e_state;
    /* printf( "\n\n"); */
    /* /\* count the state types *\/ */
    /* while((next_state = curr_state -> next) != NULL){ */
    /*   printf( "state %d, type = %d\n", curr_state -> state_idx, curr_state -> type); */
    /*   curr_state = next_state; */
    /* } */
  }

  /* count the state types */
  curr_state = curr_info_node -> root_e_state;
  while((next_state = curr_state -> next) != NULL){
    if ((curr_state -> type) == 1) {
      (curr_info_node -> n_gs) += 1;
    }
    else if ((curr_state -> type) == 2) {
      (curr_info_node -> n_is) += 1;
    }
    else if ((curr_state -> type) == 3) {
      (curr_info_node -> n_fs) += 1;
    }
    curr_state = next_state;
  }

  for (j=0; j<3; j++) {
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
parse_input_molcas (char * fn_infile) {

  int j,k,l,m,i,n,k_its; /* control loop variables */
  int mode; /* string matching mode flag */
  int match_start;
  int match_end;
  int last_int;
  int next_int;
  int idx_from;
  int idx_to;
  int extr_i;
  int n_states,n_trans;
  float extr_f;

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

  match_ln[0] = 0;
  match_ln[1] = (n_states+1)*256;
  match_ln[2] = (n_states)*256;
  match_ln[3] = match_ln[2] + ((n_trans+1)*256);

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

  num_idxs2[0] = 1;
  num_idxs2[1] = 2;
  num_idxs2[2] = 9;
  int n_idxs2 = 3;

  for (j=0; j<n_lookup_str; j+=2) {

    l = 0; /* index for string matches */
    m = 0;

    match_start = match_ln[j];
    match_end = match_ln[j+1];

    fseek(fp_tmpdata, match_start, 0);

    k_its = match_end-match_start;

    for (k=0; k<k_its; k++) {

      if ((c = fgetc(fp_tmpdata)) == EOF) {
        break;
      }

      str_buf[l] = (char)c;

      if ((str_buf[l] == '\n') && (l > 0)) { /* dont send blank lines */

        if ((j == 0) && (isempty(str_buf,l) != 1)) { /* extract energy eigenvalues and state indexes */

          get_numsl(str_buf,num_idxs1,l,n_idxs1,&e_eigval[m]);
          /* printf( "e_eigval[%d] = %le\n", m, e_eigval[m]); */
          /* fflush(stdout); */
          m++;
        }

        if ((j == 2) && (isempty(str_buf,l) != 1)) { /* extract transition moments and transition indexes */

          get_numsl(str_buf,num_idxs2,l,n_idxs2,&trans_idxs[0][m],\
                    &trans_idxs[1][m],&t_mom[m]);
          /* printf( "to %le from %le, %le \n", trans_idxs[0][m], trans_idxs[1][m], t_mom[m]); */
              /* fprintf(stderr, "\n\n=======Valgrind eject point=======\n\n"); */
              /* exit(1); */
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
  /* fflush(stdout); */
  /* for (k=0; k<n_trans; k++) { */
  /*   printf( "%le   %le   %le   %le   %le\n", parsed_input[0][k], parsed_input[1][k], parsed_input[2][k], parsed_input[3][k], parsed_input[4][k]); */
  /* } */

  fclose(fp_infile);
  fclose(fp_tmpdata);

  free(num_idxs1);
  free(num_idxs2);
  free(str_buf);
  free(e_eigval);
  free(trans_idxs);
  free(t_mom);

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
        printf( "found molcas.\n" );
        parse_input_molcas(fn_infile);
        printf( "processed molcas\n" );
        break; /* for */
      }

      else {
        printf( "format not found in input_formats.h\n" );
      }
    }
  }

  printf( "file processed\n" );
  return EXIT_SUCCESS;
}
