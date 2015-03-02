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

#define BUF_SIZE 256

int n_info_nodes; /* counter for number of  */
double * input_data[4];
int ** state_indices;

info_node info_node_root;



int
sort_states (int n_args,
             ...){
  int j;

  /* some hard-coded threshold values that will get replaced by the call
   values */
  float thrsh_vals[3] = {0.5, 0.0001, 0.0001};
  float tmp_thrsh;

  va_list argv;
  va_start(argv, n_args);


  for (j=0; j<n_args; j++) {
    tmp_thrsh = (float)va_arg(argv, double); /* grab the next vararg */
    thrsh_vals[j] = tmp_thrsh;
  }

  /* first, we screen for the initial states. */

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

  /* group_states(fs); */

  /* free up parsed_input, since it now lives in the tree */
  for (j=0; j<5; j++) {
    free(pi[j]);
  }

  free(pi);
  return 0;
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

  return 0;
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
  info_node curr_info_node = info_node_root;
  info_node next_info_node;

  e_state end_state; /* pointer to the last state in the llist */
  e_state next_state;
  e_state curr_state;

  /* outer loop state pointers */
  e_state next_state2;
  e_state curr_state2;

  int * tmp_idxs;
  double * tmp_tmoms;
  double * tmp_evals;
  double * e_vals;

  if((tmp_idxs = malloc(n_states*sizeof(int))) == NULL ){
    fprintf(stderr, "parse_input:function init_state_ll, malloc: failed \
to allocate memory for \"idxs_to\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  if((tmp_tmoms = malloc(n_states*sizeof(double))) == NULL ){
    fprintf(stderr, "parse_input:function init_state_ll, malloc: failed \
to allocate memory for \"idxs_to\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  if((tmp_evals = malloc(n_states*sizeof(double))) == NULL ){
    fprintf(stderr, "parse_input:function init_state_ll, malloc: failed \
to allocate memory for \"idxs_to\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  if((e_vals = malloc(n_states*sizeof(double))) == NULL ){
    fprintf(stderr, "parse_input:function init_state_ll, malloc: failed \
to allocate memory for \"idxs_to\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  if((groups = malloc(3*sizeof(int*))) == NULL ){
    fprintf(stderr, "parse_input:function init_state_ll, malloc: failed \
to allocate memory for \"idxs_to\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  for (j=0; j<3; j++) {
    if((groups[j] = malloc((n_states+1)*sizeof(int))) == NULL ){
      fprintf(stderr, "parse_input:function init_state_ll, malloc: failed \
to allocate memory for \"idxs_to\"\n");
      printf( "program terminating due to the previous error.\n");
      exit(1);
    }
  }

  /* loop over the info node list and find the one with the right string
   identifier */

  while(strstr((curr_info_node -> str_id),id) == NULL) {

    next_info_node = curr_info_node -> next;
    curr_info_node = next_info_node;

    if (curr_info_node == NULL) {
      fprintf(stderr, "parse_input.c: set_state_ll, no info node can be found\
 with a str_id == %s\n", id);
      printf( "program terminating due to the previous error.\n");
      exit(1);
    }
  }

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
      next_state = curr_state -> next;
      curr_state = next_state;

      if (((curr_state -> next) == NULL) && ((j+1) <= n_trans)) {
        fprintf(stderr, "parse_input.c:set_state_node, the list of electronic states ended before all input data was transfered.\n");
        printf( "program terminating due to the previous error.\n");
        exit(1);
      }

      from_last = from;
      k=0;
      j--;
    }
    else {
      tmp_idxs[k] = parsed_input[1][j];
      tmp_evals[k] = parsed_input[3][j];
      tmp_tmoms[k] = parsed_input[4][j];
      k++;  /* increase counter for transitions counted */
    }
  }

  /* the last node in the llist is the highest energy intermediate state */
  curr_state -> type = 3;
  curr_info_node -> bw_sum = tmp_sum;

  /* loop over the electronic state llist and set the rel_bw property
     of each state, also group the ground states */
  /* curr_state = (curr_info_node -> root_e_state); */

  /* curr_state -> type = 1; */

  /* first look up all states that do not transition to the lowest
   energy ground state. these are ground states. */
  /* int * gs_idxs = malloc((curr_info_node -> n_states)*sizeof(int)); */

  /* for (j=0; j<(curr_info_node -> n_states); j++) { */
  /*   gs_idxs[j] = j+1; */
  /*   printf( "%d", gs_idxs[j]); */
  /* } */
  /* printf( "\n" ); */
  /* int * idxs_res = getintinint(((curr_info_node -> root_e_state)->idxs_to),gs_idxs,(curr_info_node -> n_states)); */

  /* for (j=0; j<(curr_info_node -> root_e_state)->n_tfrom; j++) { */
  /*   printf( "%d", idxs_res[j]); */
  /* } */
  /* printf( "\n" ); */

  /* curr_state = curr_info_node -> root_e_state; */
  /* for (j=0; j<=n_states; j++){ */
  /*   if (intinint(idxs_res,(curr_state->state_idx),(curr_info_node -> root_e_state)->n_tfrom) == 0) { */
  /*     curr_state -> type = 1; */
  /*     printf( "state %d is %d\n",curr_state -> state_idx, 1 ); */
  /*   } */

  /*   next_state = curr_state -> next; */
  /*   curr_state = next_state; */
  /*   if ((curr_state -> next)->next == NULL){ */
  /*     break; */
  /*   } */
  /* } */

  /* free(idxs_res); */
  /* free(gs_idxs); */

  k_meansl(e_vals, groups, l);

  end_state = curr_state -> last;
  curr_state2 = end_state;

  j = 0;
  for (k=1; k<=groups[j][0]; k++) {

    sorted_idx = groups[j][k];
    g = &groups[j][1];
    curr_state = curr_info_node -> root_e_state;

    while((next_state = curr_state -> next) != NULL){

      curr_state -> rel_bw = (curr_state -> bw)/tmp_sum;
      to = (curr_state -> state_idx);

      if (to == sorted_idx) {
        /* store a pointer to the "to" state idxs_to array and check so that every foundstate is not inside of it as well */
        ito = (curr_state -> idxs_to);
        tmp_ntfrom = curr_state -> n_tfrom;
        tmp_type = curr_state -> type;
        printf( "\nto = %d, type = %d\n", to, tmp_type);
        /* sleep(1); */

        /* loop over the llist once more and find all states that do not have
           any transitions to the selected state */
        curr_state = curr_info_node -> root_e_state;
        while((next_state = curr_state -> next) != NULL){
          tmp_idx = curr_state -> state_idx;
          printf( "from %d to %d, type %d\n", to, tmp_idx, curr_state -> type);
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
                } else if (tmp_type == 1) {
                  curr_state -> type = 1;
                } else if ((tmp_type == 2) && (curr_state -> type == 0)) {
                  curr_state -> type = 2;
                }
                printf( "foundstate %d\n", curr_state -> state_idx);
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

  while((next_state = curr_state -> next) != NULL){
   printf( "type = %d\n", curr_state -> type);
    curr_state = next_state;
  }

  for (j=0; j<3; j++) {
    free(groups[j]);
  }
  free(groups);
  fprintf(stderr, "\n\n=======Valgrind eject point=======\n\n");
  exit(1);
  /* trim off any nodes at the end of the list that hasnt gotten allocated. */
  free(tmp_idxs);
  free(tmp_tmoms);
  free(tmp_evals);
  free(e_vals);

  return 0;
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

  /* construct a linked-list tree of the parsed input */
  /* load the parsed input into the linked list of electronic states */
  /* set_state_ll(parsed_input,fn_infile); */

  /* the input data is processed into the info_node linked list struct
   so lets free up space for that */
  /* for (j=0; j<5; j++) { */
  /*   free(parsed_input[j]); */
  /* } */
  /* free(parsed_input); */
  /* sort_states(parsed_input); */

  printf( "file processed\n" );
  return 0;
}
