#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <signal.h>
#include <string.h>
#include "std_char_ops.h"
#include "std_num_ops.h"
#include "get_numsl.h"
#include "parse_input.h"
#include "input_formats.h"
#include "sci_const.h"
#include "state.h"
/* function

   * synopsis:

   * algorithm:

   * input:

   * output:

   * side-effects:

   */
int
ISINSIDE(double v,
         double r1,
         double r2) {

  if ((v>=r1) && (r2>=v)) {
    return 1;
  } else {
    return 0;
  }
}


#define BUF_SIZE 256

int nt;
int * idxs_map;
double ** parsed_input;
double tmax_d, tmax_q, e0;

int
check_pi (){
  int j = 0;

  /* check so that every state in PI can be reached with teh get_i function */
  while (1) {
    if ((int)parsed_input[0][j] != (int)parsed_input[0][get_i((int)parsed_input[0][j])]) {
      return 0;
    }

    j = get_inext((int)parsed_input[0][j]);

    if (parsed_input[0][j] == -1) {
      break;
    }
  }
  return 1;
}

void
set_tmax(){
  int j = 0;
  tmax_d = tmax_q = -1;

  while (1) {
    if (parsed_input[5][j] == 1) {
      if ( parsed_input[4][j] > tmax_d) {
        tmax_d = parsed_input[4][j];
      }
    }
    else{
      if (parsed_input[4][j] > tmax_q) {
        tmax_q = parsed_input[4][j];
      }
    }
    j++;
    if (parsed_input[0][j] == -1) {
      break;
    }
  }
}

void
state2s(int idx){
  int j = get_i(idx);

  fprintf(stderr, "\n\n================================\n");
  fprintf(stderr, "=======content of state %d=======\n",idx);
  while ((int)parsed_input[0][j] == idx) {
    printf( "\n%le   %le   %le   %le   %le  %le\n",
            parsed_input[0][j],
            parsed_input[1][j],
            parsed_input[2][j],
            parsed_input[3][j],
            parsed_input[4][j],
            parsed_input[5][j]);
    /* printf( "index data %d %d %d\n", idxs_map[last_i-1], idxs_map[last_i-1]+j,(int)parsed_input[0][idxs_map[last_i-1]+j]); */
    j++;
  }
  fprintf(stderr, "\n\n=======content of state %d=======\n",idx);
  fprintf(stderr,     "================================\n\n");

}

void
pi2s () {
  int j = 0;
  fprintf(stderr, "\n\n===========================\n");
  fprintf(stderr, "=======content of PI=======\n\n");
  while (1) {
    printf( "\n%le   %le   %le   %le   %le  %le\n",
            parsed_input[0][j],
            parsed_input[1][j],
            parsed_input[2][j],
            parsed_input[3][j],
            parsed_input[4][j],
            parsed_input[5][j]);

    /* printf( "index data %d %d %d\n", idxs_map[last_i-1], idxs_map[last_i-1]+j,(int)parsed_input[0][idxs_map[last_i-1]+j]); */
    j++;
    if (parsed_input[0][j] == -1) {
      break;
    }
  }
  fprintf(stderr, "\n\n=======content of PI=======\n");
  fprintf(stderr, "===========================\n\n");

}

double
get_efrom (int from) {

  int j = 0;
  while (1) {
    if ((int)parsed_input[0][j] == from) {
      return parsed_input[2][j];
    }
    j++;
    if ((int)parsed_input[0][j] == -1) {
      break;
    }
  }
  fprintf(stderr, "\n\nget_t:ERROR\n");
  printf( "program terminating due to the previous error.\n");
  fprintf(stderr, "\n\n=======Valgrind eject point=======\n\n");
  exit(1);
}

double
get_eto (int to) {

  int j = 0;
  while (1) {
    if ((int)parsed_input[1][j] == to) {
      return parsed_input[3][j];
    }
    j++;
    if (parsed_input[0][j] == -1) {
      break;
    }
  }

  fprintf(stderr, "\n\nget_t:ERROR\n");
  printf( "program terminating due to the previous error.\n");
  fprintf(stderr, "\n\n=======Valgrind eject point=======\n\n");
  exit(1);
}

double
get_t (int from,
       int to) {

  int j = 0;

  while (1) {
    if (((int)parsed_input[0][j] == from) &&
      ((int)parsed_input[1][j] == to)){
      return parsed_input[4][j];
    }
    j++;
    if ((int)parsed_input[0][j] == -1) {
      break;
    }
  }

  fprintf(stderr, "\n\nget_t:ERROR\n");
  printf( "program terminating due to the previous error.\n");
  fprintf(stderr, "\n\n=======Valgrind eject point=======\n\n");
  exit(1);
}

int
get_i (int from
       ) {
  int last_i = (int)parsed_input[0][0];
  int j = 0;
  while (last_i != -1) {

    if ((int)parsed_input[0][j] == from) {
      return j;
    }
    j++;
    last_i = (int)parsed_input[0][j];
  }
  /* fprintf(stderr, "\n\nget_i:ERROR, tried to get state %d\n", from); */
  /* printf( "program terminating due to the previous error.\n"); */
  /* fprintf(stderr, "\n\n=======Valgrind eject point=======\n\n"); */
  /* exit(1); */
  return -1;
}

int
get_inext (int from
           ) {
  int last_i = (int)parsed_input[0][0];
  int j      = 0;
  while (last_i != -1) {

    if ((int)parsed_input[0][j] == from) {
      /* if (parsed_input[0][j]== 4) { */
      /*   printf( "snipe2!\n" ); */
      /*   fprintf(stderr, "\n\n=======Valgrind eject point=======\n\n"); */
      /*   exit(1); */
      /* } */
      /* if (from == 4) { */
      /*   printf( "snipe1! %d\n",parsed_input[0][j] ); */
      /*   fprintf(stderr, "\n\n=======Valgrind eject point=======\n\n"); */
      /*   exit(1); */
      /* } */
      break;
    }
    j++;
    last_i = (int)parsed_input[0][j];
  }

  if ((int)parsed_input[0][j] != from) {
    fprintf(stderr, "\n\nget_inext:ERROR, cant get state %d\n", from);
    printf( "program terminating due to the previous error.\n");
    fprintf(stderr, "\n\n=======Valgrind eject point=======\n\n");
    exit(1);
  }


  while(1){
    if ((int)parsed_input[0][j] != from) {
      return j;
    }
    j++;
  }

  fprintf(stderr, "\n\nget_inext:ERROR\n");
  printf( "program terminating due to the previous error.\n");
  fprintf(stderr, "\n\n=======Valgrind eject point=======\n\n");
  exit(1);
}

/* mdda_s * */
/* screen_states (char * fn_relpath, */
/*                double * state_t, */
/*                double * state_er */
/*                ){ */
/*   printf( "  -screening states\n\n" ); */

/*   info_node inode = get_inode(fn_relpath); */

/*   int j,k; /\* looping variables *\/ */
/*   int n_states = inode -> n_states; */
/*   int gs_idx; /\* ground state index *\/ */
/*   int is_idx; /\* intermediate state index *\/ */

/*   int n_gfs = inode -> n_gfs; */
/*   int n_is = inode -> n_is; */

/*   /\* number of screen GS, IS, and FS *\/ */
/*   int n_sgs = 0; */
/*   int n_sis = 0; */
/*   int n_tmpsis = 0; /\* number of IS for a specific GS *\/ */
/*   int n_sfs = 0; */

/*   int * tmp_idxs; */

/*   /\* arrays used for storing the screened..  *\/ */
/*   mdda_s * igs; /\* .. ground states *\/ */
/*   mdda_s * iis; /\* .. intermediate states *\/ */

/*   /\* some hard-coded threshold values that will get replaced by the call */
/*    values *\/ */
/*   double bw_sum = inode -> bw_sum; */
/*   double mt_is, mt_fs; */
/*   double * tmp_trans; */

/*   estate root_state = inode -> root_e_state; */
/*   estate curr_state = root_state; */
/*   estate next_state; */
/*   estate gs_bm; /\* ground state bookmark *\/ */
/*   estate is_bm; /\* intermediate state bookmark *\/ */

/* /\* initialize the data structures used to store the screened states *\/ */
/*   igs = mdda_init(); */
/*   iis = mdda_init(); */

/*   /\* for a description of the screening algorithm below, look up this function */
/*      in the parse_input.h file *\/ */
/*   for (j=0; j<(inode -> mom_types[0]); j++) { */
/*       printf( "%le %le\n", (inode -> mt)[0][j],(inode -> mt)[1][j]); */
/*   } */

/*   n_sgs = 0; */
/*   while((next_state = curr_state -> next) != NULL){ */

/*     gs_idx = curr_state -> state_idx; */
/*     if ((curr_state -> type) == 1) { /\* if we found a ground state *\/ */
/*       /\* printf( "GS: %d %le\n", gs_idx,((curr_state ->  bw)/bw_sum)); *\/ */
/*       /\* perform stage 1 of the screening process *\/ */
/*       /\* if (((curr_state ->  bw)/(inode -> bw_sum)) >= state_t[0]) { *\/ */
/*       if ((((curr_state ->  bw)/bw_sum) >= state_t[1]) || (curr_state -> state_idx == 1)) { */
/*         /\* we found a ground state with high enough boltzmann weight *\/ */
/*         gs_idx = curr_state -> state_idx; */
/*         n_sgs++; */
/*         mdda_set(igs, 0, 0, n_sgs); */
/*         mdda_set(igs, 0, n_sgs, gs_idx ); */
/*         /\* printf( "found one GS at %d!\n\n", gs_idx ); *\/ */
/*         /\* printf( "\ngs[%d] = %d %le\n", n_sgs,gs_idx,((curr_state ->  bw)/bw_sum)); *\/ */
/*         /\* store a pointer to the final state matrix in the ground state root node *\/ */
/*         /\* if (gs_idx != 1) { *\/ */
/*         /\*   fprintf(stderr, "\n\n=======Valgrind eject point=======\n\n"); *\/ */
/*         /\*   exit(1); *\/ */
/*         /\* } *\/ */
/*         sleep(1); */
/*         /\* store a bookmark of where we found the last ground state *\/ */
/*         gs_bm = next_state; */

/*         /\* perform stage 2 of the screening process *\/ */
/*         n_tmpsis = 0; */
/*         tmp_trans = curr_state -> t_moms; */
/*         tmp_idxs = curr_state -> idxs_to; */
/*         for (j=0; j<(curr_state -> n_tfrom); j++) { */
/*           mt_is = (inode -> mt)[0][(curr_state -> ttypes[j])]; */
/*           printf( "  IS:j = %d/%d, sidx = %d, trans = %le, reltrans = %le, thrsh = %le\n, max = %le", j, (curr_state -> n_tfrom), tmp_idxs[j], tmp_trans[j], (tmp_trans[j]/mt_is), state_t[1]), mt_is; */
/*           /\* if IS transition has a transition moment above threshold *\/ */
/*           if ((tmp_trans[j]/mt_is) >= state_t[2]) { */
/*             is_bm = curr_state; */
/*             is_idx = tmp_idxs[j]; */
/*             printf( "  found one IS at %d!\n\n", is_idx ); */

/*            /\* locate the candidate intermediate state in the list and varify */
/*                that the state is indeed an IS (has type == 2) *\/ */

/*             /\* varify that it has been marked as an intermediate state *\/ */
/*             if ((curr_state = get_state_sil(inode, is_idx)) != NULL){ */
/*               /\* from the ground state found above, we also found an intermediate */
/*                  state with high enough relative transition moment *\/ */

/*               n_tmpsis++; */
/*               mdda_set(igs, n_sgs, 0, n_tmpsis); */
/*               mdda_set(igs, n_sgs, n_tmpsis, is_idx); */

/*               /\* printf( "  is[%d] = %d %le\n", n_sis, mdda_get(igs, n_sgs,n_sis), (curr_state -> t_moms)[j]/mt_is); *\/ */

/*               n_sfs = 0; */

/*               /\* before this is done, check so that the specific is_idx isnt */
/*                  in the iis already, in which case we would find its index in */
/*                  the 0th column of that mdda *\/ */
/*               if (mdda_intinint((iis->root), is_idx) == 0) { */
/*                 n_sis++; */
/*                 mdda_set(iis, 0, 0, n_sis); */
/*                 mdda_set(iis, 0, n_sis, is_idx); */

/*                 /\* perform stage 3 of the screening process *\/ */
/*                 /\* find final state transitions that pass stage 3 screening *\/ */
/*                 for (k=0; k < (curr_state -> n_tfrom); k++) { */
/*                   /\* printf( "    is %d -> FS:j = %d, sidx = %d, trans = %le, reltrans = %le, thrsh = %le\n",is_idx, k, ((curr_state -> idxs_to)[k]), ((curr_state -> t_moms)[k]), (((curr_state -> t_moms)[k])/mt_fs), state_t[3]); *\/ */

/*                   mt_fs = (inode -> mt[1][(curr_state -> ttypes[k])-1]); */
/*                   /\* screen the final state *\/ */
/*                   if ((((curr_state -> t_moms)[k])/mt_fs) >= state_t[3]) { */
/*                     /\* from the intermediate state found above, we also found a */
/*                        final state transition with high enough relative */
/*                        transition moment *\/ */
/*                     /\* printf( "    found one FS at %d!\n\n", mdda_get(iis, n_sis, n_sfs) ); *\/ */
/*                     n_sfs++; */
/*                     mdda_set(iis, n_sis, 0, n_sfs); */
/*                     mdda_set(iis, n_sis, n_sfs, (curr_state -> idxs_to)[k]); */

/*                     /\* printf("    is %d -> fs[%d][%d] = %d %le\n",is_idx, n_sis, n_sfs, mdda_get(iis, n_sis, n_sfs), (((curr_state -> t_moms)[k])/mt_fs)); *\/ */
/*                   }/\*  else { *\/ */
/*                   /\*   printf( "    is %d -> FS not accepted.\n",is_idx ); *\/ */
/*                   /\* } *\/ */
/*                 } */
/*               /\* if we looped over all final state transitions for this specific */
/*                  intermediate state and didnt find any FS transitions of high */
/*                  enough intensity, remove that intermediate state from the list *\/ */

/*                 /\* COMMENT OUT THIS BLOCK IF YOU WANT ELASTIC TRANSITIONS *\/ */
/*                 /\* if (n_sfs == 0) { *\/ */
/*                 /\*   /\\* printf( "\nis %d: %d %d %d\n",is_idx, n_sgs, n_tmpsis, n_sis); *\\/ *\/ */
/*                 /\*   /\\* mdda2s(igs); *\\/ *\/ */
/*                 /\*   /\\* mdda2s(iis); *\\/ *\/ */
/*                 /\*   mdda_set(iis, 0, 0, n_sis-1); *\/ */
/*                 /\*   mdda_set(iis, 0, n_sis, 0); *\/ */
/*                 /\*   mdda_set(iis, n_sis, 0, 0); *\/ */
/*                 /\*   mdda_set(iis, n_sis, n_sfs, 0); *\/ */

/*                 /\*   mdda_set(igs, n_sgs, 0, n_tmpsis-1); *\/ */
/*                 /\*   mdda_set(igs, n_sgs, n_tmpsis, 0); *\/ */
/*                 /\*   /\\* mdda_set(igs, n_sgs, n_tmpsis, n_tmpsis); *\\/ *\/ */
/*                 /\*   /\\* mdda_show(igs); *\\/ *\/ */

/*                 /\*   n_sis--; *\/ */
/*                 /\*   n_tmpsis--; *\/ */
/*                 /\*   /\\* mdda2s(igs); *\\/ *\/ */
/*                 /\*   /\\* mdda2s(iis); *\\/ *\/ */
/*                 /\*   /\\* printf( "\nstate treshold: %le\n", (inode -> mt_is)*state_t[2]); *\\/ *\/ */
/*                 /\*   /\\* e_state2s(get_state_si(inode, gs_idx),2); *\\/ *\/ */

/*                 /\*   /\\* printf( "\nstate treshold: %le\n", (inode -> mt_fs)*state_t[3]); *\\/ *\/ */
/*                 /\*   /\\* e_state2s(get_state_si(inode, is_idx),1); *\\/ *\/ */
/*                 /\*   /\\* fprintf(stderr, "\n\n=======Valgrind eject point=======\n\n"); *\\/ *\/ */
/*                 /\*   /\\* exit(1); *\\/ *\/ */
/*                 /\* } *\/ */
/*               } */
/*             } */

/*             curr_state = is_bm; /\* reset the state to the last initial state *\/ */
/*           } else { */
/*             printf( "  IS not accepted.\n" ); */
/*           } */
/*         } */
/*         /\* jump back to the last ground state that was found in the list *\/ */
/*         next_state = gs_bm; */
/*       } /\* else { *\/ */
/*       /\* printf( "unacceptable bw: gs[%d] = %d %le\n", n_sgs,gs_idx,((curr_state ->  bw))); *\/ */
/*       /\* } *\/ */
/*     } */

/*     curr_state = next_state; */
/*   } */

/*   /\* check so that every IS index in the igs is also in the iis array *\/ */
/*   if (isiniis(igs,iis)) { */
/*     fprintf(stderr, "Error: parse_input.c, function screen_states, while \ */
/* screening the transitions, the arrays storing the state indices got corrupted \ */
/*  and some intermediate states will be excluded from the RIXS map. this is not \ */
/*  acceptable. \n"); */
/*     printf( "program terminating due to the previous error.\n"); */
/*     exit(EXIT_FAILURE); */
/*   } */
/*   /\* mdda_show(igs); *\/ */

/*   mdda_set_branch(igs, iis, 0, 0 ); */
/*   mdda_show(igs); */

/*   /\* fprintf(stderr, "\n\n=======Valgrind eject point=======\n\n"); *\/ */
/*   /\* exit(1); *\/ */
/*   return (igs->root); */
/* } */

/* int */
/* init_data_branch(double * state_er, /\* state energy ranges *\/ */
/*                  double ** pi, /\* parsed input *\/ */
/*                  int * mom, */
/*                  int ns, /\* n states *\/ */
/*                  int nt, /\* n transitions *\/ */
/*                  char * fs /\* input file name string *\/ */
/*                  ){ */
/*   int j; /\* looping variables *\/ */

/*   /\* initialize the linked list structure for the electronic states *\/ */
/*   init_estate_list(fs, mom, ns, nt); */

/*   /\* use the extracted input data to set the state of the linked list *\/ */
/*   set_estate_list(state_er, pi, mom, ns, nt, fs); */

/*   /\* free up parsed_input, since it now lives in the tree *\/ */
/*   for (j=0; j<5; j++) { */
/*     free(pi[j]); */
/*   } */

/*   free(pi); */

/*   return EXIT_SUCCESS; */
/* } */

int
parse_input_tmp (double * state_er,
                 char * fn_tmpdata
                 ) {

  int j,k,l,m,i,n,k_its; /* control loop variables */
  int mode; /* string matching mode flag */
  int match_start;
  int match_end;
  int idx_from;
  int idx_to;
  int tmp_idx2;
  int n_states;
  int n_trans;
  int last_i;

  int n_proc;

  int delim_count,trs_type;

  int * num_idxs1;
  int * num_idxs2;
  int * trs_types;
  int * proc_idx;

  double tmp_idx = 0;
  double maxr    = get_maxl(state_er,state_er[0]);

  double *  e_eigval;
  int *     idxs_eigval;
  double ** pi_buf;
  double ** trans_idxs;
  double *  t_mom;

  FILE * fp_tmpdata;

  int n_matchfound  = 0;
  int match_vals[2] = {0,0}; /* place to store the indexes of the lines containing the
                              matches */

  int c; /* temporary char for storing input file characters */
  /* we are looking for four strings in the output file: one at the beginning
     of each data block, and one at the end. */
  char * str_buf = malloc(BUF_SIZE*2);
  const char DAT_DELIM[32] =  "============ data block finish\n";

  /* open the input file */
  if((fp_tmpdata = fopen(fn_tmpdata, "r")) == NULL) {
    fprintf(stderr,"parse_input: unable to open the input file %s.\n",fn_tmpdata);
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  k = 0; /* index for tmp_string */
  l = 0; /* index for lookup string */

  /* storage for the energy eigenvalues */
  if((e_eigval = malloc(5000000*sizeof(double))) == NULL ){
    fprintf(stderr, "parse_input_molcas, malloc: failed to allocate memory for\
 \"e_eigval\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }
  if((trs_types = malloc(5000000*sizeof(double))) == NULL ){
    fprintf(stderr, "parse_input_molcas, malloc: failed to allocate memory for\
 \"e_eigval\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }
  if((idxs_eigval = malloc(5000000*sizeof(int))) == NULL ){
    fprintf(stderr, "parse_input_molcas, malloc: failed to allocate memory for\
 \"idxs_eigval\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }
  for (j=0; j<5000000; j++) {
    idxs_eigval[j] = -1;
  }
  /* storage for the transition moments */
  if((t_mom = malloc(5000000*sizeof(double))) == NULL ){
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
    if((trans_idxs[j] = malloc(5000000*sizeof(double))) == NULL ){
      fprintf(stderr, "parse_input_molcas, malloc: failed to allocate memory for\
 \"e_eigval\"\n");
      printf( "program terminating due to the previous error.\n");
      exit(1);
    }
  }

  if((num_idxs1 = malloc(2*sizeof(int))) == NULL ){
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

  num_idxs1[0] = 0;
  num_idxs1[1] = 1;
  int n_idxs1  = 2;

  num_idxs2[0] = 0;
  num_idxs2[1] = 1;
  num_idxs2[2] = 2;
  int n_idxs2  = 3;

  j = l = m = 0;

  n_states    = 0;
  n_trans     = 0;
  delim_count = 0;
  trs_type    = 1;

  /* now that data structures of the rights size has memory allocated for
     them, start reading data from the temporary file */
  /* rewind(fp_tmpdata); */
  while ((c = fgetc(fp_tmpdata)) != EOF) {

    str_buf[l] = (char)c;

    if ((str_buf[l]                  == '\n') && (l > 0)) { /* dont send blank lines */
      if (strcmp(DAT_DELIM, str_buf) <= 0) {
        delim_count++;
        l = 0; /* index for string matches */
        j=2;
        /* printf( "Found the delimiter\n" ); */
        /* sleep(1); */
        if (delim_count == 2) {
          /* read quadrupole transition */
          trs_type                    = 2;
        }
      }
      else{
        if ((j == 0) && (isempty(str_buf,l) != 1)) {

          /* extract energy eigenvalues and state indexes */

          get_numsl(str_buf,num_idxs1,l,n_idxs1,&tmp_idx,&e_eigval[n_states]);

          /* printf( "pre%d\n", j); */
          if ((fabs(e_eigval[n_states]-e_eigval[0])*AUTOEV) \
              < (maxr + 100 )) {
            /* printf( "e_eigval[%d]  = %le, %d\n", n_states+1, e_eigval[n_states],tmp_idx ); */
            idxs_eigval[(int)tmp_idx] = n_states;
            n_states++;

          }
          /* printf( "post%d\n", j); */
        }
        else if ((j == 2) && (isempty(str_buf,l) != 1)) {

          /* extract transition moments and transition indexes */
          get_numsl(str_buf,num_idxs2,l,n_idxs2,&trans_idxs[0][n_trans],\
                    &trans_idxs[1][n_trans],&t_mom[n_trans]);
          trs_types[n_trans] = trs_type;

          if (!((idxs_eigval[(int)(trans_idxs[0][n_trans])] == -1) || (idxs_eigval[(int)(trans_idxs[1][n_trans])] == -1)) &&
              (((fabs(e_eigval[idxs_eigval[(int)(trans_idxs[1][n_trans])]]-e_eigval[0])*AUTOEV) \
                < maxr) &&                \
               ((fabs(e_eigval[idxs_eigval[(int)(trans_idxs[0][n_trans])]]-e_eigval[0])*AUTOEV) \
                < maxr)))
            {
              /* printf( "%le %le %le\n", trans_idxs[0][n_trans], trans_idxs[1][n_trans], t_mom[n_trans]); */
              n_trans++;
          }
        }
        l=0; /* reset the buffer write head to start reading a the next line */
      }
    }
    else{
      l++;
    }
  }
  /* fprintf(stderr, "\n\n=======Valgrind eject point=======\n\n"); */
  /* exit(1); */
  /* allocate space for the "parsed input matrix" that will be filled with data
     in the remaining sections of this function */

  if((proc_idx = malloc(n_states*sizeof(int))) == NULL ){
    fprintf(stderr, "parse_input_molcas, malloc: failed to allocate memory for\
 \"proc_idx\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  if((idxs_map = malloc(n_states*sizeof(int))) == NULL ){
    fprintf(stderr, "parse_input_molcas, malloc: failed to allocate memory for\
 \"idxs_map\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  if((parsed_input = malloc(6*sizeof(double*))) == NULL ){
    fprintf(stderr, "parse_input_molcas, malloc: failed to allocate memory for\
 \"input_data\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  for (j=0; j<6; j++) {
    if((parsed_input[j] = malloc((n_trans+1)*sizeof(double))) == NULL ){
      fprintf(stderr, "parse_input_molcas, malloc: failed to allocate memory \
for pointers in \"input_data\"\n");
      printf( "program terminating due to the previous error.\n");
      exit(1);
    }
  }
  parsed_input[0][n_trans] = -1; /* end of list */

  if((pi_buf = malloc(6*sizeof(double*))) == NULL ){
    fprintf(stderr, "parse_input_molcas, malloc: failed to allocate memory for\
 \"input_data\"\n");
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  for (j=0; j<6; j++) {
    if((pi_buf[j] = malloc(n_trans*sizeof(double))) == NULL ){
      fprintf(stderr, "parse_input_molcas, malloc: failed to allocate memory \
for pointers in \"input_data\"\n");
      printf( "program terminating due to the previous error.\n");
      exit(1);
    }
  }

  /* finally, store the data in the parsed_input matrix for the parse_input function  */
  n_proc      = 0;
  idxs_map[0] = 1;
  last_i = 1;
  j           = k = l = 0;
  e0 = e_eigval[0];

  while (l<n_trans-1) {
    if ((j+l) == (n_trans)) {
      /* dont read beyond the last value in PI */
      /* for (k=0; k<j; k++) { */
      /*   printf( "%le %le %le %le %le\n", pi_buf[0][k],pi_buf[1][k],pi_buf[2][k],pi_buf[3][k],pi_buf[4][k],pi_buf[5][k]); */
      /* } */
      /* break; */
      pi_buf[0][j] = -1;
      idx_from = -1;
    } else {

      /* if ((trs_types[j+l] == 1) && */
      /*     ISINSIDE(e_eigval[idxs_eigval[(int)trans_idxs[1][j+l]]],state_er[1],state_er[2]) && */
      /*     ISINSIDE(e_eigval[idxs_eigval[(int)trans_idxs[0][j+l]]],state_er[3],state_er[4])) */
      /*   { */
      /* if (((trans_idxs[1][j+l] == 572) || */
      /*     (trans_idxs[1][j+l]      == 573)) && */
      /*     (trs_types[j+l] == 1) */
      /*     ) { */
      if ((ISINSIDE((e_eigval[idxs_eigval[(int)trans_idxs[1][j+l]]]-e0)*AUTOEV,state_er[3],state_er[4]) &&
           (ISINSIDE((e_eigval[idxs_eigval[(int)trans_idxs[0][j+l]]]-e0)*AUTOEV,state_er[5],state_er[6]))) &&
          (trs_types[j+l]             == 1)
          ) {
        idx_from = trans_idxs[1][j+l];
        idx_to = trans_idxs[0][j+l];
        /* printf( "%d %d %d\n", (int)trans_idxs[0][j+l], (int)trans_idxs[1][j+l], (int)trs_types[j+l]); */
        /* fprintf(stderr, "\n\n ======= Valgrind eject point=======\n\n"); */
        /* exit(1); */
      }
      else {
        idx_from = trans_idxs[0][j+l];
        idx_to   = trans_idxs[1][j+l];
      }
      pi_buf[0][j] = idx_from;
      pi_buf[1][j] = idx_to;
      pi_buf[2][j] = e_eigval[idxs_eigval[idx_from]];
      pi_buf[3][j] = e_eigval[idxs_eigval[idx_to]];
      pi_buf[4][j] = t_mom[j+l];
      pi_buf[5][j] = trs_types[j+l];
    }
    if (idx_from != last_i) {

      /* printf( "\n\n === idxs_from = %d, last_i    = %d === \n\n", idx_from, last_i); */
      /* we have read all transitions for a state */
      /* check if the last_i has already been processed */
      if (intinint(proc_idx,last_i,n_proc)          == 0) {
        /* printf( "\n\n                          ==== adding data for state %d====\n\n", last_i ); */
        m                                            = l;
        while ((int)pi_buf[0][m-l] != idx_from) {
          /* if (idx_from == -1) { */
          /*   printf( "PI_BUF,from                  = %d, parsed = %d,m = %d\n\n",(int)pi_buf[0][m],(int)pi_buf[1][m],m); */
          /* } */
          parsed_input[0][m]                         = pi_buf[0][m-l];
          parsed_input[1][m]                         = pi_buf[1][m-l];
          parsed_input[2][m]                         = pi_buf[2][m-l];
          parsed_input[3][m]                         = pi_buf[3][m-l];
          parsed_input[4][m]                         = pi_buf[4][m-l];
          parsed_input[5][m] = pi_buf[5][m-l];
          m++;
        }

        /* printf( "\n\n====added data for state %d====\n\n",last_i ); */
        proc_idx[n_proc++] = last_i;

      }
      else{
        /* printf( "\n\n====splicing data for state %d ==== \n\n", last_i ); */
        /* else, use fwdsplice to splice the buffer into PI */

        /* fwdsplice(pi_buf,parsed_input,l,l+j,6,j); */
        /* printf( "\n\n===============splice state %d ..",last_i ); */

        /* fflush(stdout); */
        /* printf( "\n\n\n\ngetting inext\n\n\n\n" ); */
        tmp_idx2 = get_inext(last_i);
        /* printf( "\n\n\n\ngot inext\n\n\n\n" ); */
        /* if (last_i == 32) { */
        /*   printf( "\n\nfound four\n\n\n" ); */
        /*   tmp_idx2 = get_i(32); */
        /*   printf( "\n\n\n\ngot istate\n\n\n\n" ); */
        /* state2s(3); */
        /* state2s(4); */
        /*   printf( "\n\ngot state\n\n\n" ); */
        /*   fprintf(stderr, "\n\n ======= Valgrind eject point=======\n\n"); */
        /*   exit(1); */
        /* } */

        fwdsplice(pi_buf,parsed_input,tmp_idx2,l,j,6);
        /* state2s(3); */
        /* state2s(4); */
        /* fprintf(stderr, "\n\n ======= Valgrind eject point=======\n\n"); */
        /* exit(1); */
        /* printf( "..done  ===============\n\n" ); */
        /* fflush(stdout); */
        /* state2s(4); */
        /* sleep(2); */
        /* fprintf(stderr, "\n\n=======Valgrind eject point=======\n\n"); */
        /* exit(1); */
        /* printf( "====\n\nspliced data for state %d====\n\n",last_i ); */
      }
      /* if (last_i == 5) { */
      /*   fprintf(stderr, "\n\n=======Valgrind eject point=======\n\n"); */
      /*   exit(1); */
      /* } */

      l += j;
      j = 0;
      /* idxs_map[(int)(idx_from)-1] = l; */
      /* printf( "%d %d\n", (int)(idx_from)-1, l-1); */
      /* if (idx_from == 7) { */
      /*   pi2s(); */
      /*   /\*         printf( "%le\n",         get_t(3,572)); *\/ */
      /*   /\* printf( "%le\n",         get_t(3,573)); *\/ */
      /*   /\* printf( "%le\n",         get_t(1,3)); *\/ */
      /*   /\* printf( "%le\n",         get_t(1,5)); *\/ */
      /*   /\* printf( "%le\n",         get_t(1,15)); *\/ */
      /*   /\* printf( "%le\n",         get_t(1,16)); *\/ */
      /*   fprintf(stderr, "\n\n=======Valgrind eject point=======\n\n"); */
      /*   exit(1); */

      /* } */
    } else {
      j++;
    }
    last_i = idx_from;

  }

  parsed_input[0][l] = -1;

  nt = l;

  if (check_pi() == 0) {
    fprintf(stderr, "\n\nparse_input.c, function parse_input: input matrix integrity check failure.\n");
    printf( "program terminating due to the previous error.\n");
    exit(EXIT_FAILURE);
  }

  /* state2s(2); */
  /* state2s(3); */
  /* state2s(4); */
  /* state2s(5); */
  /* pi2s(); */
  /* printf( "%d %d %d\n", get_i(1), get_i(2), get_i(3)); */

  set_tmax();

  fclose(fp_tmpdata);

  free(str_buf);
  free(num_idxs1);
  free(num_idxs2);

  free(e_eigval);

  free(t_mom);

  free(idxs_map);
  free(proc_idx);

  for (j=0; j<6; j++) {
    free(pi_buf[j]);
  }
  free(pi_buf);

  for (j=0; j<2; j++) {
    free(trans_idxs[j]);
  }
  free(trans_idxs);

  return EXIT_SUCCESS;

}

int
parse_input (double * state_er,
             char * fn_relpath, /* name of input file */
             int len_fn){

  int j, k, l; /* looping variables */
  int n_states,n_trans;

  char format[BUF_SIZE] = {0};

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
  j = len_fn;
  while(fn_relpath[--j] != '.'){};

  /* loop over the input file name and extract the ending */
  for (k=0; j<len_fn; j++) {
    format[k] = fn_relpath[j];
    k++;
  }
  format[k] = '\0';

  parse_input_tmp(state_er, fn_relpath);

  return EXIT_SUCCESS;
}
