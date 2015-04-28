#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <signal.h>
#include <string.h>
#include "std_char_ops.h"
#include "std_num_ops.h"
#include "appc.h"
#include "get_numsl.h"
#include "parse_input.h"
#include "input_formats.h"
#include "sci_const.h"
#include "dynarray.h"
#include "structs.h"
#include "info_node.h"
#include "estate.h"

#define BUF_SIZE 256

info_node root_inode;

mdda_s *
screen_states (char * fn_infile,
               double * state_t,
               double * state_er
               ){
  printf( "  -screening states\n\n" );

  info_node inode = get_inode(fn_infile);

  int j,k; /* looping variables */
  int n_states = inode -> n_states;
  int gs_idx; /* ground state index */
  int is_idx; /* intermediate state index */

  int n_gfs = inode -> n_gfs;
  int n_is = inode -> n_is;

  /* number of screen GS, IS, and FS */
  int n_sgs = 0;
  int n_sis = 0;
  int n_tmpsis = 0; /* number of IS for a specific GS */
  int n_sfs = 0;
  int n_tmpsfs; /* number of FS for a specific IS */

  int * tmp_idxs;

  /* arrays used for storing the screened..  */
  mdda_s * igs; /* .. ground states */
  mdda_s * iis; /* .. intermediate states */

  /* some hard-coded threshold values that will get replaced by the call
   values */
  double bw_sum = inode -> bw_sum;

  double * tmp_trans;

  estate root_state = inode -> root_e_state;
  estate curr_state = root_state;
  estate next_state;
  estate gs_bm; /* ground state bookmark */
  estate is_bm; /* intermediate state bookmark */

/* initialize the data structures used to store the screened states */
  igs = mdda_init();
  iis = mdda_init();

  /* for a description of the screening algorithm below, look up this function
     in the parse_input.h file */

  n_sgs = 0;
  while((next_state = curr_state -> next) != NULL){

    gs_idx = curr_state -> state_idx;
    if ((curr_state -> type) == 1) { /* if we found a ground state */
      /* printf( "GS: %d %le\n", gs_idx,((curr_state ->  bw)/bw_sum)); */
      /* perform stage 1 of the screening process */
      /* if (((curr_state ->  bw)/(inode -> bw_sum)) >= state_t[0]) { */
      if ((((curr_state ->  bw)/bw_sum) >= state_t[1]) || (curr_state -> state_idx == 1)) {
        /* we found a ground state with high enough boltzmann weight */
        gs_idx = curr_state -> state_idx;
        n_sgs++;
        mdda_set(igs, 0, 0, n_sgs);
        mdda_set(igs, 0, n_sgs, gs_idx );
        /* printf( "found one GS at %d!\n\n", gs_idx ); */
        /* printf( "\ngs[%d] = %d %le\n", n_sgs,gs_idx,((curr_state ->  bw)/bw_sum)); */
        /* store a pointer to the final state matrix in the ground state root node */

        /* sleep(1); */
        /* store a bookmark of where we found the last ground state */
        gs_bm = next_state;

        /* perform stage 2 of the screening process */
        n_tmpsis = 0;
        tmp_trans = curr_state -> t_moms;
        tmp_idxs = curr_state -> idxs_to;
        for (j=0; j<(curr_state -> n_tfrom); j++) {
          /* printf( "  IS:j = %d/%d, sidx = %d, trans = %le, reltrans = %le, thrsh = %le\n", j, (curr_state -> n_tfrom), tmp_idxs[j], tmp_trans[j], (tmp_trans[j]/(inode -> mt_is)), state_t[1]); */
          /* if IS transition has a transition moment above threshold */
          if ((tmp_trans[j]/(inode -> mt_is)) >= state_t[2]) {
            is_bm = curr_state;
            is_idx = tmp_idxs[j];
            /* printf( "  found one IS at %d!\n\n", is_idx ); */

           /* locate the candidate intermediate state in the list and varify
               that the state is indeed an IS (has type == 2) */

            /* varify that it has been marked as an intermediate state */
            if (((curr_state = get_state_si(inode, is_idx)) != NULL)\
                && ((curr_state -> type) == 2)){
              /* from the ground state found above, we also found an intermediate
                 state with high enough relative transition moment */

              n_tmpsis++;
              mdda_set(igs, n_sgs, 0, n_tmpsis);
              mdda_set(igs, n_sgs, n_tmpsis, is_idx);

              /* printf( "  is[%d] = %d %le\n", n_sis, mdda_get(igs, n_sgs,n_sis), (curr_state -> t_moms)[j]/(inode -> mt_is)); */

              n_sfs = 0;

              /* before this is done, check so that the specific is_idx isnt
                 in the iis already, in which case we would find its index in
                 the 0th column of that mdda */
              if (mdda_intinint((iis->root), is_idx) == 0) {
                n_sis++;
                mdda_set(iis, 0, 0, n_sis);
                mdda_set(iis, 0, n_sis, is_idx);

                /* perform stage 3 of the screening process */
                /* find final state transitions that pass stage 3 screening */
                for (k=0; k < (curr_state -> n_tfrom); k++) {
                  /* printf( "    FS:j = %d, sidx = %d, trans = %le, reltrans = %le, thrsh = %le\n", k, ((curr_state -> idxs_to)[k]), ((curr_state -> t_moms)[k]), (((curr_state -> t_moms)[k])/(inode -> mt_fs)), state_t[2]); */
                  if ((((curr_state -> t_moms)[k])/(inode -> mt_fs)) >= state_t[3]) {
                    /* from the intermediate state found above, we also found a
                       final state transition with high enough relative
                       transition moment */
                    /* printf( "    found one FS!\n\n" ); */
                    n_sfs++;
                    mdda_set(iis, n_sis, 0, n_sfs);
                    mdda_set(iis, n_sis, n_sfs, (curr_state -> idxs_to)[k]);

                    /* printf("    fs[%d][%d] = %d %le\n", n_sis, n_sfs, mdda_get(iis, n_sis, n_sfs), (((curr_state -> t_moms)[k])/(inode -> mt_fs))); */
                  }/*  else { */
                  /* printf( "    FS not accepted.\n" ); */
                  /* } */
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
          }/*  else { */
          /*   printf( "  IS not accepted.\n" ); */
          /* } */
        }
        /* jump back to the last ground state that was found in the list */
        next_state = gs_bm;
      } /* else { */
      /* printf( "unacceptable bw: gs[%d] = %d %le\n", n_sgs,gs_idx,((curr_state ->  bw))); */
      /* } */
    }

    curr_state = next_state;
  }
  /* mdda_show(igs); */
  /* fprintf(stderr, "\n\n=======Valgrind eject point=======\n\n"); */
  /* exit(1); */
  mdda_set_branch(igs, iis, 0, 0 );

  return (igs->root);
}

int
init_data_branch(double * state_er, /* state energy ranges */
                 double ** pi, /* parsed input */
                 int ns, /* n states */
                 int nt, /* n transitions */
                 char * fs /* input file name string */
                 ){
  int j; /* looping variables */

  /* initialize the linked list structure for the electronic states */
  init_estate_list(fs, ns, nt);

  /* use the extracted input data to set the state of the linked list */
  set_estate_list(state_er, pi, ns, nt, fs);

  /* free up parsed_input, since it now lives in the tree */
  for (j=0; j<5; j++) {
    free(pi[j]);
  }

  free(pi);

  return EXIT_SUCCESS;
}

int
parse_input_molcas (double * state_er,
                    char * fn_infile
                    ) {

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

  int c; /* temporary char for storing input file characters */
  /* we are looking for four strings in the output file: one at the beginning
     of each data block, and one at the end. */
  char * str_buf = malloc(BUF_SIZE*2);
  const char DAT_DELIM[32] =  "============ data block finish\n";

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

  if((fp_tmpdata = fopen("/home/kimchi/dev/smap/tmp/molcas_data.tmp", "w+")) == NULL) {
    fprintf(stderr,"parse_input: unable to open the output file %s.\n",fn_infile);
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  k = 0; /* index for tmp_string */
  l = 0; /* index for lookup string */
  m = 0; /* index for string matches */
  mode = 0; /* start of in string search mode */

  /* read the Molcas input file */
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
            m++;
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
              match_vals[n_matchfound];

              /* write the extracted input data to a temporary file */
              fprintf(fp_tmpdata,DAT_DELIM, n_matchfound );
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

  num_idxs1[0] = 1;
  int n_idxs1 = 1;

  num_idxs2[0] = 0;
  num_idxs2[1] = 1;
  num_idxs2[2] = 6;
  int n_idxs2 = 3;
  j = l = m = 0;

  /* now that data structures of the rights size has memory allocated for
     them, start reading data from the temporary file */
  rewind(fp_tmpdata);
  while ((c = fgetc(fp_tmpdata)) != EOF) {

    str_buf[l] = (char)c;

    if ((str_buf[l] == '\n') && (l > 0)) { /* dont send blank lines */
      /* printf( "%s",str_buf ); */
      if (strcmp(DAT_DELIM, str_buf) <= 0) {
        l = 0; /* index for string matches */
        m = 0;
        j=2;
        /* printf( "Found the delimiter\n" ); */
        /* sleep(1); */
      }
      else{
        if ((j == 0) && (isempty(str_buf,l) != 1)) {

          /* extract energy eigenvalues and state indexes */
          get_numsl(str_buf,num_idxs1,l,n_idxs1,&e_eigval[m]);
          /* printf( "e_eigval[%d] = %le\n", m+1, e_eigval[m]); */

          m++;
        }
        else if ((j == 2) && (isempty(str_buf,l) != 1)) {
          /* extract transition moments and transition indexes */

          get_numsl(str_buf,num_idxs2,l,n_idxs2,&trans_idxs[0][m],\
                    &trans_idxs[1][m],&t_mom[m]);
          if (trans_idxs[1][m] >= 1000) {
            /* printf( "to %le from %le, %le \n", trans_idxs[0][m], trans_idxs[1][m], t_mom[m]); */
            /* sleep(1); */
            /* exit(1); */
          }

          m++;
        }
        l=0; /* reset the buffer write head to start reading a the next line */
      }
    }
    else{
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

  return init_data_branch(state_er, parsed_input,n_states,n_trans,fn_infile);

}

int
parse_input (double * state_er,
             char * fn_infile, /* name of input file */
             int len_fn){

  int j, k, l; /* looping variables */
  int n_states,n_trans;

  char format[BUF_SIZE] = {0};
  estate root_e_state;

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
        parse_input_molcas(state_er, fn_infile);

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
