int
parse_molout (double * state_er,
              int * mom,
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

  const int n_lookup_str = 3; /* number of strings used for searching the input file */
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

  const char * lookup_str[3] = {s1,s2,s3,s4};

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

    /* keep extracting characters from the input data until an entire line
       has been stored in the temporary str_buf buffer */
    if (str_buf[k] == '\n') {

      /* check every line for a matching substring */
      /* mode = 1 and a line match means that we reached
         the end of this data block */
      if (strstr(str_buf,lookup_str[l])) {

        /* we found the first substring, the data we're looking for is
           inside the coming block of text. switch to mode 1.*/
        if (mode == 0) {
          l++;
          mode = bin_flip(mode);
        } else {
          mode = bin_flip(mode);
        }
      }
      if (mode == 1) {
        fprintf(fp_tmpdata,DAT_DELIM, n_matchfound );
      }
    }
  }
  k = 0;
}
