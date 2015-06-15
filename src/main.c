#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include "std_char_ops.h"
#include "get_numsl.h"
#include "smap.h"
#include "input_formats.h"
#include "parse_input.h" /* for input_data array, as well as other
                            calculation-specific variables */
#include "smap_cfg.h"
#define BUF_SIZE 256

int
parse_molout (char * fn_relpath,
              char * fn_infile,
              int len_infile
              ) {

  int j,k,l,m; /* control loop variables */
  int mode; /* string matching mode flag */

  FILE * fp_tmpdata;
  FILE * fp_relpath;

  const int n_lookup_str = 6; /* number of strings used for searching the input file */

  int c; /* temporary char for storing input file characters */
  /* we are looking for four strings in the output file: one at the beginning
     of each data block, and one at the end. */
  char * str_buf = malloc(BUF_SIZE*5);
  const char DAT_DELIM[32] =  "============ data block finish\n";

  const char s1[26] = "        Relative EVac(au)";
  const char s2[51] = "Weights of the five most important spin-orbit-free";
  const char s3[40] = "Dipole transition strengths (SO states)";
  const char s4[44] = "Quadrupole transition strengths (SO states)";

  const char s5[36] = "         To  From     Osc. strength";
  const char s6[8] = "the end";

  const char * lookup_str[6] = {s1,s2,s3,s4,s5,s6};

  char outpath[10] = "../output";
  char tmpformat[5] = ".tmp";
  int len_op = 9;
  int len_tf = 4;

  char * tmp_fpstr = malloc(len_op+len_tf+len_infile+1);


  for (j=0; j<len_op; j++) {
    tmp_fpstr[j] = outpath[j];
  }

  tmp_fpstr[j++] = '/';

  /* for (k=0; j<len_op+len_infile; j++) { */
  for (k=0; fn_infile[k] != '\0'; j++) {
    tmp_fpstr[j] = fn_infile[k++];
  }

  for (k=0; tmpformat[k] != '\0'; j++) {
    tmp_fpstr[j] = tmpformat[k++];
  }
  tmp_fpstr[j] = '\0';

  /* open the input file */
  if((fp_relpath = fopen(fn_relpath, "r")) == NULL) {
    fprintf(stderr,"parse_input: unable to open the input file %s.\n",fn_relpath);
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }
  if((fp_tmpdata = fopen(tmp_fpstr, "w+")) == NULL) {
    fprintf(stderr,"parse_input: unable to open the output file %s.\n",tmp_fpstr);
    printf( "program terminating due to the previous error.\n");
    exit(1);
  }

  k = 0; /* index for tmp_string */
  l = 0; /* index for lookup string */
  mode = 0; /* start of in string search mode */

  /* read the Molcas input file */
  for (j=0; ((c = fgetc(fp_relpath)) != EOF); j++, k++) {
    str_buf[k] = (char)c;

    /* keep extracting characters from the input data until an entire line
       has been stored in the temporary str_buf buffer */
    if (str_buf[k] == '\n') {
      if (mode == 1) {
        if ((isanyalpha(str_buf,k) == 0) &&
            (isdashes(str_buf,k) == 0) &&
            (isempty(str_buf,k) == 0)){
          for (m=0; m<=k; m++) {
            fputc(str_buf[m],fp_tmpdata);
          }
          /* printf( "%s\n",str_buf ); */
        } /* else { */
        /*   printf( "%s\n",str_buf ); */
        /*   sleep(1); */
        /* } */
      }
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
          l++;
          mode = bin_flip(mode);
          fprintf(fp_tmpdata,DAT_DELIM );
        }
      }
      k = 0;
    }
  }
  fclose(fp_tmpdata);
  fclose(fp_relpath);
}


int
main (int argc, char * argv[]) {
  printf( "\n\n smap calculation initiated.\n\n" );
  int j,k,l; /* iteration variables */
  int n_t; /* number of numbers read from input flag */
  int n_er; /* number of numbers eigenstate energy values provided */
  int n_res;
  int len_fn = 0;
  int len_mn = 0;
  int dbg_flag = 0;
  double * res = malloc(2*sizeof(double));
  /* arrays for storing input file name data */
  char * input_sbuff = malloc(BUF_SIZE);
  char * fn_infile;
  char * fn_relpath;
  char method[5] = "test";
  char * num_buf;

  char format[BUF_SIZE] = {0};

  /* double state_t[3] = {0.005, 0.001, 0.001}; */

  /* use first element to specify number of values stored in matrix */
  /* threshold values for the three or six different state types */
  double * state_t = malloc(7*sizeof(double));
  /* initial/final and intermediate state energy ranges */
  double * state_er = malloc(9*sizeof(double));

  state_t[0] = 0;
  state_er[0] = 0;

  /* char fn_relpath[BUF_SIZE] = {0}; */

  /* process the input arguments */
  if (argc == 1) {
    printf("No command line arguments provided. To read the help documentation, provide the argument \"-h\". Program terminating.\n");
    fflush(stdout);
  }

  while (argc > 1 && (argv[1][0] == '-')) {

    switch(argv[1][1]){

    case 'd' :
      dbg_flag = 1;
      break;

    case 'h' :
      printf("some help was requested. \n");
      exit(1);
      break;

    case 'i' :
      /* printf( "processing input file: " ); */

      k = 0;
      /* extract the file name */
      for (j=3;argv[1][j] != '\0'; j++) {
        if (argv[1][j] == '/') {
          k = j+1;
        }
      }

      fn_infile = malloc(j-k+1);

      for (l=0,j=k; argv[1][j] != '\0'; j++){

        if (argv[1][j] == '.') {
          j--;
          break;
        }

        fn_infile[l] = argv[1][j];
        l++;
      }

      for (j=3; argv[1][j] != '\0'; j++) {
        input_sbuff[j-3] = argv[1][j];
        /* printf( "%c", input_sbuff[j-3]); */
      }

      /* printf( "\n" ); */
      len_fn = j-3;
      fn_relpath = malloc(len_fn+1);

      for (j=0; j<len_fn; j++) {
        fn_relpath[j] = input_sbuff[j];
        /* printf( "%c", fn_relpath[j]); */
      }
      /* printf( "\n" ); */
      fn_relpath[len_fn] = '\0';

      break;

    /* case 'm' : */
    /*   /\* the user specified a method to be used for calculating the scattering map *\/ */
    /*   for (j=3; argv[1][j] != '\0'; j++) { */
    /*     input_sbuff[j-3] = argv[1][j]; */
    /*   } */
    /*   argv[1] + j; */
    /*   len_mn = j-3; */
    /*   method = malloc(len_mn+1); */

    /*   for (j=0; j<len_mn; j++) { */
    /*     method[j] = input_sbuff[j]; */
    /*   } */
    /*   method[len_mn] = '\0'; */
    /*   printf( "got the method: %s\n", method); */

    /*   break; */

    case 'r' :
      n_res = 0;
      /* the user specified a method to be used for calculating the scattering map */
      for (j=3,k=0; argv[1][j] != '\0'; j++) {

        input_sbuff[k++] = argv[1][j];
        /* if ((argv[1][j] == ',') || ((argv[1][j+1] == '\0')||(argv[1][j+1] == '-'))) { */

        if ((argv[1][j] == ',') || (argv[1][j+1] == '\0')) {
          /* k--; /\* we dont want the comma or null sent to atof *\/ */
          num_buf = malloc(k+1);

          for (l=0; l<k; l++) {
            num_buf[l] = input_sbuff[l];
            /* printf( "%c", num_buf[l]); */
          }
          num_buf[k] = '\0';

          /* printf( "\n" ); */

          res[n_res] = atof(num_buf);
          free(num_buf);
          num_buf = NULL;
          n_res++;
          /* if ((n_res > 2)) { */
          /*   /\* switch on the dipole+quadrupole flag  *\/ */
          /*   break; */
          /* } */
          /* if ((n_res > 8)) { */
          /*   break; */
          /* } */

          k = 0;
        }
      }

      break;

    case 'e' :

      n_er = 1;
      /* the user specified a method to be used for calculating the scattering map */
      for (j=3,k=0; argv[1][j] != '\0'; j++) {

        input_sbuff[k++] = argv[1][j];
        /* if ((argv[1][j] == ',') || ((argv[1][j+1] == '\0')||(argv[1][j+1] == '-'))) { */

        if ((argv[1][j] == ',') || (argv[1][j+1] == '\0')) {
          /* k--; /\* we dont want the comma or null sent to atof *\/ */
          num_buf = malloc(k+1);

          for (l=0; l<k; l++) {
            num_buf[l] = input_sbuff[l];
            /* printf( "%c", num_buf[l]); */
          }
          num_buf[k] = '\0';

          /* printf( "\n" ); */

          state_er[n_er] = atof(num_buf);
          free(num_buf);
          num_buf = NULL;
          n_er++;
          /* if ((n_er > 4)) { */
          /*   /\* switch on the dipole+quadrupole flag  *\/ */
          /*   break; */
          /* } */
          /* if ((n_er > 8)) { */
          /*   break; */
          /* } */

          k = 0;
        }
      }
      state_er[0] = n_er-1;

      break;

    case 't' :
      n_t = 1;
      /* the user specified a method to be used for calculating the scattering map */
      for (j=3,k=0; argv[1][j] != '\0'; j++) {

        input_sbuff[k++] = argv[1][j];
        /* if ((argv[1][j] == ',') || ((argv[1][j+1] == '\0')||(argv[1][j+1] == '-'))) { */

        if ((argv[1][j] == ',') || (argv[1][j+1] == '\0')) {
          /* k--; /\* we dont want the comma or null sent to atof *\/ */
          num_buf = malloc(k+1);

          for (l=0; l<k; l++) {
            num_buf[l] = input_sbuff[l];
            /* printf( "%c", num_buf[l]); */
          }
          num_buf[k] = '\0';

          /* printf( "\n" ); */

          state_t[n_t] = atof(num_buf);
          free(num_buf);
          num_buf = NULL;

          if ((n_t++) > 3) {
            break;
          }

          k = 0;
        }
      }
      state_t[0] = n_t;
      break;

    default :
      printf("Can not process the flag \"-%c\". To read the help documentation, \
 call tau with the flag \"-h\". Program terminating.\n",(char)argv [1][1]);
      exit(1);
    }
    argv++;
    argv++;
    argc--;
    argc--;
  }

  j = len_fn;
  while(fn_relpath[--j] != '.'){};

  for (k=0; j<=len_fn; j++) {
    format[k] = fn_relpath[j];
    k++;
  }
  format[k] = '\0';

  if (strcmp(format,MOLCAS_FORMAT) <= 0) {
    /* reduce the molcas output to a temp file */
    printf( "  - reducing the molcas logfile to a suitable input format ..");
    parse_molout(fn_relpath, fn_infile, len_fn-1);
    printf( " done.\n" );

    goto dealloc;
  }

  if (len_fn == 0) {
    fprintf(stderr, "\n\Error: smap.c, main: you didnt provide the path to an\
 input file.");
    printf( "program terminating due to the previous error.\n");
    exit(EXIT_FAILURE);
  } else {
    /* extract the needed data from the input */
    printf( "  - extracting data from the input file ..");
    fflush(stdout);
    if (parse_input(state_er, fn_relpath, len_fn+1)) {
      fprintf(stderr, "smap.c, main: unable to parse the input data \
contained in %s.\n",fn_relpath);
      printf( "program terminating due to the previous error.\n");
      exit(EXIT_FAILURE);
    } else {
      printf( " done.\n" );
    }
  }

  printf( "\n executing smap with the following..\n\n" );
  printf( "  - data contained in the input file:\n    %s\n\n", fn_relpath);
  printf( "  - threshold values:\n    " );

  for (j=1; j<n_t; j++) {
    printf( "%le, ", state_t[j]);
  }
  printf( "\n\n" );

  printf( "  - state energy intervals values:\n    " );

  for (j=1; j<n_er; j++) {
    printf( "%le, ", state_er[j]);
  }

  printf( "\n\n" );
  printf( " execution progress:\n\n");

  calc_smap_m(fn_relpath,fn_infile,len_fn-1,state_er, state_t, res);
  write_log(state_er, state_t, res, fn_relpath, fn_infile,len_fn-1, 5);

 dealloc:

  free(fn_relpath);
  free(fn_infile);
  free(input_sbuff);
  free(state_er);
  free(state_t);
  printf( "\n smap successfully executed.\n" );
  printf( " program terminating.\n\n" );

  return 0;
}
