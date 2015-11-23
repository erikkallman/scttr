/* Copyright (C) 2015 Erik Källman */
/* This file is part of the scttr program. */

/* scttr is free software: you can redistribute it and/or modify */
/* it under the terms of the GNU Lesser General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or */
/* (at your option) any later version. */

/* scttr is distributed in the hope that it will be useful, */
/* but without any warranty; without even the implied warranty of */
/* merchantability or fitness for a particular purpose. See the */
/* GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License */
/* along with scttr, found in the "license" subdirectory of the root */
/* directory of the scttr program. */
/* If not, see <http://www.gnu.org/licenses/>. */
/**
   * @file scttr_io.c
   * @author Erik Källman
   * @date November 2015
   * @brief Defines all functions that relate directly to the manipulation or
   * storage of user i/o.
   */
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#include <string.h>
#include <fcntl.h>
#include <errno.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "dyn_array.h"
#include "sci_const.h"
#include "formats.h"
#include "iquicks.h"
#include "std_char_ops.h"
#include "std_num_ops.h"
#include "get_nums.h"
#include "spectrum.h"
#include "spectrum_s.h"
#include "inp_node_s.h"
#include "metadata_s.h"
#include "transitions.h"
#include "scttr_cfg.h"

struct inp_node *root_inp;
int n_inp = 0;

/* suffixes for the output and input files */
const char *dat_sfx = ".dat";
const char *plot_sfx = ".gp";
const char *log_sfx  = ".txt";
const char *bin_sfx  = ".bin";
const char *tmp_sfx  = ".tmp";

struct metadata *
init_md ()
{
  struct metadata *new_md;

  if((new_md = malloc(sizeof(struct metadata))) == NULL ) {
    fprintf(stderr, "\n\nscttr_io.c, function init_md: failed to allocate memory for \"new_md\"\n");
    printf("program terminating due to the previous error.\n\n");
    exit(1);
  }

  return new_md;
}

int
free_md (struct metadata *md)
{
  free(md -> outpath);
  free(md -> inpath);
  free(md -> inp_fn);
  free(md);
  return 1;
}

struct inp_node *
get_inp (char *id)
{
  struct inp_node *curr_inp = root_inp;
  struct inp_node *next_inp;

  while(strstr((curr_inp -> md -> inp_fn), id) == NULL) {

    next_inp = curr_inp -> next;
    curr_inp = next_inp;
    if (curr_inp == NULL) {
      fprintf(stderr, "\n\nscttr_io.c, function get_inp: no inp_node node can be found with a str_id == %s\n"
              , id);
      printf("program terminating due to the previous error.\n\n");
      exit(1);
    }
  }

  return curr_inp;
}

struct inp_node *
init_inp (struct metadata *md)
{
  /* in order not to have to look through the input node linked list, store a
     static pointer to the last node in the lists */
  static struct inp_node *last_inp;

  struct inp_node *new_inp;

  if((new_inp = malloc(sizeof(struct inp_node))) == NULL ) {
    fprintf(stderr, "\n\nscttr_io.c, function init_inp: failed to allocate memory for \"new_inp\"\n");
    printf("program terminating due to the previous error.\n\n");
    exit(1);
  }

  new_inp -> md = md;

  if (n_inp == 0) {
    new_inp -> last = NULL;
    root_inp = new_inp;
    last_inp = root_inp;
  }
  else {
    new_inp -> last = last_inp;
    last_inp -> next = new_inp;
  }

  new_inp -> next = NULL;
  new_inp -> idx = n_inp++;
  new_inp -> root_spec = NULL;

  return new_inp;
}

int
free_inp (struct inp_node *inp)
{
  int j;

  if (n_inp > 1) {
    inp -> last -> next = inp -> next;
    inp -> next -> last = inp -> last;
  }

  for (j = 0; j < 6; j++) {
    free(inp -> trs[j]);
  }

  free(inp -> trs);
  free(inp -> idx_map);
  free_md(inp -> md);
  free_all_specs(inp);
  free(inp);

  return 1;
}

int
parse_input_bin (struct inp_node *inp, char *bin_fpstr)
{
  int j;
  int pi_xdim, pi_ydim;

  printf("\n      processing binary file corresponding to the provided input file ..");
  fflush(stdout);

  FILE *fp_bin = fopen(bin_fpstr, "rb");

  if ( fread(&pi_ydim, sizeof(int), 1, fp_bin) != 1) {
    fprintf(stderr, "\n\nscttr_io.c, function parse_input_bin: unable to read pi_ydim from binary output file \n");
    printf("program terminating due to the previous error.\n\n");
    exit(1);
  }

  pi_xdim = 6;
  inp -> n_trans = pi_ydim;

  if((inp -> trs = malloc(pi_xdim * sizeof(double *))) == NULL ) {
    fprintf(stderr, "\n\nscttr_io, function parse_input_bin: failed to allocate memory for \"input_data\"\n");
    printf("program terminating due to the previous error.\n\n");
    exit(1);
  }

  for (j = 0; j < pi_xdim; j++) {
    if((inp -> trs[j] = malloc((pi_ydim + 1) * sizeof(double))) == NULL ) {
      fprintf(stderr, "\n\nscttr_io.c, function parse_input_bin: failed to allocate memory for pointers in \"input_data\"\n");
      printf("program terminating due to the previous error.\n\n");
      exit(1);
    }
  }

  for (j = 0; j < 6; j++) {
    if (fread(inp -> trs[j], sizeof(double), inp -> n_trans, fp_bin)
        != inp -> n_trans) {
      fprintf(stderr, "\n\nscttr_io.c, function parse_input_bin: unable to read the PI matrix\n");
      printf("program terminating due to the previous error.\n\n");
      fclose(fp_bin);
      exit(1);
    }
  }

  inp -> e0 = inp -> trs[2][j];

  for (j = 0; j < pi_ydim; j++) {
    if (inp -> trs[2][j] < inp -> e0) {
      inp -> e0 = inp -> trs[2][j];
    }
  }

  inp -> trs[0][pi_ydim] = -1;

  if (fclose(fp_bin) != 0) {
    fprintf(stderr, "\n\nscttr_io.c, function parse_input_bin: unable to close file:\n%s\n"
            , bin_fpstr);
    printf("program terminating due to the previous error.\n\n");
    exit(1);
  }
  printf(" done\n");
  return 1;
}

int
parse_molout (struct inp_node *inp, char *fn_relpath, char *tmp_fpstr)
{
  printf("\n      parsing the molcas .log file .. \n");
  int j,k,l,m;

  /* the mode flag determining if the function should be searching for
     strings with numbers (1) or strings matching those in s1-3 (0) */
  int mode;
  int string_flag = 0;
  int n_match = 0;

  FILE *fp_tmpdata;
  FILE *fp_relpath;

  int c; /* temporary char for storing input file characters */

  char *s1 = NULL;
  char *s2 = NULL;
  char *s3 = NULL;

  char *str_buf = malloc(BUF_SIZE * 5);
  if (inp -> md -> so_enrg == 0) {

    /* read spin-orbit data */
    s1 = malloc(35);
    s2 = malloc(40);
    s3 = malloc(44);

    s1 = "Eigenvalues of complex Hamiltonian";
    s2 = "Dipole transition strengths (SO states)";
    s3 = "Quadrupole transition strengths (SO states)";
  }
  else if(inp -> md -> so_enrg == 1) {

    /* read spin-orbit free data */
    s1 = malloc(19);
    s2 = malloc(29);
    s3 = malloc(32);

    s1 = "SPIN-FREE ENERGIES";
    s2 = "Dipole transition strengths:";
    s3 = "Quadrupole transition strengths:";
  }

  const char *lookup_str[3] = {s1, s2, s3};

  if((fp_relpath = fopen(fn_relpath, "r")) == NULL) {
    fprintf(stderr, "\n\nscttr_io.c, function parse_molout: unable to open the input file %s.\n"
            , fn_relpath);
    printf("program terminating due to the previous error.\n\n");
    exit(1);
  }
  if((fp_tmpdata = fopen(tmp_fpstr, "w+")) == NULL) {
    fprintf(stderr, "\n\nscttr_io.c, function parse_molout: unable to open the output file %s.\n"
            , tmp_fpstr);
    printf("program terminating due to the previous error.\n\n");
    exit(1);
  }

  mode = 0;
  string_flag = 0;

  for (j = 0, k = 0; ((c = fgetc(fp_relpath)) != EOF) && (n_match != 3); j++, k++) {
    if (j % 100000 == 0) {
      printf("        %.2f%%\r", (((float)(j * sizeof(char))
                                   / (float)inp -> md -> sz_inp) * 100));
    }

    str_buf[k] = (char)c;

    /* keep extracting characters from the input data until an entire line
       has been stored in the temporary str_buf buffer */
    if (str_buf[k] == '\n') {
      if (mode == 1) {

        /* check if the string that was read only contains numbers */
        if ((isanyalpha(str_buf, k) == 0) &&
            (isdashes(str_buf, k) == 0) &&
            (isempty(str_buf, k) == 0)) {
          if (string_flag == 0) {

            /* we are now reading numbers. accept no strings */
            string_flag =  BIN_FLIP(string_flag);
          }
          for (m = 0; m <= k; m++) {
            fputc(str_buf[m], fp_tmpdata);
          }
        }

        /* if we find a flag while string_flag ==1 and mode ==1, we jhave
           read beyond the table */
        else if(string_flag == 1) {
          string_flag = BIN_FLIP(string_flag);

          /* something was read that was not a number.
             flip back to searching for new line matches */
          mode = BIN_FLIP(mode);
          n_match++;
        }
      }

      /* check every line for a matching substring */
      /* mode = 1 and a line match means that we reached
         the end of this data block */
      else if((l = strscmp(str_buf, lookup_str, 3)) != -1) {

        /* we found the first substring, the data we're looking for is
           inside the coming table of numbers of text. switch to mode 1.*/
        fprintf(fp_tmpdata,"%s\n", lookup_str[l]);
        mode = BIN_FLIP(mode);
      }
      k = 0;
    }
  }
  printf("          100%%\r");
  if ((fclose(fp_tmpdata) != 0)
      && (fclose(fp_relpath) != 0)) {
    fprintf(stderr, "\n\nscttr_io.c, function parse_molout: unable to close files:\n%s\n%s\n"
            , tmp_fpstr,fn_relpath);
    printf("program terminating due to the previous error.\n\n");
    exit(1);
  }
  free(str_buf);
  printf("\n      done. \n");
  return 1;
}

int
parse_input_tmp (struct inp_node *inp, char *fn_tmpdata)
{
  printf("\n      parsing the tmp file .. \n");
  int j, k, l, m, j_test; /* control loop variables */
  int idx_from, idx_to; /* index in each transition from state x to y  */
  int tmp_idx2;
  int n_states, n_trans;
  int last_i;
  int n_proc; /* number of states processed at any given time in the
               execution of the function */
  int trs_type; /* type  */

  int n_idxs1 = 2;
  int n_idxs2 = 3;

  int *num_idxs1;
  int *num_idxs2;
  int *trs_types;
  int *idxs_eigval;
  int *proc_idx;

  double *state_er = inp -> md -> state_er;
  double tmp_idx = 0;
  double maxr = get_maxl(state_er, state_er[0]);
  double from_state_en, to_state_en, tmp_state_en;

  double *e_eigval;
  double *t_mom;

  double **trs_buf;
  double **trans_idxs;

  /* variables used later to approximate the amount of input data */
  int sz_tmp;
  struct stat st = {0};

  FILE *fp_tmpdata;

  int c; /* temporary char for storing input file characters */
  /* we are looking for four strings in the output file: one at the beginning
     of each data block, and one at the end. */
  char *str_buf = malloc(BUF_SIZE * 5);
  char *s1 = NULL;
  char *s2 = NULL;
  char *s3 = NULL;

  if (inp -> md -> so_enrg == 0) {
    /* read spin-orbit data */
    s1 = "Eigenvalues of complex Hamiltonian";
    s2 = "Dipole transition strengths (SO states)";
    s3 = "Quadrupole transition strengths (SO states)";
  }
  else if(inp -> md -> so_enrg == 1) {
    /* read spin-orbit free data */
    s1 = "SPIN-FREE ENERGIES";
    s2 = "Dipole transition strengths";
    s3 = "Quadrupole transition strengths";
  }

  /* create a pointer to the three data block beginners s1,s3, */
  const char *lookup_str[3] = {s1,s2,s3};

  stat(fn_tmpdata, &st);
  sz_tmp = (int)st.st_size; /* file size in bytes */
  sz_tmp = sz_tmp / 34; /* shortest line = 34 chars */

  /* open the input file */
  if((fp_tmpdata = fopen(fn_tmpdata, "r")) == NULL) {
    fprintf(stderr, "\n\nscttr_io.c, function parse_input_tmp: unable to open the input file %s.\n"
            , fn_tmpdata);
    printf("program terminating due to the previous error.\n\n");
    exit(1);
  }

  k = 0; /* index for tmp_string */
  l = 0; /* index for lookup string */

  /* storage for the energy eigenvalues */
  if((e_eigval = malloc(sz_tmp * sizeof(double))) == NULL ) {
    fprintf(stderr, "\n\nscttr_io.c, function parse_input_tmp: failed to allocate memory for \"e_eigval\"\n");
    printf("program terminating due to the previous error.\n\n");
    exit(1);
  }

  if((trs_types = malloc(sz_tmp * sizeof(double))) == NULL ) {
    fprintf(stderr, "\n\nscttr_io.c, function parse_input_tmp: failed to allocate memory for \"trs_types\"\n");
    printf("program terminating due to the previous error.\n\n");
    exit(1);
  }

  if((idxs_eigval = malloc(sz_tmp * sizeof(int))) == NULL ) {
    fprintf(stderr, "\n\nscttr_io.c, function parse_input_tmp: failed to allocate memory for \"idxs_eigval\"\n");
    printf("program terminating due to the previous error.\n\n");
    exit(1);
  }

  for (j = 0; j < sz_tmp; j++) {
    idxs_eigval[j] = -1;
  }

  /* storage for the transition moments */
  if((t_mom = malloc(sz_tmp * sizeof(double))) == NULL ) {
    fprintf(stderr, "\n\nscttr_io.c, function parse_input_tmp: failed to allocate memory for \"t_mom\"\n");
    printf("program terminating due to the previous error.\n\n");
    exit(1);
  }

  /* storage for the transition indexes, column 1 is from a state
     column 2 is to state index */
  if((trans_idxs = malloc(2 * sizeof(double *))) == NULL ) {
    fprintf(stderr, "\n\nscttr_io.c, function parse_input_tmp: failed to allocate memory for \"trans_idxs\"\n");
    printf("program terminating due to the previous error.\n\n");
    exit(1);
  }

  for (j = 0; j < 2; j++) {
    if((trans_idxs[j] = malloc(sz_tmp * sizeof(double))) == NULL ) {
      fprintf(stderr, "\n\nscttr_io.c, function parse_input_tmp: failed to allocate memory for \"trans_idxs[%d]\"\n"
              , j);
      printf("program terminating due to the previous error.\n\n");
      exit(1);
    }
  }

  if((num_idxs1 = malloc(2 * sizeof(int))) == NULL ) {
    fprintf(stderr, "\n\nscttr_io.c, function parse_input_tmp: failed to allocate memory for \"num_idxs1\"\n");
    printf("program terminating due to the previous error.\n\n");
    exit(1);
  }

  if((num_idxs2 = malloc(3 * sizeof(int))) == NULL ) {
    fprintf(stderr, "\n\nscttr_io.c, function parse_input_tmp: failed to allocate memory for \"num_idxs2\"\n");
    printf("program terminating due to the previous error.\n\n");
    exit(1);
  }

  num_idxs1[0] = 0;
  num_idxs1[1] = 1;

  num_idxs2[0] = 0;
  num_idxs2[1] = 1;
  num_idxs2[2] = 2;

  j = j_test = l = m = 0;

  n_states = 0;
  n_trans = 0;
  trs_type = 1;

  /* now that data structures of the rights size has memory allocated for
     them, start reading data from the temporary file */
  while ((c = fgetc(fp_tmpdata)) != EOF) {

    str_buf[l] = (char)c;

    if ((str_buf[l]  == '\n') && (l > 0)) { /* dont send blank lines */
      str_buf[l + 1] = '\0';
      if ((j_test = strscmp(str_buf, lookup_str, 3)) != -1) {
        j = j_test;
        trs_type = j;
      }
      else {
        if ((j == 0) && (isempty(str_buf, l) != 1)) {

          /* extract energy eigenvalues and state indexes */
          get_nums(str_buf, num_idxs1, l, n_idxs1, &tmp_idx
                   , &e_eigval[n_states]);
          if ((fabs(e_eigval[n_states] - e_eigval[0]) * AUTOEV) \
              < maxr) {
            idxs_eigval[(int)tmp_idx - 1] = n_states + 1;
            n_states++;
          }
        }
        else if ((j > 0) && (isempty(str_buf, l) != 1)) {
          if (n_trans == 0) { /* pre-process the eigenvalues */

            /* sort the states in energy */
            iquicks(e_eigval, idxs_eigval, 0, n_states - 1, n_states);
            inp -> e0 = e_eigval[0];

            /* adjust n_states so that it accounts for states not to be read
             due to being outside of the input range */
            m = 0;
            for (k = 0; k < n_states; k++) {
              if ((e_eigval[k] - inp -> e0) * AUTOEV < maxr) {
                m++;
              }
            }
            n_states = m;
          }
          /* extract transition moments and transition indexes */
          get_nums(str_buf, num_idxs2, l, n_idxs2, &trans_idxs[0][n_trans],
                   &trans_idxs[1][n_trans], &t_mom[n_trans]);
          trs_types[n_trans] = trs_type;

          from_state_en = get_wi(e_eigval, idxs_eigval
                                 , (int)(trans_idxs[1][n_trans]), n_states);

          to_state_en = get_wi(e_eigval,idxs_eigval
                               , (int)(trans_idxs[0][n_trans]), n_states);

          if ((((int)to_state_en != -1) && ((int)from_state_en != -1))
              && (((to_state_en-inp -> e0) * AUTOEV < maxr)
                  && ((from_state_en-inp -> e0) * AUTOEV < maxr))){
            n_trans++;
          }
        }
      }
      /* reset the buffer write head to start reading a the next line */
      l = 0;
    }
    else{
      l++;
    }
  }

  /* allocate space for the "parsed input matrix" that will be filled with data
     in the remaining sections of this function */

  if((proc_idx = malloc(n_states * sizeof(int))) == NULL ) {
    fprintf(stderr, "\n\nscttr_io.c, function parse_input_tmp: failed to allocate memory for \"proc_idx\"\n");
    printf("program terminating due to the previous error.\n\n");
    exit(1);
  }

  if((inp -> trs = malloc(6 * sizeof(double *))) == NULL ) {
    fprintf(stderr, "\n\nscttr_io.c, function parse_input_tmp: failed to allocate memory for \"inp -> trs\"\n");
    printf("program terminating due to the previous error.\n\n");
    exit(1);
  }

  for (j = 0; j < 6; j++) {
    if((inp -> trs[j] = malloc((n_trans + 1) * sizeof(double))) == NULL ) {
      fprintf(stderr, "\n\nscttr_io.c, function parse_input_tmp: failed to allocate memory for pointers in \"inp -> trs[%d]\"\n"
              , j);
      printf("program terminating due to the previous error.\n\n");
      exit(1);
    }
  }
  inp -> trs[0][n_trans] = -1; /* end of list */

  if((trs_buf = malloc(6 * sizeof(double *))) == NULL ) {
    fprintf(stderr, "\n\nscttr_io.c, function parse_input_tmp: failed to allocate memory for \"trs_buf\"\n");
    printf("program terminating due to the previous error.\n\n");
    exit(1);
  }

  for (j = 0; j < 6; j++) {
    if((trs_buf[j] = malloc((n_trans + 1) * sizeof(double))) == NULL ) {
      fprintf(stderr, "\n\nscttr_io.c, function parse_input_tmp: failed to allocate memory for pointers in \"trs_buf[%d]\"\n"
              , j);
      printf("program terminating due to the previous error.\n\n");
      exit(1);
    }
  }

  /* finally, store the data in the trs matrix for the parse_input function  */
  n_proc = 0;
  last_i = 1;
  j = k = l = 0;

  /* l = index of line in the data arrays */
  /* j = index of data read into trs_buf, reset upon reading the transitions
     from a new state */

  while (l < n_trans - 1) {
    if (l % 1 == 0) {
      printf("        %.2f%%\r", (((float)l / (float)n_trans) * 100));
    }
    if ((j + l) == (n_trans)) {
      /* dont read beyond the last value in PI */
      trs_buf[0][j] = -1;
      idx_from = -1;
    }
    else {
      from_state_en = get_wi(e_eigval, idxs_eigval
                             , (int)(trans_idxs[0][j + l]), n_states);
      to_state_en = get_wi(e_eigval, idxs_eigval
                           , (int)(trans_idxs[1][j + l]), n_states);

      if (inrange((to_state_en - inp -> e0) * AUTOEV, state_er[3]
                  , state_er[4])
          && (inrange((from_state_en - inp -> e0)*AUTOEV,state_er[5]
                      ,state_er[6]))
          && (((int)state_er[1] != (int)state_er[5])
              && ((int)state_er[2] != (int)state_er[6]))
          ) {

        idx_from = trans_idxs[1][j + l];
        idx_to = trans_idxs[0][j + l];

        tmp_state_en = from_state_en;
        from_state_en = to_state_en;
        to_state_en = tmp_state_en;
      }
      else {

        idx_from = trans_idxs[0][j + l];
        idx_to   = trans_idxs[1][j + l];
      }

      trs_buf[0][j] = idx_from;
      trs_buf[1][j] = idx_to;
      trs_buf[2][j] = from_state_en;
      trs_buf[3][j] = to_state_en;
      trs_buf[4][j] = t_mom[j + l];
      trs_buf[5][j] = trs_types[j + l];
    }
    if (idx_from != last_i) {

      tmp_idx2 = get_inext(inp -> trs, last_i);
      /* we have read all transitions for a state */
      /* check if the last_i has already been processed */
      if (intinint(proc_idx, last_i, n_proc) == -1) {

        m = l;
        while ((int)trs_buf[0][m - l] != idx_from) {

          inp -> trs[0][m] = trs_buf[0][m - l];
          inp -> trs[1][m] = trs_buf[1][m - l];
          inp -> trs[2][m] = trs_buf[2][m - l];
          inp -> trs[3][m] = trs_buf[3][m - l];
          inp -> trs[4][m] = trs_buf[4][m - l];
          inp -> trs[5][m] = trs_buf[5][m - l];

          m++;
        }

        inp -> trs[0][m] = -1;
        proc_idx[n_proc++] = last_i;
      }
      else{

        tmp_idx2 = get_inext(inp -> trs, last_i);
        fwdsplice(trs_buf, inp -> trs, tmp_idx2, l, j, 6);
        fflush(stdout);
      }
      l += j;
      j = 0;

    } else {
      j++;
    }
    last_i = idx_from;
  }

  inp -> trs[0][l] = -1;
  inp -> n_trans = l;
  printf("          100%%\r");

  if (fclose(fp_tmpdata) != 0) {
    fprintf(stderr, "\n\nscttr_io.c, function parse_input_tmp: unable to close file:\n%s\n"
            , fn_tmpdata);
    printf("program terminating due to the previous error.\n\n");
    exit(1);
  }

  free(e_eigval);
  free(trs_types);
  free(idxs_eigval);
  free(t_mom);

  for (j = 0; j < 2; j++) {
    free(trans_idxs[j]);
  }
  free(trans_idxs);

  free(num_idxs1);
  free(num_idxs2);
  free(proc_idx);

  for (j = 0; j<6; j++) {
    free(trs_buf[j]);
  }
  free(trs_buf);
  free(str_buf);

  printf("\n      done.\n");
  return 1;
}

int
parse_input (struct inp_node *inp)
{
  int j;
  int rc;                       /* return code */

  FILE *fp_tmp;
  FILE *fp_bin;

  struct metadata *md = inp -> md;

  char *bin_fpstr = concs(3,md -> outpath,md -> inp_fn, bin_sfx);
  char *tmp_fpstr = concs(3,md -> outpath,md -> inp_fn, tmp_sfx);
  char *format = md -> inp_sfx;
  char *inpath = md -> inpath;

  fp_bin = fopen(bin_fpstr, "r");
  if (fp_bin != NULL) {

    /* process the binary file instead */
    fclose(fp_bin);
    parse_input_bin(inp, bin_fpstr);
  }
  else {

    if (strcmp(format, MOLCAS_FORMAT) <= 0) {

      /* reduce the molcas output to a temp file */
      parse_molout(inp, inpath, tmp_fpstr );

      fp_tmp = fopen(inpath, "r");
      if (fp_tmp == NULL) {

        fprintf(stderr, "\n\nscttr_io.c, function parse_input: unable to locate file %s for processing.\n"
                ,inpath);
        printf("program terminating due to the previous error.\n\n");
        fclose(fp_tmp);
        exit(1);
      } else {

        parse_input_tmp(inp, tmp_fpstr);
      }
    }
    else {

      /* the tmp file path was provided in the input */
      parse_input_tmp(inp, inpath);
    }

    if ((md -> state_er[1] == md -> state_er[5])
        && (md -> state_er[2] == md -> state_er[6])) {

      count_states(inp);
      add_eltrans(inp);
    }

    /* the only way the parse_input_tmp function can get called is if there is if
       no binary file present. we can therefore safely write to the binary file
       without checking if it already exists. */
    if ((fp_bin = fopen(bin_fpstr, "wb")) == NULL) {
      fprintf(stderr, "\n\nscttr_io.c, function parse_input: unable to open the binary file used to store the trs matrix: %s\n"
              , bin_fpstr);
      printf("program terminating due to the previous error.\n\n");
      exit(1);
    }

    fwrite((const void*) & (inp -> n_trans), sizeof(int), 1, fp_bin);
    fflush(fp_bin);
    fclose(fp_bin);

    if ((fp_bin = fopen(bin_fpstr, "ab")) == NULL) {
      fprintf(stderr, "\n\nscttr_io.c, function parse_input: unable to open the binary file used to store the trs matrix: %s\n"
              , bin_fpstr);
      printf("program terminating due to the previous error.\n\n");
      exit(1);
    }

    for (j = 0; j < 6; j++) {
      if (fwrite(inp -> trs[j], sizeof(double), inp -> n_trans, fp_bin) != inp -> n_trans) {
        fprintf(stderr, "\n\nscttr_io.c, function parse_input: function parse_input_bin: unable to write the trs matrix\n");
        printf("program terminating due to the previous error.\n\n");
        fclose(fp_bin);
        exit(1);
      }
      fflush(fp_bin);
    }

    if (fclose(fp_bin) != 0) {
      fprintf(stderr, "\n\nscttr_io.c, function parse_input: unable to close files:\n%s\n"
              , bin_fpstr);
      printf("program terminating due to the previous error.\n\n");
      exit(1);
    }
  }

  if ((rc = set_root_spec(inp)) != 0) {

    fprintf(stderr, "\n\nscttr_io.c, function parse_input: set_root_spec() returned non-zero integer.\n");

    if(rc >= 0) {
      fprintf(stderr, "\nstate %d occured multiple times in the input matrix. Symmetric transitions were erronously added, most likely a failure in the fwdsplice function.\n"
              ,rc);
    }
    printf("program terminating due to the previous error.\n\n");
    exit(EXIT_FAILURE);
  }

  free(bin_fpstr);
  free(tmp_fpstr);

  return 1;
}

int
write_spec (struct inp_node *inp,
            struct spectrum *spec)
{
  int x, y;

  FILE *fp_smat_out;
  char *dat_fpstr = concs(3, inp -> md -> outpath,
                          inp -> md -> inp_fn, dat_sfx);

  if((fp_smat_out = fopen(dat_fpstr, "w"))==NULL) {
    fprintf(stderr, "\n\nscttr_io.c, function write_spec: unable to open the output file %s.\n"
            , dat_fpstr);
    printf("program terminating due to the previous error.\n\n");
    exit(1);
  }

  printf("  - writing RIXS map to file: %s ..", dat_fpstr);
  for (x = 0; x < spec -> n_elx; x++) {
    for (y = 0; y < spec -> n_ely; y++) {
      spec -> s_mat[x][y] = spec -> s_mat[x][y] / spec -> sfac;
      fprintf(fp_smat_out,"%le %le %le\n", (spec -> omega_x[x][y]) * AUTOEV, spec -> omega_y[x][y]
              * AUTOEV, spec -> s_mat[x][y]);
      fflush(fp_smat_out);
    }
    fprintf(fp_smat_out,"\n");
    fflush(fp_smat_out);
  }

  printf(" done.\n");
  if (fclose(fp_smat_out) != 0) {
    fprintf(stderr, "\n\nscttr_io.c, function write_spec:: unable to close some of the files:\n%s\n",
            dat_fpstr);
    printf("program terminating due to the previous error.\n\n");
    exit(1);
  }

  free(dat_fpstr);
  return 1;
}

int
write_plotscript (struct inp_node *inp,
                  struct spectrum *spec)
{

  if (spec == NULL) {
    fprintf(stderr, "scctr_io.c, function write_plotscript: variable spec not defined. \n");
    printf( "program terminating due to the previous error.\n\n");
    exit(1);
  }

  FILE *fp_plot_in;
  FILE *fp_plot_out;
  char *plot_fpstr = concs(3, inp -> md -> outpath,
                           inp -> md -> inp_fn, plot_sfx);

  /* open the placeholder file */
  if((fp_plot_in = fopen(PLOT_TEMPLATE_PTH, "r"))==NULL) {
    fprintf(stderr, "\n\nscttr_io.c, function write_plotscript: unable to open the output file %s.\n"
            ,PLOT_TEMPLATE_PTH);
    printf("program terminating due to the previous error.\n\n");
    exit(1);
  }

  if((fp_plot_out=fopen(plot_fpstr, "w"))==NULL) {
    fprintf(stderr, "\n\nscttr_io.c, function write_plotscript: unable to open the output file %s.\n"
            , plot_fpstr);
    printf("program terminating due to the previous error.\n\n");
    exit(1);
  }

  struct metadata *md = inp -> md;

  size_t len;
  ssize_t read;

  char *line;

  line = NULL;
  len = 0;

  /* construct a gnuplot script from the energy ranges */
  while((read = getline(&line, &len, fp_plot_in)) != -1) {
    fprintf(fp_plot_out, "%s", line);
  }

  fprintf(fp_plot_out, "set output \"./%s.png\"\n", md -> inp_fn);
  fprintf(fp_plot_out, "set title \"%s\" font \"Helvetica,40\" offset 0,1\n"
          , md -> inp_fn);
  fprintf(fp_plot_out, "splot [%le:%le][%le:%le] \"./%s.dat\" u (($1+xshift)/1):2:($3*sc) with pm3d title \"\"",
          spec -> omega_x[0][0] * AUTOEV,
          spec -> omega_x[spec -> n_elx - 1][0] * AUTOEV,
          spec -> omega_y[0][0] * AUTOEV,
          spec -> omega_y[0][spec -> n_ely - 1] * AUTOEV,
          md -> inp_fn);

  if (fclose(fp_plot_in) != 0) {
    fprintf(stderr, "calc_spec.c, function calc_spec: unable to close ../src/plot_template\n");
    printf("program terminating due to the previous error.\n\n");
    exit(1);
  }

  if (fclose(fp_plot_out) != 0) {
    fprintf(stderr, "calc_spec.c, function calc_spec: unable to close some of the files:\n%s\n",
            plot_fpstr);
    printf("program terminating due to the previous error.\n\n");
    exit(1);
  }

  free(line);
  free(plot_fpstr);
  return 1;
}