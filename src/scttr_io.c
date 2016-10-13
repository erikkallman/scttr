/* Copyright (C) 2015 Erik Källman */
/* This file is part of the scttr program. */

/* scttr is free software: you can redistribute it and/or modify */
/* it under the terms of the GNU Lesser General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or */
/* (at your option) any later version. */

/* scttr is distributed in the hope that it will be useful, */
/* but without any warranty; without even the implied warranty of */
/* merchantability or fitness for a particular purpose. See the */
/* GNU Lesser General Public License for more details. */

/* You should have received a copy of the GNU Lesser General Public License */
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
#include <fcntl.h>
#include <errno.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <omp.h>
#include "timing.h"
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
#include "glob_time.h"

double total_t;

struct inp_node *root_inp;
int n_inp = 0;

/* suffixes for the output and input files */
const char *dat_sfx = ".dat";
const char *plot_sfx = ".gp";
const char *log_sfx  = ".txt";
const char *time_sfx  = "_time.txt";
const char *stick_sfx  = "_stick.txt";
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
  free(md -> inp_sfx);
  free(md -> state_er);
  free(md -> state_t);
  free(md -> res);

  free(md -> gx);
  free(md -> gy);
  free(md -> lx);
  free(md -> ly);

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
  inp -> bw_sum = 0;

  printf("\n      processing binary file corresponding to the provided input file ..");
  fflush(stdout);

  FILE *fp_bin = fopen(bin_fpstr, "rb");

  if ( fread(&pi_ydim, sizeof(int), 1, fp_bin) != 1) {
    fprintf(stderr, "\n\nscttr_io.c, function parse_input_bin: unable to read pi_ydim from binary output file \n");
    printf("program terminating due to the previous error.\n\n");
    exit(1);
  }

  if ( fread(&inp -> bw_sum, sizeof(double), 1, fp_bin) != 1) {
    fprintf(stderr, "\n\nscttr_io.c, function parse_input_bin: unable to read bw_sum from binary output file \n");
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
    if ((int)fread(inp -> trs[j], sizeof(double), inp -> n_trans, fp_bin)
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
  /* printf("\n" ); */
  /* for (j = 0; j < pi_ydim; j++) { */
  /*   printf("%le %le\n", (inp -> trs[2][j] - inp -> e0)*AUTOEV, (inp -> trs[3][j] - inp -> e0)*AUTOEV); */
  /*     /\* inp -> bw_sum += inp *\/ */
  /*     /\* inp -> e0 = inp -> trs[2][j]; *\/ */
  /* } */
  /* printf("\n" ); */
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

  int j,k,l,m;

  /* the mode flag determining if the function should be searching for
     strings with numbers (1) or strings matching those in s1-3 (0) */
  int mode;
  int string_flag = 0;
  int n_match = 0;

  FILE *fp_tmpdata;
  FILE *fp_relpath;

  int c; /* temporary char for storing input file characters */

  char *ltime = malloc(20);
  char *s1 = NULL;
  char *s2 = NULL;
  char *s3 = NULL;

  char *str_buf = malloc(BUF_SIZE * 5);
  printf("\n      parsing the molcas .log file (%s) ..", get_loctime(ltime));
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

  if ((fclose(fp_tmpdata) != 0)
      && (fclose(fp_relpath) != 0)) {
    fprintf(stderr, "\n\nscttr_io.c, function parse_molout: unable to close files:\n%s\n%s\n"
            , tmp_fpstr,fn_relpath);
    printf("program terminating due to the previous error.\n\n");
    exit(1);
  }
  free(str_buf);
  printf(" done (%s).\n", get_loctime(ltime));
  free(ltime);
  return 1;
}

/* int */
/* parse_input_tmp_el (struct inp_node *inp, char *fn_tmpdata) */
/* { */
/*   int j, k, l, m, j_test; /\* control loop variables *\/ */
/*   int idx_from, idx_to; /\* index in each transition from state x to y  *\/ */
/*   int tmp_idx2; */
/*   int n_states, n_trans; */
/*   int last_i; */
/*   int n_proc; /\* number of states processed at any given time in the */
/*                execution of the function *\/ */
/*   int trs_type; /\* type  *\/ */

/*   int n_idxs1 = 2; */
/*   int n_idxs2 = 3; */

/*   char *ltime = malloc(20); */

/*   int *num_idxs1; */
/*   int *num_idxs2; */
/*   int *trs_types; */
/*   int *idxs_eigval; */
/*   int *proc_idx; */

/*   double *state_er = inp -> md -> state_er; */
/*   double tmp_idx = 0; */
/*   double maxr = get_maxl(state_er, state_er[0]); */
/*   double from_state_en, to_state_en, tmp_state_en; */

/*   double *e_eigval; */
/*   double *t_mom; */

/*   double **trs_buf; */
/*   double **trans_idxs; */

/*   /\* variables used later to approximate the amount of input data *\/ */
/*   int sz_tmp; */
/*   struct stat st = {0}; */

/*   FILE *fp_tmpdata; */

/*   int c; /\* temporary char for storing input file characters *\/ */
/*   /\* we are looking for four strings in the output file: one at the beginning */
/*      of each data block, and one at the end. *\/ */
/*   char *str_buf = malloc(BUF_SIZE * 5); */
/*   char *s1 = NULL; */
/*   char *s2 = NULL; */
/*   char *s3 = NULL; */

/*   printf("      parsing the tmp file (%s) ..", get_loctime(ltime)); */

/*   if (inp -> md -> so_enrg == 0) { */
/*     /\* read spin-orbit data *\/ */
/*     s1 = "Eigenvalues of complex Hamiltonian"; */
/*     s2 = "Dipole transition strengths (SO states)"; */
/*     s3 = "Quadrupole transition strengths (SO states)"; */
/*   } */
/*   else if(inp -> md -> so_enrg == 1) { */
/*     /\* read spin-orbit free data *\/ */
/*     s1 = "SPIN-FREE ENERGIES"; */
/*     s2 = "Dipole transition strengths"; */
/*     s3 = "Quadrupole transition strengths"; */
/*   } */

/*   /\* create a pointer to the three data block beginners s1,s3, *\/ */
/*   const char *lookup_str[3] = {s1,s2,s3}; */

/*   stat(fn_tmpdata, &st); */
/*   sz_tmp = (int)st.st_size; /\* file size in bytes *\/ */
/*   sz_tmp = sz_tmp / 34; /\* shortest line = 34 chars *\/ */

/*   /\* open the input file *\/ */
/*   if((fp_tmpdata = fopen(fn_tmpdata, "r")) == NULL) { */
/*     fprintf(stderr, "\n\nscttr_io.c, function parse_input_tmp_el_el: unable to open the input file %s.\n" */
/*             , fn_tmpdata); */
/*     printf("program terminating due to the previous error.\n\n"); */
/*     exit(1); */
/*   } */

/*   k = 0; /\* index for tmp_string *\/ */
/*   l = 0; /\* index for lookup string *\/ */

/*   /\* storage for the energy eigenvalues *\/ */
/*   if((e_eigval = malloc(sz_tmp * sizeof(double))) == NULL ) { */
/*     fprintf(stderr, "\n\nscttr_io.c, function parse_input_tmp_el: failed to allocate memory for \"e_eigval\"\n"); */
/*     printf("program terminating due to the previous error.\n\n"); */
/*     exit(1); */
/*   } */

/*   if((trs_types = malloc(sz_tmp * sizeof(double))) == NULL ) { */
/*     fprintf(stderr, "\n\nscttr_io.c, function parse_input_tmp_el: failed to allocate memory for \"trs_types\"\n"); */
/*     printf("program terminating due to the previous error.\n\n"); */
/*     exit(1); */
/*   } */

/*   if((idxs_eigval = malloc(sz_tmp * sizeof(int))) == NULL ) { */
/*     fprintf(stderr, "\n\nscttr_io.c, function parse_input_tmp_el: failed to allocate memory for \"idxs_eigval\"\n"); */
/*     printf("program terminating due to the previous error.\n\n"); */
/*     exit(1); */
/*   } */

/*   for (j = 0; j < sz_tmp; j++) { */
/*     idxs_eigval[j] = -1; */
/*   } */

/*   /\* storage for the transition moments *\/ */
/*   if((t_mom = malloc(sz_tmp * sizeof(double))) == NULL ) { */
/*     fprintf(stderr, "\n\nscttr_io.c, function parse_input_tmp_el: failed to allocate memory for \"t_mom\"\n"); */
/*     printf("program terminating due to the previous error.\n\n"); */
/*     exit(1); */
/*   } */

/*   /\* storage for the transition indexes, column 1 is from a state */
/*      column 2 is to state index *\/ */
/*   if((trans_idxs = malloc(2 * sizeof(double *))) == NULL ) { */
/*     fprintf(stderr, "\n\nscttr_io.c, function parse_input_tmp_el: failed to allocate memory for \"trans_idxs\"\n"); */
/*     printf("program terminating due to the previous error.\n\n"); */
/*     exit(1); */
/*   } */

/*   for (j = 0; j < 2; j++) { */
/*     if((trans_idxs[j] = malloc(sz_tmp * sizeof(double))) == NULL ) { */
/*       fprintf(stderr, "\n\nscttr_io.c, function parse_input_tmp_el: failed to allocate memory for \"trans_idxs[%d]\"\n" */
/*               , j); */
/*       printf("program terminating due to the previous error.\n\n"); */
/*       exit(1); */
/*     } */
/*   } */

/*   if((num_idxs1 = malloc(2 * sizeof(int))) == NULL ) { */
/*     fprintf(stderr, "\n\nscttr_io.c, function parse_input_tmp_el: failed to allocate memory for \"num_idxs1\"\n"); */
/*     printf("program terminating due to the previous error.\n\n"); */
/*     exit(1); */
/*   } */

/*   if((num_idxs2 = malloc(3 * sizeof(int))) == NULL ) { */
/*     fprintf(stderr, "\n\nscttr_io.c, function parse_input_tmp_el: failed to allocate memory for \"num_idxs2\"\n"); */
/*     printf("program terminating due to the previous error.\n\n"); */
/*     exit(1); */
/*   } */

/*   num_idxs1[0] = 0; */
/*   num_idxs1[1] = 1; */

/*   num_idxs2[0] = 0; */
/*   num_idxs2[1] = 1; */
/*   num_idxs2[2] = 2; */

/*   j = j_test = l = m = 0; */

/*   n_states = 0; */
/*   n_trans = 0; */
/*   trs_type = 1; */

/*   /\* now that data structures of the rights size has memory allocated for */
/*      them, start reading data from the temporary file *\/ */
/*   while ((c = fgetc(fp_tmpdata)) != EOF) { */

/*     str_buf[l] = (char)c; */

/*     if ((str_buf[l]  == '\n') && (l > 0)) { /\* dont send blank lines *\/ */
/*       str_buf[l + 1] = '\0'; */
/*       if ((j_test = strscmp(str_buf, lookup_str, 3)) != -1) { */
/*         j = j_test; */
/*         trs_type = j; */
/*       } */
/*       else { */
/*         if ((j == 0) && (isempty(str_buf, l) != 1)) { */

/*           /\* extract energy eigenvalues and state indexes *\/ */
/*           get_nums(str_buf, num_idxs1, l, n_idxs1, &tmp_idx */
/*                    , &e_eigval[n_states]); */
/*           if ((fabs(e_eigval[n_states] - e_eigval[0]) * AUTOEV) \ */
/*               < maxr) { */
/*             idxs_eigval[(int)tmp_idx - 1] = n_states + 1; */
/*             n_states++; */
/*           } */
/*         } */
/*         else if ((j > 0) && (isempty(str_buf, l) != 1)) { */
/*           if (n_trans == 0) { /\* pre-process the eigenvalues *\/ */

/*             /\* sort the states in energy *\/ */
/*             iquicks(e_eigval, idxs_eigval, 0, n_states - 1, n_states); */
/*             inp -> e0 = e_eigval[0]; */

/*             /\* adjust n_states so that it accounts for states not to be read */
/*              due to being outside of the input range *\/ */
/*             m = 0; */
/*             for (k = 0; k < n_states; k++) { */
/*               if ((e_eigval[k] - inp -> e0) * AUTOEV < maxr) { */
/*                 m++; */
/*               } */
/*             } */
/*             n_states = m; */

/*           } */
/*           /\* extract transition moments and transition indexes *\/ */
/*           get_nums(str_buf, num_idxs2, l, n_idxs2, &trans_idxs[0][n_trans], */
/*                    &trans_idxs[1][n_trans], &t_mom[n_trans]); */
/*           trs_types[n_trans] = trs_type; */

/*           from_state_en = get_wi(e_eigval, idxs_eigval */
/*                                  , (int)(trans_idxs[1][n_trans]), n_states); */

/*           to_state_en = get_wi(e_eigval,idxs_eigval */
/*                                , (int)(trans_idxs[0][n_trans]), n_states); */

/*           if ((((int)to_state_en != -1) && ((int)from_state_en != -1)) */
/*               && (((to_state_en-inp -> e0) * AUTOEV < maxr) */
/*                   && ((from_state_en-inp -> e0) * AUTOEV < maxr))){ */
/*             n_trans++; */
/*           } */
/*         } */
/*       } */
/*       /\* reset the buffer write head to start reading a the next line *\/ */
/*       l = 0; */
/*     } */
/*     else{ */
/*       l++; */
/*     } */
/*   } */

/*   /\* allocate space for the "parsed input matrix" that will be filled with data */
/*      in the remaining sections of this function *\/ */

/*   if((proc_idx = malloc(n_states * sizeof(int))) == NULL ) { */
/*     fprintf(stderr, "\n\nscttr_io.c, function parse_input_tmp_el: failed to allocate memory for \"proc_idx\"\n"); */
/*     printf("program terminating due to the previous error.\n\n"); */
/*     exit(1); */
/*   } */

/*   if((inp -> trs = malloc(6 * sizeof(double *))) == NULL ) { */
/*     fprintf(stderr, "\n\nscttr_io.c, function parse_input_tmp_el: failed to allocate memory for \"inp -> trs\"\n"); */
/*     printf("program terminating due to the previous error.\n\n"); */
/*     exit(1); */
/*   } */

/*   for (j = 0; j < 6; j++) { */
/*     if((inp -> trs[j] = malloc((n_trans + 1) * sizeof(double))) == NULL ) { */
/*       fprintf(stderr, "\n\nscttr_io.c, function parse_input_tmp_el: failed to allocate memory for pointers in \"inp -> trs[%d]\"\n" */
/*               , j); */
/*       printf("program terminating due to the previous error.\n\n"); */
/*       exit(1); */
/*     } */
/*   } */
/*   inp -> trs[0][n_trans] = -1; /\* end of list *\/ */

/*   if((trs_buf = malloc(6 * sizeof(double *))) == NULL ) { */
/*     fprintf(stderr, "\n\nscttr_io.c, function parse_input_tmp_el: failed to allocate memory for \"trs_buf\"\n"); */
/*     printf("program terminating due to the previous error.\n\n"); */
/*     exit(1); */
/*   } */

/*   for (j = 0; j < 6; j++) { */
/*     if((trs_buf[j] = malloc((n_trans + 1) * sizeof(double))) == NULL ) { */
/*       fprintf(stderr, "\n\nscttr_io.c, function parse_input_tmp_el: failed to allocate memory for pointers in \"trs_buf[%d]\"\n" */
/*               , j); */
/*       printf("program terminating due to the previous error.\n\n"); */
/*       exit(1); */
/*     } */
/*   } */

/*   /\* finally, store the data in the trs matrix for the parse_input function  *\/ */
/*   n_proc = 0; */
/*   last_i = 1; */
/*   j = k = l = 0; */

/*   /\* l = index of line in the data arrays *\/ */
/*   /\* j = index of data read into trs_buf, reset upon reading the transitions */
/*      from a new state *\/ */

/*   while (l < n_trans - 1) { */
/*     /\* printf("THE DATA: to %d from %d\n ", (int)(trans_idxs[0][j+l]),(int)(trans_idxs[1][j+l])); *\/ */
/*     if ((j + l) == (n_trans)) { */
/*       /\* dont read beyond the last value in PI *\/ */
/*       trs_buf[0][j] = -1; */
/*       idx_from = -1; */
/*     } */
/*     else { */
/*       from_state_en = get_wi(e_eigval, idxs_eigval */
/*                              , (int)(trans_idxs[0][j + l]), n_states); */
/*       to_state_en = get_wi(e_eigval, idxs_eigval */
/*                            , (int)(trans_idxs[1][j + l]), n_states); */

/*       if (inrange((to_state_en - inp -> e0) * AUTOEV, state_er[3] */
/*                   , state_er[4]) */
/*           && (inrange((from_state_en - inp -> e0)*AUTOEV,state_er[5] */
/*                       ,state_er[6])) */
/*           && (((int)state_er[1] != (int)state_er[5]) */
/*               && ((int)state_er[2] != (int)state_er[6])) */
/*           ) { */

/*         idx_from = trans_idxs[1][j + l]; */
/*         idx_to = trans_idxs[0][j + l]; */

/*         tmp_state_en = from_state_en; */
/*         from_state_en = to_state_en; */
/*         to_state_en = tmp_state_en; */
/*       } */
/*       else { */

/*         idx_from = trans_idxs[0][j + l]; */
/*         idx_to   = trans_idxs[1][j + l]; */
/*       } */

/*       trs_buf[0][j] = idx_from; */
/*       trs_buf[1][j] = idx_to; */
/*       trs_buf[2][j] = from_state_en; */
/*       trs_buf[3][j] = to_state_en; */
/*       trs_buf[4][j] = t_mom[j + l]; */
/*       trs_buf[5][j] = trs_types[j + l]; */
/*     } */
/*     if (idx_from != last_i) { */

/*       tmp_idx2 = get_inext(inp -> trs, last_i); */
/*       /\* we have read all transitions for a state *\/ */
/*       /\* check if the last_i has already been processed *\/ */
/*       if (intinint(proc_idx, last_i, n_proc) == -1) { */

/*         m = l; */
/*         while ((int)trs_buf[0][m - l] != idx_from) { */

/*           inp -> trs[0][m] = trs_buf[0][m - l]; */
/*           inp -> trs[1][m] = trs_buf[1][m - l]; */
/*           inp -> trs[2][m] = trs_buf[2][m - l]; */
/*           inp -> trs[3][m] = trs_buf[3][m - l]; */
/*           inp -> trs[4][m] = trs_buf[4][m - l]; */
/*           inp -> trs[5][m] = trs_buf[5][m - l]; */

/*           m++; */
/*         } */

/*         inp -> trs[0][m] = -1; */
/*         proc_idx[n_proc++] = last_i; */
/*       } */
/*       else{ */

/*         tmp_idx2 = get_inext(inp -> trs, last_i); */
/*         fwdsplice(trs_buf, inp -> trs, tmp_idx2, l, j, 6); */
/*         fflush(stdout); */
/*       } */
/*       l += j; */
/*       j = 0; */

/*     } else { */
/*       j++; */
/*     } */
/*     last_i = idx_from; */
/*   } */

/*   inp -> trs[0][l] = -1; */
/*   inp -> n_trans = l; */

/*   if (fclose(fp_tmpdata) != 0) { */
/*     fprintf(stderr, "\n\nscttr_io.c, function parse_input_tmp_el: unable to close file:\n%s\n" */
/*             , fn_tmpdata); */
/*     printf("program terminating due to the previous error.\n\n"); */
/*     exit(1); */
/*   } */

/*   if (remove(fn_tmpdata) != 0) { */
/*     fprintf(stderr, "\n\nscttr_io.c, function parse_input_tmp_el: unable to remove temporary file:\n%s\n" */
/*             , fn_tmpdata); */
/*     printf("program terminating due to the previous error.\n\n"); */
/*     exit(1); */
/*   } */

/*   free(e_eigval); */
/*   free(trs_types); */
/*   free(idxs_eigval); */
/*   free(t_mom); */

/*   for (j = 0; j < 2; j++) { */
/*     free(trans_idxs[j]); */
/*   } */
/*   free(trans_idxs); */

/*   free(num_idxs1); */
/*   free(num_idxs2); */
/*   free(proc_idx); */

/*   for (j = 0; j<6; j++) { */
/*     free(trs_buf[j]); */
/*   } */
/*   free(trs_buf); */
/*   free(str_buf); */

/*   printf(" done (%s).\n", get_loctime(ltime)); */
/*   free(ltime); */
/*   return 1; */
/* } */

int
parse_input_tmp_el (struct inp_node *inp, char *fn_tmpdata)
{
  int j; /* the write head for the trs_buf */
  int k, l, m, j_test; /* control loop variables */
  int idx_to, idx_from; /* index in each transition from state x to y  */
  int tmp_idx2;
  int tmp_i = 0;
  int first = 0; /*  == 1 means dipole transitions have been read */
  int proci = 0;
  int n_states, n_trans;
  int last_i;
  int n_proc; /* number of states processed at any given time in the
                 execution of the function */
  int nsc = 0; /* number of screened out ground states */
  int trs_type; /* type  */

  int n_idxs1 = 2;
  int n_idxs2 = 3;
  int rc;
  int idx_comp;/* what transition index to compare to */

  /* rh (readhead): the place from which data is being read in the data arrays
     defined above */
  int READHEAD;

  /* wh (writehead): the place to where we are writing data in the trs matrix */
  int WRITEHEAD;

  char *ltime = malloc(20);

  int *num_idxs1;
  int *num_idxs2;
  int *trs_types;
  int *idxs_eigval;
  int *proc_idx;

  float bw_thrsh = inp -> md -> state_t[2];

  double tmp_bw;

  double *state_er = inp -> md -> state_er;
  double tmp_idx = 0;
  double maxr = get_maxl(state_er, state_er[0]);
  double to_state_en, from_state_en, tmp_state_en;

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

  printf("      parsing the tmp file (%s) ..", get_loctime(ltime));

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

  inp -> bw_sum = 0;

  /* now that data structures of the rights size has memory allocated for
     them, start reading data from the temporary file */
  last_i = 0;
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
          if ((fabs(e_eigval[n_states] - e_eigval[0]) * (double)AUTOEV) \
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

          to_state_en = get_wi(e_eigval, idxs_eigval
                               , (int)(trans_idxs[0][n_trans]), n_states);

          from_state_en = get_wi(e_eigval,idxs_eigval
                                 , (int)(trans_idxs[1][n_trans]), n_states);

          if ((((int)from_state_en != -1) && ((int)to_state_en != -1))
              && (((from_state_en - inp -> e0) * AUTOEV < maxr)
                  && ((to_state_en - inp -> e0) * AUTOEV < maxr))){


            /* make sure that we only count the ground state boltzmann weights
             once and in the right energy range */
            if ((last_i != (int)(trans_idxs[0][n_trans]))
                && (inrange((to_state_en - inp -> e0) * AUTOEV, state_er[1]
                            , state_er[2])
                    && inrange((from_state_en - inp -> e0)*AUTOEV,state_er[3]
                               ,state_er[4]))) {
              inp -> bw_sum += get_boltzw((to_state_en - inp -> e0) * (double)AUTOEV);
              last_i = trans_idxs[0][n_trans];
            }
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
  last_i = trans_idxs[0][0];
  j = k = l = 0;

  /* l = index of line in the data arrays */
  /* j = index of data read into trs_buf, reset upon reading the transitions
     from a new state */

  READHEAD = WRITEHEAD = 0;
  while (READHEAD <= n_trans) {

    from_state_en = get_wi(e_eigval, idxs_eigval
                         , (int)(trans_idxs[0][READHEAD]), n_states);
    to_state_en = get_wi(e_eigval, idxs_eigval
                           , (int)(trans_idxs[1][READHEAD]), n_states);

    if (inrange((from_state_en - inp -> e0) * AUTOEV, state_er[1]
                , state_er[2])
        && inrange((to_state_en - inp -> e0)*AUTOEV,state_er[3]
                   ,state_er[4]))
      {         /* triggered on any transition */
        idx_to = trans_idxs[1][READHEAD];
        idx_from = trans_idxs[0][READHEAD];

        proci = 1;
        idx_comp = idx_to;
        first = 0;
      }
    else {
      proci = 0;
    }


    if (READHEAD == n_trans) {
      /* Make sure also the last set of transitions are added. */
      /* dont read beyond the last value in PI and mark the end of the array */
      READHEAD++;
      proci = 1;
      idx_comp = -last_i;
    }

    if (proci) { /* a transition was found that fit the energy ranges */

      if (idx_comp != last_i) {
        /* all transitions for a specific "to" state have been found */

        tmp_idx2 = get_inext(inp -> trs, last_i);

        /* we have read all transitions for a state */
        /* check if the last_i has already been processed or if it is
           the last state in the list */
        if ((intinint(proc_idx, last_i, n_proc) == -1)
            || (tmp_idx2 == WRITEHEAD - 1)
            || (tmp_idx2 == -1)) {

          m = 0;
          for (m = 0; m < j ; m++, WRITEHEAD++) {

            inp -> trs[0][WRITEHEAD] = trs_buf[0][m];
            inp -> trs[1][WRITEHEAD] = trs_buf[1][m];
            inp -> trs[2][WRITEHEAD] = trs_buf[2][m];
            inp -> trs[3][WRITEHEAD] = trs_buf[3][m];
            inp -> trs[4][WRITEHEAD] = trs_buf[4][m];
            inp -> trs[5][WRITEHEAD] = trs_buf[5][m];
          }

          inp -> trs[0][WRITEHEAD] = -1;
          proc_idx[n_proc++] = last_i;

        }
        else {
          /* The to state transition has to be spliced into the matrix since
             it has already had transitions from it registered previously  */
          tmp_idx2 = get_inext(inp -> trs, last_i);
          rc = fwdsplice(trs_buf, inp -> trs, tmp_idx2, WRITEHEAD, j, 6);

          if (rc == 1) {

            fflush(stdout);
            fprintf(stderr, "\n\n\n Error! scttr_io.c, function parse_input_tmp: the fwdsplice function return a non-zero integer. Printing debug information\n");
            printf("input arguments: %d, %d, %d\n", tmp_idx2, WRITEHEAD, j);
            printf("additional variables: last_i = %d, idx_to = %d processed = %d\n", last_i, idx_to, intinint(proc_idx, last_i, n_proc));
            printf("Reading transition %d/%d\n", READHEAD, n_trans);
            printf( "program terminating due to the previous error.\n");
            printf(" trs matrix dump: \n" );
            for (l = 0; l < WRITEHEAD; l++) {
              printf("%le %le %le %le %le\n", inp -> trs[0][l],inp -> trs[1][l], inp -> trs[2][l],inp -> trs[3][l],inp -> trs[4][l]);
            }
            for (l = 0; l < n_proc; l++) {
              printf("proc_idx = %d\n", proc_idx[l]);
            }
            fflush(stdout);
            exit(1);
          }

          WRITEHEAD += j;
        }

        j = 0;
        last_i = idx_comp;
        READHEAD--;
      }
      else {

        trs_buf[0][j] = idx_to;
        trs_buf[1][j] = idx_from;
        trs_buf[2][j] = to_state_en;
        trs_buf[3][j] = from_state_en;
        trs_buf[4][j] = t_mom[READHEAD];
        trs_buf[5][j] = trs_types[READHEAD];

        j++;

      }
    }
    READHEAD++;
  }

  inp -> n_trans = WRITEHEAD;

  if (fclose(fp_tmpdata) != 0) {
    fprintf(stderr, "\n\nscttr_io.c, function parse_input_tmp: unable to close file:\n%s\n"
            , fn_tmpdata);
    printf("program terminating due to the previous error.\n\n");
    exit(1);
  }

  if (remove(fn_tmpdata) != 0) {
    fprintf(stderr, "\n\nscttr_io.c, function parse_input_tmp: unable to remove temporary file:\n%s\n"
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

  printf(" done (%s).\n", get_loctime(ltime));
  free(ltime);

  return 1;
}

int
parse_input_tmp (struct inp_node *inp, char *fn_tmpdata)
{

  int j; /* the write head for the trs_buf */
  int k, l, m, j_test; /* control loop variables */
  int idx_to, idx_from; /* index in each transition from state x to y  */
  int tmp_idx2;
  int tmp_i = 0;
  int first = 0; /*  == 1 means dipole transitions have been read */
  int proci = 0;
  int n_states, n_trans;
  int last_i;
  int n_proc; /* number of states processed at any given time in the
                 execution of the function */
  int nsc = 0; /* number of screened out ground states */
  int trs_type; /* type  */

  int n_idxs1 = 2;
  int n_idxs2 = 3;
  int rc;
  int idx_comp;/* what transition index to compare to */

  /* rh (readhead): the place from which data is being read in the data arrays
     defined above */
  int READHEAD;

  /* wh (writehead): the place to where we are writing data in the trs matrix */
  int WRITEHEAD;

  char *ltime = malloc(20);

  int *num_idxs1;
  int *num_idxs2;
  int *trs_types;
  int *idxs_eigval;
  int *proc_idx;

  float bw_thrsh = inp -> md -> state_t[2];

  double tmp_bw;

  double *state_er = inp -> md -> state_er;
  double tmp_idx = 0;
  double maxr = get_maxl(state_er, state_er[0]);
  double to_state_en, from_state_en, tmp_state_en;

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

  printf("      parsing the tmp file (%s) ..", get_loctime(ltime));

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

  inp -> bw_sum = 0;

  /* now that data structures of the rights size has memory allocated for
     them, start reading data from the temporary file */
  last_i = 0;
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
          if ((fabs(e_eigval[n_states] - e_eigval[0]) * (double)AUTOEV) \
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

          to_state_en = get_wi(e_eigval, idxs_eigval
                               , (int)(trans_idxs[0][n_trans]), n_states);

          from_state_en = get_wi(e_eigval,idxs_eigval
                                 , (int)(trans_idxs[1][n_trans]), n_states);

          if ((((int)from_state_en != -1) && ((int)to_state_en != -1))
              && (((from_state_en - inp -> e0) * AUTOEV < maxr)
                  && ((to_state_en - inp -> e0) * AUTOEV < maxr))){


            /* make sure that we only count the ground state boltzmann weights
             once and in the right energy range */
            if ((last_i != (int)(trans_idxs[0][n_trans]))
                && (inrange((to_state_en - inp -> e0) * AUTOEV, state_er[1]
                            , state_er[2])
                    && inrange((from_state_en - inp -> e0)*AUTOEV,state_er[3]
                               ,state_er[4])
                    && (trs_type == inp -> abs))) {
              inp -> bw_sum += get_boltzw((to_state_en - inp -> e0) * (double)AUTOEV);
              last_i = trans_idxs[0][n_trans];
            }
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
  last_i = trans_idxs[0][0];
  j = k = l = 0;

  /* l = index of line in the data arrays */
  /* j = index of data read into trs_buf, reset upon reading the transitions
     from a new state */
  /* printf("NEGHSSHRH %d\n", n_trans); */
  READHEAD = WRITEHEAD = 0;
  while (READHEAD <= n_trans) {

    to_state_en = get_wi(e_eigval, idxs_eigval
                         , (int)(trans_idxs[0][READHEAD]), n_states);
    from_state_en = get_wi(e_eigval, idxs_eigval
                           , (int)(trans_idxs[1][READHEAD]), n_states);

    /* printf("energies: %le %le \n", (to_state_en - inp -> e0) * AUTOEV, (from_state_en - inp -> e0) * AUTOEV); */

    if (inrange((to_state_en - inp -> e0) * AUTOEV, state_er[1]
                , state_er[2])
        && ((inrange((from_state_en - inp -> e0)*AUTOEV,state_er[3]
                     ,state_er[4]))
            || (inrange((from_state_en - inp -> e0)*AUTOEV,state_er[1]
                        ,state_er[2])))
        && (trs_types[READHEAD] == inp -> abs))
      {         /* triggered on absorption transitions */

        /* fprintf(stderr, "\n\n=======foundquadValgrind eject point=======\n\n"); */
        /* exit(1); */
        idx_to = trans_idxs[1][READHEAD];
        idx_from = trans_idxs[0][READHEAD];

        tmp_state_en = to_state_en;
        to_state_en = from_state_en;
        from_state_en = tmp_state_en;
        proci = 1;
        idx_comp = idx_to;
        first = 0;
      }
    else if (inrange((to_state_en - inp -> e0) * AUTOEV, state_er[5]
                     , state_er[6])
             && ((inrange((from_state_en - inp -> e0)*AUTOEV,state_er[3]
                          ,state_er[4]))
                 || (inrange((from_state_en - inp -> e0)*AUTOEV,state_er[1]
                             ,state_er[2])))
             /* && (inrange((from_state_en - inp -> e0)*AUTOEV,state_er[3] */
             /*             ,state_er[4])) */
             && (trs_types[READHEAD] == inp -> ems))
      { /* triggered on emission transitions */
        /* fprintf(stderr, "\n\n=======founddipoleValgrind eject point=======\n\n"); */
        /* exit(1); */
        idx_to = trans_idxs[0][READHEAD];
        idx_from = trans_idxs[1][READHEAD];
        first = 1;
        proci = 1;
        idx_comp = idx_to;
      }
    else {
      proci = 0;
    }

    /* fflush(stdout); */
    /* /\* printf("\nNUM NUM %d %d %d\n", trs_types[READHEAD], trs_types[READHEAD] == inp -> ems, trs_types[READHEAD] == inp -> abs ); *\/ */
    /* fflush(stdout); */

    if (READHEAD == n_trans) {
      /* Make sure also the last set of transitions are added. */
      /* dont read beyond the last value in PI and mark the end of the array */
      READHEAD++;
      proci = 1;
      idx_comp = -last_i;
    }

    if (proci) { /* a transition was found that fit the energy ranges */

      if (idx_comp != last_i) {
        /* all transitions for a specific "to" state have been found */

        tmp_idx2 = get_inext(inp -> trs, last_i);

        /* we have read all transitions for a state */
        /* check if the last_i has already been processed or if it is
           the last state in the list */
        if ((intinint(proc_idx, last_i, n_proc) == -1)
            || (tmp_idx2 == WRITEHEAD - 1)
            || (tmp_idx2 == -1)) {

          m = 0;
          for (m = 0; m < j ; m++, WRITEHEAD++) {

            inp -> trs[0][WRITEHEAD] = trs_buf[0][m];
            inp -> trs[1][WRITEHEAD] = trs_buf[1][m];
            inp -> trs[2][WRITEHEAD] = trs_buf[2][m];
            inp -> trs[3][WRITEHEAD] = trs_buf[3][m];
            inp -> trs[4][WRITEHEAD] = trs_buf[4][m];
            inp -> trs[5][WRITEHEAD] = trs_buf[5][m];
          }

          inp -> trs[0][WRITEHEAD] = -1;
          proc_idx[n_proc++] = last_i;

        }
        else {
          /* The to state transition has to be spliced into the matrix since
             it has already had transitions from it registered previously  */
          tmp_idx2 = get_inext(inp -> trs, last_i);
          rc = fwdsplice(trs_buf, inp -> trs, tmp_idx2, WRITEHEAD, j, 6);

          if (rc == 1) {

            fflush(stdout);
            fprintf(stderr, "\n\n\n Error! scttr_io.c, function parse_input_tmp: the fwdsplice function return a non-zero integer. Printing debug information\n");
            printf("input arguments: %d, %d, %d\n", tmp_idx2, WRITEHEAD, j);
            printf("additional variables: last_i = %d, idx_to = %d processed = %d\n", last_i, idx_to, intinint(proc_idx, last_i, n_proc));
            printf("Reading transition %d/%d\n", READHEAD, n_trans);
            printf( "program terminating due to the previous error.\n");
            printf(" trs matrix dump: \n" );
            for (l = 0; l < WRITEHEAD; l++) {
              printf("%le %le %le %le %le\n", inp -> trs[0][l],inp -> trs[1][l], inp -> trs[2][l],inp -> trs[3][l],inp -> trs[4][l]);
            }
            for (l = 0; l < n_proc; l++) {
              printf("proc_idx = %d\n", proc_idx[l]);
            }
            fflush(stdout);
            exit(1);
          }

          WRITEHEAD += j;
        }

        j = 0;
        last_i = idx_comp;
        READHEAD--;
      }
      else {

        trs_buf[0][j] = idx_to;
        trs_buf[1][j] = idx_from;
        trs_buf[2][j] = to_state_en;
        trs_buf[3][j] = from_state_en;
        trs_buf[4][j] = t_mom[READHEAD];
        trs_buf[5][j] = trs_types[READHEAD];

        j++;

      }
    }
    READHEAD++;
  }

  inp -> n_trans = WRITEHEAD;

  if (fclose(fp_tmpdata) != 0) {
    fprintf(stderr, "\n\nscttr_io.c, function parse_input_tmp: unable to close file:\n%s\n"
            , fn_tmpdata);
    printf("program terminating due to the previous error.\n\n");
    exit(1);
  }

  if (remove(fn_tmpdata) != 0) {
    fprintf(stderr, "\n\nscttr_io.c, function parse_input_tmp: unable to remove temporary file:\n%s\n"
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

  printf(" done (%s).\n", get_loctime(ltime));
  free(ltime);

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
  char *bin_fpstr;
  char *tmp_fpstr = concs(3,md -> outpath,md -> inp_fn, tmp_sfx);
  char *format = md -> inp_sfx;
  char *inpath = md -> inpath;

  /* did the user provide a binary file as input? */
  if (strcmp(format, bin_sfx) <= 0) {
    bin_fpstr = concs(1, md -> inpath);
  } else {
    bin_fpstr =  concs(3, md -> outpath, md -> inp_fn, bin_sfx);
  }

  fp_bin = fopen(bin_fpstr, "r");
  /* if (fp_bin != NULL) { */
  if (NULL) {
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
      }
      else {
        if (inp -> el == 1) {
          parse_input_tmp_el(inp, tmp_fpstr);
        }
        else if(inp -> el == 0) {
          parse_input_tmp(inp, tmp_fpstr);
        }
      }
    }

    /* if (inp -> el == 1) { */
      printf("\nPRE adding\n ");
      for (j = 0; j < inp->n_trans; j++) {
        printf("%d %d %le %le %le\n", (int)inp -> trs[0][j], (int)inp -> trs[1][j], (inp -> trs[2][j] - inp -> e0) * AUTOEV, (inp -> trs[3][j] - inp -> e0) * AUTOEV, inp -> trs[4][j]);
      }
      /* count_states(inp); */
      add_eltrans(inp);
      fflush(stdout);
      printf("\nPOST adding\n ");
      printf("%d, %d\n", inp->n_gfs, inp->n_is );

      for (j = 0; j < inp->n_trans; j++) {
        printf("%d %d %le %le %le\n", (int)inp -> trs[0][j], (int)inp -> trs[1][j], (inp -> trs[2][j] - inp -> e0) * AUTOEV, (inp -> trs[3][j] - inp -> e0) * AUTOEV, inp -> trs[4][j]);
      }
      /* fprintf(stderr, "\n\n=======Valgrind eject point=======\n\n"); */
      /* exit(1); */


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
    fwrite((const void*) & (inp -> bw_sum), sizeof(double), 1, fp_bin);
    fflush(fp_bin);
    fclose(fp_bin);

    if ((fp_bin = fopen(bin_fpstr, "ab")) == NULL) {
      fprintf(stderr, "\n\nscttr_io.c, function parse_input: unable to open the binary file used to store the trs matrix: %s\n"
              , bin_fpstr);
      printf("program terminating due to the previous error.\n\n");
      exit(1);
    }

    for (j = 0; j < 6; j++) {
      if ((int)fwrite(inp -> trs[j], sizeof(double), inp -> n_trans, fp_bin) != inp -> n_trans) {
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
  /* printf("\nPOST adding\n "); */
  /* for (j = 0; (int)inp -> trs[0][j] != -1; j++) { */
  /*   printf("%d %d %le %le %le %d\n", (int)inp -> trs[0][j], (int)inp -> trs[1][j], (inp -> trs[2][j] - inp -> e0) * AUTOEV, (inp -> trs[3][j] - inp -> e0) * AUTOEV, inp -> trs[4][j], (int)inp -> trs[5][j]); */
  /* } */
  /* fprintf(stderr, "\n\n=======Valgrind eject point=======\n\n"); */
  /* exit(1); */
  /* fprintf(stderr, "\n\n=======Valgrind eject point=======\n\n"); */
  /* exit(1); */

  /* if (inp -> el == 1) { */
    /* rc = set_root_spec_el(inp); */
  /* } */
  /* else if(inp -> el == 0) { */
  rc = set_root_spec(inp);
  /* } */

  if (rc != 0) {

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

  return 0;
}

int
write_spec (struct inp_node *inp, struct spectrum *spec)
{
  int j, k, x, y;

  FILE *fp_smat_out;
  char *dat_fpstr = concs(3, inp -> md -> outpath,
                          inp -> md -> inp_fn, dat_sfx);

  double emin_x = spec -> emin_x;
  double emin_y = spec -> emin_y;
  double de_x = inp -> md -> res[0] / AUTOEV;
  double de_y = inp -> md -> res[1] / AUTOEV;

  if((fp_smat_out = fopen(dat_fpstr, "w"))==NULL) {
    fprintf(stderr, "\n\nscttr_io.c, function write_spec: unable to open the output file %s.\n"
            , dat_fpstr);
    printf("program terminating due to the previous error.\n\n");
    exit(1);
  }

  printf("  - writing RIXS map to file: %s ..", dat_fpstr);

  for (j = 0, x = 0, y = 0; x < spec -> n_elx/* j < spec -> npr_tot */; j++) {
    for (k = 0; (k < spec -> prsz) && (x < spec -> n_elx); k++, y++) {
      fprintf(fp_smat_out,"%le %le %le\n", ((emin_x + (x * de_x)))
              * AUTOEV, (emin_y + (y * de_y)) * AUTOEV,
              spec -> s_mat[j][k] / spec -> sfac);
      fflush(fp_smat_out);
      if (y == spec -> n_ely-1) {
        x++;
        y = 0;
        fprintf(fp_smat_out,"\n");
        fflush(fp_smat_out);
      }
    }
  }

  printf(" done.\n");
  if (fclose(fp_smat_out) != 0) {
    fprintf(stderr, "\n\nscttr_io.c, function write_spec:: unable to close some of the files:\n%s\n",
            dat_fpstr);
    printf("program terminating due to the previous error.\n\n");
    exit(1);
  }

  free(dat_fpstr);
  return 0;
}

int
write_plotscript (struct inp_node *inp, struct spectrum *spec)
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
          spec -> emin_x * AUTOEV,
          spec -> emax_x * AUTOEV,
          spec -> emin_y * AUTOEV,
          spec -> emax_y * AUTOEV,
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
  return 0;
}

int
write_timings (struct inp_node *inp)
{
  int nth = 0;
  env2int("OMP_NUM_THREADS", &nth);
  FILE *fp_time_out;
  char *time_fpstr = concs(3, inp -> md -> outpath,
                           inp -> md -> inp_fn, time_sfx);
  double st;

  if((fp_time_out = fopen(time_fpstr, "a"))==NULL) {
    fprintf(stderr, "\n\nscttr_io.c, function write_timings: unable to open the output file %s.\n"
            , time_fpstr);
    printf("program terminating due to the previous error.\n\n");
    exit(1);
  }

  printf("  - writing execution timing data to file: %s ..", time_fpstr);
  st = omp_get_wtime() - serial_t - para_t;
  total_t = st + para_t;
  fprintf(fp_time_out,"#sfac = %le\n", get_spec(inp,2) -> sfac);
  fprintf(fp_time_out,"Ncores\tTserial\tTparallel\tTtotal\n" );
  fprintf(fp_time_out,"%d ", nth);
  fprintf(fp_time_out,"%le ", st);
  fprintf(fp_time_out,"%le ", para_t);
  fprintf(fp_time_out,"%le\n", total_t);

  if (fclose(fp_time_out) != 0) {
    fprintf(stderr, "\n\nscttr_io.c, function write_spec:: unable to close the file:\n%s\n",
            time_fpstr);
    printf("program terminating due to the previous error.\n\n");
    exit(1);
  }

  free(time_fpstr);
  printf(" done.\n");
  return 0;
}

int
write_sticks (struct inp_node *inp, struct spectrum *spec, struct metadata *md)
{

  int j, k, l; /* iteration variables */

  char *stick_fpstr = concs(3, inp -> md -> outpath,
                          inp -> md -> inp_fn, stick_sfx);
  FILE *fp_stick_out;

  if((fp_stick_out = fopen(stick_fpstr, "w"))==NULL) {
    fprintf(stderr, "\n\nscttr_io.c, function write_sticks: unable to open the output file %s.\n"
            , stick_fpstr);
    printf("program terminating due to the previous error.\n\n");
    exit(1);
  }

  printf("  - writing stick intensities to file: %s ..", stick_fpstr);

  double tmom_gi; /* .. of excitiation (x-axis)*/
  double tmom_if; /* .. of energy transfer (y-axis)*/
  double de_gi, de_if; /* energy eigenvalue differences */
  /* double tmom_max = 0; */
  double bw; /* boltzmann weight */

  double complex tmp;  /* accumulator used in the Kramers-Heisenberg formula */

  double ** tr = spec -> trs_red;

  for (l = 0; tr[l][1] != -1;) {
    de_if = tr[l][2];
    tmom_if = tr[l][3];

    while((int)tr[++l][1] == 0) { /* loop over ground to intermediate
                                     transitions */
      bw = tr[l][0];
      de_gi = tr[l][2];
      tmom_gi = tr[l][3];

      tmp = tmom_gi * tmom_if * bw;
      tmp += 0*I;
      tmp = fabsc(tmp);

      /* if (fabsc(tmp) > tmom_max) { */
      /*   tmom_max = tmp; */
      /* } */
      /* printf("%le %le %le\n",de_gi * AUTOEV, (de_gi - de_if) * AUTOEV, creal(tmp)); */
      fprintf(fp_stick_out, "%le %le %le\n", de_gi * AUTOEV, (de_gi + de_if)
              * AUTOEV, creal(tmp));
    }
    fprintf(fp_stick_out, "\n");
  }

  /* printf("\nstick scaling factor: %le\n", tmom_max); */

  if (fclose(fp_stick_out) != 0) {
    fprintf(stderr, "\n\nscttr_io.c, function write_sticks:: unable to close the file:\n%s\n",
            stick_fpstr);
    printf("program terminating due to the previous error.\n\n");
    exit(1);
  }

  free(stick_fpstr);
  printf(" done.\n");
  return 0;
}
