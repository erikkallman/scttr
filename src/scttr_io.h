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
   * @file scttr_io.h
   * @author Erik Källman
   * @date November 2015
   * @brief All functions that relate directly to the manipulation or storage
   * of user i/o has an interface defined in this header file.
   */
#ifndef SCTTR_IO_H
#define SCTTR_IO_H
#include "spectrum_s.h"
#include "inp_node_s.h"
#include "metadata_s.h"
#include "dyn_array.h"
#include "inp_node_s.h"
#include "dyn_array.h"

/* suffixes for the output and input files */
extern const char *dat_sfx;
extern const char *plot_sfx;
extern const char *log_sfx;
extern const char *time_sfx;
extern const char *bin_sfx;
extern const char *tmp_sfx;

/**
   * @brief Allocates memory for a new metadata struct.
   * @returns @p new_md the newly created metadata struct.
   */
struct metadata *
init_md ();

/**
   * @brief Frees up the memory previously allocated for a given metadata struct.
   * @param md the metadata struct to have its memory freed.
   * @returns 1 if successful.
   * @note side effects: exits with EXIT_FAILURE upon failed memory allocation.
   */
int
free_md (struct metadata *md);

/**
   * @brief Using the root input node, this function loops of the linked list
   * of input nodes until it finds the one with the right @p id.
   * @param idx defines the unique string identifying this input
   * node in the linked list of nodes.
   * @returns curr_inp a pointer to the input node with an input file name
   * corresponding to @p id.
   */
struct inp_node *
get_inp (char *idx);

/**
   * @brief Allocates memory for a new inp_node struct and appends it to the
   * list of input nodes, connected to the global root_spec node.
   * @param md the metadata used to initialize the input node.
   * @returns @p new_inp the newly allocated inp_node struct.
   * @warning Without an associated metadata struct having been defined prior
   * to the initialization of the input node, most of the functions in the
   * program will be incompatible with the input node.
   */
struct inp_node *
init_inp (struct metadata *md);

/**
   * @brief Frees up the memory previously allocated for a given input node
   * struct, while keeping the linked list structure intact.
   * @param inp the inp_node struct to have its memory freed.
   * @returns @p 1 if successful.
   * @note side effects: exits with EXIT_FAILURE upon failed memory allocation.
   */
int
free_inp (struct inp_node *inp);

/**
   * @brief Should the user provide a binary file as input, this function
   * parses that file.
   *
   * The first bit of the binary input file has to code for the number of
   * elements in each row of the @p trs matrix, since this value is used to
   * allocate the correct amount of memory to it.
   * and @p n_trans.
   *
   * @param inp
   * @returns @p 1 if successful.
   * @note side effects: defines and initializes variable of @p inp: @p n_trans, trs, e0.
   * @note If developing an interface to handle output files from program.
   * currently not supported, the bin format used in scttr is a good target.
   */
int
parse_input_bin (struct inp_node *inp, char *bin_fpstr);

/**
   * @brief The Molcas program (see http://www.molcas.org) generates .log files
   * that contain all data needed to generate a spectrum matrix. This function
   * enables the parsing of that output file format to a "tmp" format file.
   *
   * Essentially, this function simply reduces the log file and dumps only
   * the relevant information into a .tmp file that is readable by
   * parse_input_tmp(). The tmp format is more general and provides an easy
   * format to view exactly what data is read by the program. It also takes up
   * less space.
   * @note side effects: opens, writes to and then closes a tmp file with a name
   * corresponding to the input file.
   *
   * @param inp input node corresponding to the input file at @p fn_relpath.
   * @param fn_relpath path to Molcas output file (the input file to scttr).
   * @param tmp_fpstr path to where scttr writes the .tmp file.
   * @returns @p 1 if successful.
   * @note side effects: exits with EXIT_FAILURE upon failed memory allocation.
   */
int
parse_molout (struct inp_node *inp, char *fn_relpath, char *tmp_fpstr);

/**
   * @brief After an output file has be reduced to the tmp format, this function
   * reads and parses the data contained in that file to define and initialize
   * the @p trs matrix of the input file node @p inp provided as argument to
   * the function.
   *
   * @param inp the input node struct in which to store the read data.
   * @param fn_tmpdata path the to the tmp file previously generated.
   * @returns 1 if successful.
   * @note the function allocates temporary data buffers whos memory is defined
   * from the size of the input file, assuming the shortest line in the file
   * is 34 characters long. Should the program terminate due to memory
   * allocation problems inside this function, a tmp file that is too large
   * could be the reason.
   * @note side effects: exits with EXIT_FAILURE upon failed memory allocation.
   */
int
parse_input_tmp (struct inp_node *inp, char *fn_tmpdata);

int
parse_input_tmp_el (struct inp_node *inp, char *fn_tmpdata);

/**
   * @brief In addition to figuring out which function to use for parse the
   * provided input file given its suffix (.log, .tmp, or .bin), this function
   * also calls additional helper functions in order to complete the matrix of
   * transitions (see the @p trs variable of the inp_node struct)
   * with elastic transitions, should this be given from the energy
   * ranges provided by the user (@p [0,50,7000,7500,0,50] , for example).
   *
   * When the input has been read, and the @p trs variable defined and
   * initialized, this function also creates a map of what transitions falls
   * into what interval of the energy range provided by the user (see
   * set_root_spec() for details. Furthermore, a binary file, containing
   * the data in the @p trs variable, is written to the @p bin_fpstr path
   * (see parse_input_bin() for details on how this file is used in the pgram).
   *
   * @param inp the input node in which the @p trs variable and the root
   * spectrum node will be stored.
   * @returns @p 1 if successful.
   * @note side effects: exits with EXIT_FAILURE upon failed memory allocation.
   */
int
parse_input (struct inp_node *inp);

/**
   * @brief Provided with a given spectrum matrix (@p s_mat) found in
   * @p spec, this function writes that spectrum to an output file formated
   * such that it can be easily read and plotted in program such as Gnuplot
   * (see http://www.gnuplot.info/).
   *
   * @param inp the input node containing the strings required for the function
   * to write the spectrum to the correct output path.
   * @param spec the spectrum struct whos s_mat data is to be written to file.
   * @returns @p 1 if successful.
   * @note side effects: exits with EXIT_FAILURE upon failed write.
   */

int
write_spec (struct inp_node *inp,
            struct spectrum *spec);

/**
   * @brief This function reads the plot_template file in the src directory of
   * the program and writes it to the output path provided by the user,
   * appended with information derived from the spectrum struct provided as
   * argument to the function.
   *
   * The format of the plotting script is readable with Gnuplot
   * (see http://www.gnuplot.info/).
   *
   * @param inp the input node containing the strings required for the function
   * to write the plot script to the correct output path.
   * @param spec the spectrum struct containing data needed by the program.
   * to generate the right energy ranges in the plot, etc.
   * @returns @p 1 if successful.
   * @note side effects: exits with EXIT_FAILURE upon failed write.
   */
int
write_plotscript (struct inp_node *inp,
                  struct spectrum *spec);

int
write_timings (struct inp_node *inp);

#endif /* SCTTR_IO_H */
