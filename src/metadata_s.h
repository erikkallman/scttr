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
   * @file metadata_s.h
   * @author Erik Källman
   * @date November 2015
   * @brief This file contains the implementation of the metadata struct.
   */
#ifndef METADATA_S_H
#define METADATA_S_H

/**
   * @brief All data provided from the user input (like file paths, energy
   * ranges, and so on), besides the input file, is stored in the metadata struct.
   * This data is then used throughout the program, mostly for i/o purposes.
   */
struct metadata
{

  int sz_inp; /**< size (in bytes) of the input file */
  int so_enrg; /**< if == 1, the program reads spin-orbit energies,
                  if not, reads spin-free*/

  char *outpath; /**< path to output directory */
  char *inpath; /**< path to input directory */

  char *inp_fn; /**< input filename */
  char *inp_sfx; /**< input file suffix */

  double *state_er; /**< user-provided ranges of energy eigenvalues in the input */
  double *state_t; /**< transition intensity screening thresholds */
  double *res; /**< spectral resolution (eV) of the produced spectrum */
  double *fwhm; /**< full-width half-maximum value for broadening the peaks */
};

#endif /* METADATA_S_H */