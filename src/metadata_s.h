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
  int intf_mode; /**< A flag used by to set what type of interference theory is used in the spectrum calculation. */

  int sz_inp; /**< Size (in bytes) of the input file. */
  int so_enrg; /**< if == 1, the program reads spin-orbit energies,
                  if not, reads spin-free*/
  int v; /**< verbosity mode */
  int lorz; /**< a flag for applying lorentzian broadening or not */

  char *outpath; /**< Path to output directory. */
  char *inpath; /**< Path to input directory. */

  char *inp_fn; /**< The input filename */
  char *inp_sfx; /**< input file suffix */

  double *state_er; /**< User-provided ranges of energy eigenvalues in the input. */
  double *state_t; /**< Transition intensity screening thresholds. */
  double *res; /**< Spectral resolution (eV) of the produced spectrum. */

  /**< Full-width half-maximum value for broadening the peaks. This
  array contains sets of 3 numbers, each associating a broadening factor
  to a bound energy range in one dimension of the spectrum, as follows:
     ab[0] = number of broadenings to apply in this dimension,
     ab[n] = broadening factor,
     ab[n+1] = lower energy bound,
     ab[n+2] = upper energy bound,
  where a defines the type of broadening (gaussian or lorentzian) and b
  the dimension of the broadening (x or y), and n, index of a given
  broadening value (always >0).
  */
  double *gx;
  double *gy;
  double *lx;
  double *ly;
};

#endif /* METADATA_S_H */
