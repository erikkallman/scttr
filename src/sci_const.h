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
 * @file sci_const.h
 * @author Erik Källman
 * @date November 2015
 * @brief The sci_const header file is used to store physical constants that
 * are useful in scientific calculations.
 */
#ifndef CONST_H
#define CONST_H

#define AUTOJ 4.3597482e-18 /* atomic units to joule */
#define BCONST 1.380648813e-23 /* Boltzmann constant in J K^-1 */
#define ECHARGE 1.60217733e-19 /* electron charge */
#define AUTOEV 27.2113961317877 /* atomic units to eV conversion factor */
/* #define AUTOEV AUTOJ/ECHARGE /\* atomic units to eV conversion factor *\/  */
#define TEXP 298.15 /* experimental ambient temperature */
#define TTOEV 8.617385e-05 /* temperature to eV conversion factor */
#define PI 3.14159265359
#define CTOK 273.15
#endif /* CONST_H */
