## Introduction
The Scttr program is a Free and Open Source Software (FOSS) alternative to calculating resonant inelastic X-ray scattering spectra from results obtainable from standard quantum chemical software. It is written in C and parallelized over OpenMP, and further optimized for cache memory management. Its input is a simple interface that reads energy eigenvalues and transition moments from either a binary or text file, and uses that information to calculate the experimentally observed spectrum obtained from performing resonant inelastic x-ray scattering (RIXS) experiments. As of version 1.55, the program does this without taking destructive interference between scattering channels into account, as can be see in how the Kramers-Heisenberg formula is used in the program (see M. Lundberg et. al. 2013 for an example of how it is used and defined mathematically, and the calc_spec() function for the implementation). To do this, the user is provided with a well, documented command line interface giving the user a variety of options (spectral resolution, transition screening parameters etc.) to optimize the calculation time for a given input file.

## Compilation
CMake builds the makefile subsequently used to generate the binary for the program. Build options currently limited to building a "Release" type binary (the default), or a "Debug" type binary. These options are chosen through the DCMAKE_BUILD_TYPE parameter, which is provided by the user from command line:

cmake . -DCMAKE_BUILD_TYPE=production
or ..
cmake . -DCMAKE_BUILD_TYPE=debug

The default if no flag is provided is a production build.

These flags modify what flags are provided to the GCC compiler (the only compiler currently supported). The Debug option provides GCC with flags used to generate a binary suitable for debugging the program with Valgrind (see http://valgrind.org/) or the The Gnu Project Debugger (GDB, see https://www.gnu.org/software/gdb/), while the Release option compiles the program with o3 optimization (see the GCC manual for details https://gcc.gnu.org/).

After building the make files (in or out of source), use make to build the program and view the Doxygen documentation contained in the /doc directory of build root.

## Current list of publications using the program:
E. Källman, M. Guo, M. Delcey, R. Lindh, M. Lundberg, Modeling of hard X-ray processes targeting valence excited states in metal
complexes, Manuscript

E. Källman, M. Delcey, R. Lindh, M. Lundberg, Quantitive measure of spectral similarity in highly structured spectra: Examples from
soft x-ray spectroscopy, Manuscript.

M. Guo, E. Källman, R.V. Pinjari, R.C. Couto, L. Kragh Sørensen, R. Lindh, K. Pierloot, and M. Lundberg, Fingerprinting electronic
structure of heme iron by ab initio modeling of metal L-edge X-ray absorption spectra, J. Chem. Theory Comput. 2019, 15, 477−489.

M. Guo, E. Källman, L. Kragh Sørensen, M.G. Delcey, R.V. Pinjari, and M. Lundberg, Molecular orbital simulations of metal 1s2p
resonant inelastic X-ray scattering, J. Phys. Chem. A, 2016, 120, 5848−5855.

I offer the program up for public use without any claims of functionality and according to LGPL 3.0. In some ways it is my quet protest against the extremely poorly documented and revisioned code circulating in academia. Its not pretty, but its free and does its job efficiently.
