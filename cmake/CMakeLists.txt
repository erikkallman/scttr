# This file is part of scttr.

# scttr is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# scttr is distributed in the hope that it will be useful,
# but without any warranty; without even the implied warranty of
# merchantability or fitness for a particular purpose. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with scttr, found in the "license" subdirectory of the root
# directory of the scttr program. If not, see <http://www.gnu.org/licenses/>.

# If the user specifies -DCMAKE_BUILD_TYPE on the command line, take its
# definition and dump it in the cache along with proper documentation,
# otherwise set CMAKE_BUILD_TYPE to Debug prior to calling PROJECT()
#

cmake_minimum_required(VERSION 2.8)
project(scttr C)
enable_language(C)
set (SCTTR_VERSION_MAJOR 1)
set (SCTTR_VERSION_MINOR 5)

if ( ${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR} )
  message( FATAL_ERROR "As CMake is currently configured in-source builds not allowed. In order to build Scttr, create and/or enter the directory where you want the program to be installed, then run CMake from there. You may need to remove the CMakeCache.txt file if CMake has been ran in the directory previously without being deleted." )
endif()

IF(DEFINED CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose the type of
build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug
Release RelWithDebInfo MinSizeRel.")
ELSE()
  SET(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build,
options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release
RelWithDebInfo MinSizeRel.")
ENDIF()

set(project_root_dir "../")
include_directories("${project_root_dir}/src")
include_directories("${project_root_dir}/bin")
include_directories("${project_root_dir}/libs")
include_directories("${project_root_dir}/libs/char_ops")
include_directories("${project_root_dir}/libs/num_ops")

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
   "${project_root_dir}/cmake/scttr_cfg.h.in"
   "${project_root_dir}/src/scttr_cfg.h"
   )

message(STATUS "Adding compiler flags")

# add a target to generate API documentation with Doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/../doc/doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/doxyfile.in @ONLY)
  message(STATUS "Creating documentation inside the build directory with Doxygen.")
  add_custom_target(doc
    ALL
    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/doxyfile.in
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
endif(DOXYGEN_FOUND)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  message(STATUS "Building debug-type makefiles")
  set(CMAKE_C_FLAGS "-g -o0 -Wextra -pedantic -D_GNU_SOURCE -std=c89")
else(CMAKE_BUILD_TYPE STREQUAL "Debug")
  message(STATUS "Building release-type makefiles")
  set(CMAKE_C_FLAGS "-Wall -o3 -D_GNU_SOURCE -std=c89")
endif(CMAKE_BUILD_TYPE STREQUAL "Debug")

add_executable(scttr "${project_root_dir}/src/main.c")

add_library(scttr_src
  "${project_root_dir}/src/transitions.c"
  "${project_root_dir}/src/scttr_io.c"
  "${project_root_dir}/src/calc_spec.c"
  "${project_root_dir}/src/spectrum.c"
  )

add_library(scttr_libs
  "${project_root_dir}/libs/dyn_array.c"
  "${project_root_dir}/libs/char_ops/get_nums.c"
  "${project_root_dir}/libs/char_ops/std_char_ops.c"
  "${project_root_dir}/libs/num_ops/iquicks.c"
  "${project_root_dir}/libs/num_ops/std_num_ops.c"
  )

target_link_libraries(scttr scttr_src scttr_libs -lm)
