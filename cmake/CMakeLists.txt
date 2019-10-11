# Copyright (C) 2015 Erik Källman, Steven Vancoillie (git version extraction)
# This file is part of scttr.

# scttr is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# scttr is distributed in the hope that it will be useful,
# but without any warranty; without even the implied warranty of
# merchantability or fitness for a particular purpose. See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with scttr, found in the "license" subdirectory of the root
# directory of the scttr program. If not, see <http://www.gnu.org/licenses/>.

# If the user specifies -DCMAKE_BUILD_TYPE on the command line, take its
# definition and dump it in the cache along with proper documentation,
# otherwise set CMAKE_BUILD_TYPE to debug prior to calling PROJECT()
#

cmake_minimum_required(VERSION 2.8)
enable_language(C)
project(scttr C)

INCLUDE(CMakeForceCompiler)
CMAKE_FORCE_C_COMPILER(gcc GNU)

set (SCTTR_VERSION_MAJOR 0)
set (SCTTR_VERSION_MINOR 3)
set (TEMPLATE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../src/plot_template")
set (COMPILER_TYPE 1) # A flag to let the program know what compiler was used in the build. Add an index and description to this comment for newly supported compilers. 1 = GCC, 2 = Intel

if ( ${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR} )
  message( "In-source builds are allowed but not recommended. In order to build Scttr, create and/or enter the directory where you want the program to be installed, then run CMake from there. You may need to remove the CMakeCache.txt file if CMake has been ran in the directory previously without being deleted." )
endif()

IF(DEFINED CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose the type of
build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) debug
production.")
ELSE()
  SET(CMAKE_BUILD_TYPE production CACHE STRING "Choose the type of build,
options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) debug production.")
ENDIF()

set (PLOT_TEMPLATE_PTH "${PROJECT_SOURCE_DIR}/../src/plot_template")

message ("Detecting scttr version from git repository:")
set (SCTTR_GIT_REPO "${PROJECT_SOURCE_DIR}/../.git")


find_program (GIT "git")
find_program (DOT dot)

if(EXISTS ${DOT})
  set (DOT_EXISTS "YES")
else ()
  set (DOT_EXISTS "NO")
endif ()

mark_as_advanced(FORCE GIT)

if (EXISTS ${GIT} AND EXISTS ${SCTTR_GIT_REPO})
        # always set the HEAD as a source for configuring
        configure_file ("${SCTTR_GIT_REPO}/HEAD" git-head COPYONLY)
        # if it is a ref to a named branch, include the branch ref too
        file (READ "${SCTTR_GIT_REPO}/HEAD" HEAD_STRING)
        string (STRIP "${HEAD_STRING}" HEAD_STRING)
        string (REPLACE "ref: " "" HEAD_REF "${HEAD_STRING}")
        if (EXISTS "${SCTTR_GIT_REPO}/${HEAD_REF}")
          configure_file ("${SCTTR_GIT_REPO}/${HEAD_REF}" git-head-ref COPYONLY)
        endif ()
        # now, get the actual version from the git describe tool
        execute_process (
                COMMAND ${GIT} "describe" "--always" "--dirty" "--match" "v*"
                WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                OUTPUT_VARIABLE SCTTR_VERSION
                ERROR_VARIABLE GIT_DESCRIBE_ERROR
                RESULT_VARIABLE GIT_DESCRIBE_RC
                OUTPUT_STRIP_TRAILING_WHITESPACE
                )
        if (GIT_DESCRIBE_RC)
                message (FATAL_ERROR "Failed to run git: ${GIT_DESCRIBE_ERROR}")
        endif ()
else ()
        set (SCTTR_VERSION "unknown")
endif ()

message ("-- SCTTR_VERSION: ${SCTTR_VERSION}")
configure_file (${PROJECT_SOURCE_DIR}/scttr_cfg.h.in ${CMAKE_CURRENT_BINARY_DIR}/include/scttr_cfg.h)
configure_file (${PROJECT_SOURCE_DIR}/mainpage.h.in ${CMAKE_CURRENT_BINARY_DIR}/doc/mainpage.h)
configure_file(${PROJECT_SOURCE_DIR}/doxyfile.h.in ${CMAKE_CURRENT_BINARY_DIR}/doc/doxyfile.in)


set(project_root_dir "../")
include_directories("${project_root_dir}/src")
include_directories("${CMAKE_CURRENT_BINARY_DIR}/include")
include_directories("${project_root_dir}/bin")
include_directories("${project_root_dir}/libs")
include_directories("${project_root_dir}/libs/char_ops")
include_directories("${project_root_dir}/libs/num_ops")
include_directories("${project_root_dir}/libs/libcpuid")

# add a target to generate API documentation with Doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
  message(STATUS "Creating documentation inside the build directory with Doxygen.")
  add_custom_target(doc
    ALL
    ${DOXYGEN_EXECUTABLE}
    ${CMAKE_CURRENT_BINARY_DIR}/doc/doxyfile.in
    ${CMAKE_CURRENT_BINARY_DIR}/doc/mainpage.h
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc/
    )

endif(DOXYGEN_FOUND)

# check for openmp support
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

if(CMAKE_BUILD_TYPE STREQUAL "debug")
  message(STATUS "Building debug-type makefiles")
  set(CMAKE_C_FLAGS "-pg -fopenmp -o0 -Wall -D_GNU_SOURCE -std=gnu89")

else(CMAKE_BUILD_TYPE STREQUAL "production")
  message(STATUS "Building production-type makefiles")
  set(CMAKE_C_FLAGS "-Wall -fopenmp -o3 -D_GNU_SOURCE -std=gnu89")
endif()

add_executable(scttr "${project_root_dir}/src/main.c")

add_library(scttr_src
  "${project_root_dir}/src/transitions.c"
  "${project_root_dir}/src/scttr_io.c"
  "${project_root_dir}/src/calc_spec.c"
  "${project_root_dir}/src/spectrum.c"
  )

add_library(scttr_std_libs
  "${project_root_dir}/libs/dyn_array.c"
  "${project_root_dir}/libs/timing.c"
  "${project_root_dir}/libs/char_ops/get_nums.c"
  "${project_root_dir}/libs/char_ops/std_char_ops.c"
  "${project_root_dir}/libs/num_ops/iquicks.c"
  "${project_root_dir}/libs/num_ops/std_num_ops.c"
  )

add_library(lib_cache_opt
  "${project_root_dir}/libs/cache_opt.c"
  )

add_library(lib_cpu_opt
  "${project_root_dir}/libs/cpu_opt.c"
  "${project_root_dir}/libs/libcpuid/libcpuid.c"
  "${project_root_dir}/libs/libcpuid/libcpuid_util.c"
  "${project_root_dir}/libs/libcpuid/asm-bits.c"
  "${project_root_dir}/libs/libcpuid/recog_intel.c"
  "${project_root_dir}/libs/libcpuid/recog_amd.c"
  )

target_link_libraries(scttr scttr_src scttr_std_libs lib_cache_opt
  lib_cpu_opt -lm)
