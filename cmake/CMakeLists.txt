# If the user specifies -DCMAKE_BUILD_TYPE on the command line, take its
# definition and dump it in the cache along with proper documentation,
# otherwise set CMAKE_BUILD_TYPE to Debug prior to calling PROJECT()
#
IF(DEFINED CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose the type of
build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug
Release RelWithDebInfo MinSizeRel.")
ELSE()
   SET(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build,
options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release
RelWithDebInfo MinSizeRel.")
ENDIF()

cmake_minimum_required(VERSION 2.8)
project(sctr C)
enable_language(C)
set (SCTR_VERSION_MAJOR 1)
set (SCTR_VERSION_MINOR 5)

set(project_root_dir "../")
include_directories("${project_root_dir}/src")
include_directories("${project_root_dir}/bin")
include_directories("${project_root_dir}/include")
include_directories("${project_root_dir}/include/char_ops")
include_directories("${project_root_dir}/include/num_ops")
include_directories("${project_root_dir}/src/io")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${project_root_dir}/bin")

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
   "${project_root_dir}/cmake/sctr_cfg.h.in"
   "${project_root_dir}/src/sctr_cfg.h"
   )

message(STATUS "GCC detected - Adding compiler flags")

#   set(CMAKE_C_FLAGS "-pthread -Wextra -pedantic -std=c89")

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  message("Building debug-type binaries")
  set(CMAKE_C_FLAGS "-g -o0 -Wextra -pedantic -D_GNU_SOURCE -std=c89")
else(CMAKE_BUILD_TYPE STREQUAL "Debug")
  message("Building release-type binaries")
  set(CMAKE_C_FLAGS "-Wall -D_GNU_SOURCE -std=c89")
endif(CMAKE_BUILD_TYPE STREQUAL "Debug")
# endif(CMAKE_COMPILER_IS_GNUC)

add_executable(sctr "${project_root_dir}/src/main.c")

add_library(sctr_src
  "${project_root_dir}/src/parse_input.c"
  "${project_root_dir}/src/transitions.c"
  "${project_root_dir}/src/spectrum_info.c"
  "${project_root_dir}/src/calc_spec.c"
  )

add_library(sctr_include
  "${project_root_dir}/include/dyn_array.c"

  "${project_root_dir}/include/char_ops/get_numsl.c"
  "${project_root_dir}/include/char_ops/std_char_ops.c"

  "${project_root_dir}/include/num_ops/quicksort.c"
  "${project_root_dir}/include/num_ops/std_num_ops.c"
  )

target_link_libraries(sctr sctr_src sctr_include -lm)
