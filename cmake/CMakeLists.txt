cmake_minimum_required(VERSION 2.8)
project(smap C)

set (SMAP_VERSION_MAJOR 0)
set (SMAP_VERSION_MINOR 1)

set(project_root_dir "../")
include_directories("${project_root_dir}/src")
include_directories("${project_root_dir}/bin")
include_directories("${project_root_dir}/include")
include_directories("${project_root_dir}/src/io")
include_directories("${project_root_dir}/src/char_ops")
include_directories("${project_root_dir}/src/num_ops")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${project_root_dir}/bin")

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
   "${project_root_dir}/cmake/smap_cfg.h.in"
   "${project_root_dir}/src/smap_cfg.h"
   )

if(CMAKE_COMPILER_IS_GNUCXX)
    message(STATUS "GCC detected - Adding compiler flags")
#    set(CMAKE_C_FLAGS "-Wextra -pedantic -D_GNU_SOURCE -std=c89")
#   set(CMAKE_C_FLAGS "-pthread -Wextra -pedantic -std=c89")
# set(CMAKE_C_FLAGS "-g -D_GNU_SOURCE -std=c89")
#    set(CMAKE_C_FLAGS "-Wall -pedantic -D_GNU_SOURCE -std=c89")
set(CMAKE_C_FLAGS "-g -o0 -Wextra -pedantic -D_GNU_SOURCE -std=c89")
endif(CMAKE_COMPILER_IS_GNUCXX)

add_executable(smap "${project_root_dir}/src/main.c")

add_library(smap_io
  "${project_root_dir}/src/io/parse_input.c"
)

add_library(smap_std_inc
  "${project_root_dir}/src/num_ops/std_num_ops.c"

  "${project_root_dir}/src/char_ops/get_numsl.c"
  "${project_root_dir}/src/char_ops/std_char_ops.c"

  "${project_root_dir}/src/appc.c"
  "${project_root_dir}/src/info_node.c"
  )

add_library(smap_smap
  "${project_root_dir}/src/smap.c"
  "${project_root_dir}/src/spec.c"
  )

add_library(smap_da
  "${project_root_dir}/src/dynarray.c"
)

add_library(smap_es
  "${project_root_dir}/src/estate.c"
  "${project_root_dir}/src/num_ops/k_meansl.c"
)

target_link_libraries(smap smap_io smap_std_inc smap_smap smap_da smap_es -lm)
