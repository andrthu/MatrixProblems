if(NOT MatrixProblems_FOUND)
# Whether this module is installed or not
set(MatrixProblems_INSTALLED OFF)

# Settings specific to the module

# Package initialization
# Set prefix to source dir
set(PACKAGE_PREFIX_DIR /home/erikhs/MatrixProblems)
macro(set_and_check _var _file)
  set(${_var} "${_file}")
  if(NOT EXISTS "${_file}")
    message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
  endif()
endmacro()

#report other information
set_and_check(MatrixProblems_PREFIX "${PACKAGE_PREFIX_DIR}")
set_and_check(MatrixProblems_INCLUDE_DIRS "/home/erikhs/MatrixProblems")
set(MatrixProblems_CXX_FLAGS "-std=c++17 -O2 -mavx -funroll-loops")
set(MatrixProblems_CXX_FLAGS_DEBUG "-g")
set(MatrixProblems_CXX_FLAGS_MINSIZEREL "-Os -DNDEBUG")
set(MatrixProblems_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set(MatrixProblems_CXX_FLAGS_RELWITHDEBINFO "-O2 -g -DNDEBUG")
set(MatrixProblems_DEPENDS "dune-common;dune-istl;dune-geometry;dune-grid")
set(MatrixProblems_SUGGESTS "")
set(MatrixProblems_MODULE_PATH "/home/erikhs/MatrixProblems/cmake/modules")
set(MatrixProblems_LIBRARIES "")

# Lines that are set by the CMake build system via the variable DUNE_CUSTOM_PKG_CONFIG_SECTION


#import the target
if(MatrixProblems_LIBRARIES)
  get_filename_component(_dir "${CMAKE_CURRENT_LIST_FILE}" PATH)
  include("${_dir}/MatrixProblems-targets.cmake")
endif()
endif()
