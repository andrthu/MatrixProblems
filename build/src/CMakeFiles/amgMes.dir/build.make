# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.17

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Disable VCS-based implicit rules.
% : %,v


# Disable VCS-based implicit rules.
% : RCS/%


# Disable VCS-based implicit rules.
% : RCS/%,v


# Disable VCS-based implicit rules.
% : SCCS/s.%


# Disable VCS-based implicit rules.
% : s.%


.SUFFIXES: .hpux_make_needs_suffix_list


# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /cm/shared/ex3-modules/0.5.0/pkgs/cmake-3.17.2/bin/cmake

# The command to remove a file.
RM = /cm/shared/ex3-modules/0.5.0/pkgs/cmake-3.17.2/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/erikhs/MatrixProblems

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/erikhs/MatrixProblems/build

# Include any dependencies generated for this target.
include src/CMakeFiles/amgMes.dir/depend.make

# Include the progress variables for this target.
include src/CMakeFiles/amgMes.dir/progress.make

# Include the compile flags for this target's objects.
include src/CMakeFiles/amgMes.dir/flags.make

src/CMakeFiles/amgMes.dir/amg_meassure.cpp.o: src/CMakeFiles/amgMes.dir/flags.make
src/CMakeFiles/amgMes.dir/amg_meassure.cpp.o: ../src/amg_meassure.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/erikhs/MatrixProblems/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/CMakeFiles/amgMes.dir/amg_meassure.cpp.o"
	cd /home/erikhs/MatrixProblems/build/src && /cm/shared/ex3-modules/0.5.0/pkgs/gcc-10.1.0/bin/g++-10.1  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/amgMes.dir/amg_meassure.cpp.o -c /home/erikhs/MatrixProblems/src/amg_meassure.cpp

src/CMakeFiles/amgMes.dir/amg_meassure.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/amgMes.dir/amg_meassure.cpp.i"
	cd /home/erikhs/MatrixProblems/build/src && /cm/shared/ex3-modules/0.5.0/pkgs/gcc-10.1.0/bin/g++-10.1 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/erikhs/MatrixProblems/src/amg_meassure.cpp > CMakeFiles/amgMes.dir/amg_meassure.cpp.i

src/CMakeFiles/amgMes.dir/amg_meassure.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/amgMes.dir/amg_meassure.cpp.s"
	cd /home/erikhs/MatrixProblems/build/src && /cm/shared/ex3-modules/0.5.0/pkgs/gcc-10.1.0/bin/g++-10.1 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/erikhs/MatrixProblems/src/amg_meassure.cpp -o CMakeFiles/amgMes.dir/amg_meassure.cpp.s

# Object files for target amgMes
amgMes_OBJECTS = \
"CMakeFiles/amgMes.dir/amg_meassure.cpp.o"

# External object files for target amgMes
amgMes_EXTERNAL_OBJECTS =

src/amgMes: src/CMakeFiles/amgMes.dir/amg_meassure.cpp.o
src/amgMes: src/CMakeFiles/amgMes.dir/build.make
src/amgMes: /home/andreast/modules/dune/0.5.0/rpp-mod-2.7.1/lib/libdunegrid.a
src/amgMes: /home/andreast/modules/dune/0.5.0/rpp-mod-2.7.1/lib/libdunegeometry.a
src/amgMes: /home/andreast/modules/dune/0.5.0/rpp-mod-2.7.1/lib/libdunecommon.a
src/amgMes: /cm/shared/ex3-modules/0.5.0/pkgs/suitesparse-5.7.2/lib/libldl.so
src/amgMes: /cm/shared/ex3-modules/0.5.0/pkgs/suitesparse-5.7.2/lib/libspqr.so
src/amgMes: /cm/shared/ex3-modules/0.5.0/pkgs/suitesparse-5.7.2/lib/libumfpack.so
src/amgMes: /cm/shared/ex3-modules/0.5.0/pkgs/suitesparse-5.7.2/lib/libcholmod.so
src/amgMes: /cm/shared/ex3-modules/0.5.0/pkgs/suitesparse-5.7.2/lib/libamd.so
src/amgMes: /cm/shared/ex3-modules/0.5.0/pkgs/suitesparse-5.7.2/lib/libcamd.so
src/amgMes: /cm/shared/ex3-modules/0.5.0/pkgs/suitesparse-5.7.2/lib/libcolamd.so
src/amgMes: /cm/shared/ex3-modules/0.5.0/pkgs/suitesparse-5.7.2/lib/libccolamd.so
src/amgMes: /cm/shared/ex3-modules/0.5.0/pkgs/suitesparse-5.7.2/lib/libsuitesparseconfig.so
src/amgMes: /cm/shared/ex3-modules/0.5.0/pkgs/openblas-0.3.12/lib/libopenblas.so
src/amgMes: /cm/shared/ex3-modules/0.5.0/pkgs/openmpi-4.0.5/lib/libmpi.so
src/amgMes: /home/andreast/modules/flow/0.5.0/rpp-dune-amgcpr-branch-2/lib/libopmsimulators.a
src/amgMes: /cm/shared/ex3-modules/0.5.0/pkgs/boost-1.73.0/lib/libboost_date_time.so
src/amgMes: /cm/shared/ex3-modules/0.5.0/pkgs/suitesparse-5.7.2/lib/libldl.so
src/amgMes: /cm/shared/ex3-modules/0.5.0/pkgs/suitesparse-5.7.2/lib/libspqr.so
src/amgMes: /cm/shared/ex3-modules/0.5.0/pkgs/suitesparse-5.7.2/lib/libcholmod.so
src/amgMes: /cm/shared/ex3-modules/0.5.0/pkgs/suitesparse-5.7.2/lib/libcamd.so
src/amgMes: /cm/shared/ex3-modules/0.5.0/pkgs/suitesparse-5.7.2/lib/libcolamd.so
src/amgMes: /cm/shared/ex3-modules/0.5.0/pkgs/suitesparse-5.7.2/lib/libccolamd.so
src/amgMes: /home/andreast/fork_opm/opm-grid/build-rpp-dune-amgcpr/lib/libopmgrid.a
src/amgMes: /home/andreast/modules/flow/0.5.0/rpp-dune-amgcpr-branch-2/lib/libopmgrid.a
src/amgMes: /home/andreast/modules/dune/0.5.0/rpp-mod-2.7.1/lib/libdunegeometry.a
src/amgMes: /cm/shared/ex3-modules/0.5.0/pkgs/suitesparse-5.7.2/lib/libumfpack.so
src/amgMes: /cm/shared/ex3-modules/0.5.0/pkgs/suitesparse-5.7.2/lib/libamd.so
src/amgMes: /cm/shared/ex3-modules/0.5.0/pkgs/openblas-0.3.12/lib/libopenblas.so
src/amgMes: /usr/lib/x86_64-linux-gnu/libm.so
src/amgMes: /cm/shared/ex3-modules/0.5.0/pkgs/suitesparse-5.7.2/lib/libsuitesparseconfig.so
src/amgMes: /home/andreast/modules/parmetis/0.5.0/4.0.3/lib/libparmetis.a
src/amgMes: /cm/shared/ex3-modules/0.5.0/pkgs/metis-5.1.0/lib/libmetis.so
src/amgMes: /cm/shared/ex3-modules/0.5.0/pkgs/openmpi-4.0.5/lib/libmpi_cxx.so
src/amgMes: /cm/shared/ex3-modules/0.5.0/pkgs/openmpi-4.0.5/lib/libmpi.so
src/amgMes: /home/andreast/modules/dune/0.5.0/rpp-mod-2.7.1/lib/libdunecommon.a
src/amgMes: /home/andreast/modules/dune/0.5.0/rpp-mod-2.7.1/lib/libdunegrid.a
src/amgMes: /home/andreast/modules/flow/0.5.0/rpp-dune-amgcpr-branch-2/lib/libopmcommon.a
src/amgMes: /cm/shared/ex3-modules/0.5.0/pkgs/boost-1.73.0/lib/libboost_system.so
src/amgMes: /cm/shared/ex3-modules/0.5.0/pkgs/boost-1.73.0/lib/libboost_unit_test_framework.so
src/amgMes: /home/andreast/modules/zoltan/0.5.0/12.16/lib/libzoltan.a
src/amgMes: /home/andreast/modules/flow/0.5.0/rpp-dune-amgcpr-branch-2/lib/libopmgrid.a
src/amgMes: /home/andreast/modules/dune/0.5.0/rpp-mod-2.7.1/lib/libdunegeometry.a
src/amgMes: /cm/shared/ex3-modules/0.5.0/pkgs/suitesparse-5.7.2/lib/libumfpack.so
src/amgMes: /cm/shared/ex3-modules/0.5.0/pkgs/suitesparse-5.7.2/lib/libamd.so
src/amgMes: /cm/shared/ex3-modules/0.5.0/pkgs/openblas-0.3.12/lib/libopenblas.so
src/amgMes: /usr/lib/x86_64-linux-gnu/libm.so
src/amgMes: /cm/shared/ex3-modules/0.5.0/pkgs/suitesparse-5.7.2/lib/libsuitesparseconfig.so
src/amgMes: /home/andreast/modules/parmetis/0.5.0/4.0.3/lib/libparmetis.a
src/amgMes: /cm/shared/ex3-modules/0.5.0/pkgs/metis-5.1.0/lib/libmetis.so
src/amgMes: /cm/shared/ex3-modules/0.5.0/pkgs/openmpi-4.0.5/lib/libmpi_cxx.so
src/amgMes: /cm/shared/ex3-modules/0.5.0/pkgs/openmpi-4.0.5/lib/libmpi.so
src/amgMes: /home/andreast/modules/dune/0.5.0/rpp-mod-2.7.1/lib/libdunecommon.a
src/amgMes: /home/andreast/modules/dune/0.5.0/rpp-mod-2.7.1/lib/libdunegrid.a
src/amgMes: /home/andreast/modules/flow/0.5.0/rpp-dune-amgcpr-branch-2/lib/libopmcommon.a
src/amgMes: /cm/shared/ex3-modules/0.5.0/pkgs/boost-1.73.0/lib/libboost_system.so
src/amgMes: /cm/shared/ex3-modules/0.5.0/pkgs/boost-1.73.0/lib/libboost_unit_test_framework.so
src/amgMes: /home/andreast/modules/zoltan/0.5.0/12.16/lib/libzoltan.a
src/amgMes: /cm/shared/ex3-modules/0.5.0/pkgs/openblas-0.3.12/lib/libopenblas.so
src/amgMes: /cm/shared/ex3-modules/0.5.0/pkgs/openmpi-4.0.5/lib/libmpi_cxx.so
src/amgMes: /cm/shared/ex3-modules/0.5.0/pkgs/openmpi-4.0.5/lib/libmpi.so
src/amgMes: /home/andreast/modules/flow/0.5.0/rpp-dune-amgcpr-branch-2/lib/libopmcommon.a
src/amgMes: /cm/shared/ex3-modules/0.5.0/pkgs/boost-1.73.0/lib/libboost_system.so
src/amgMes: /cm/shared/ex3-modules/0.5.0/pkgs/boost-1.73.0/lib/libboost_unit_test_framework.so
src/amgMes: /home/andreast/modules/dune/0.5.0/rpp-mod-2.7.1/lib/libdunecommon.a
src/amgMes: /home/andreast/modules/flow/0.5.0/rpp-dune-amgcpr-branch-2/lib/libopmcommon.a
src/amgMes: /cm/shared/ex3-modules/0.5.0/pkgs/boost-1.73.0/lib/libboost_system.so
src/amgMes: /cm/shared/ex3-modules/0.5.0/pkgs/boost-1.73.0/lib/libboost_unit_test_framework.so
src/amgMes: /cm/shared/ex3-modules/0.5.0/pkgs/boost-1.73.0/lib/libboost_program_options.so
src/amgMes: /home/andreast/modules/dune/0.5.0/rpp-mod-2.7.1/lib/libdunegeometry.a
src/amgMes: /home/andreast/modules/dune/0.5.0/rpp-mod-2.7.1/lib/libdunecommon.a
src/amgMes: /cm/shared/ex3-modules/0.5.0/pkgs/suitesparse-5.7.2/lib/libldl.so
src/amgMes: /cm/shared/ex3-modules/0.5.0/pkgs/suitesparse-5.7.2/lib/libspqr.so
src/amgMes: /cm/shared/ex3-modules/0.5.0/pkgs/suitesparse-5.7.2/lib/libcholmod.so
src/amgMes: /cm/shared/ex3-modules/0.5.0/pkgs/suitesparse-5.7.2/lib/libcamd.so
src/amgMes: /cm/shared/ex3-modules/0.5.0/pkgs/suitesparse-5.7.2/lib/libcolamd.so
src/amgMes: /cm/shared/ex3-modules/0.5.0/pkgs/suitesparse-5.7.2/lib/libccolamd.so
src/amgMes: /home/andreast/modules/flow/0.5.0/rpp-dune-amgcpr-branch-2/lib/libopmsimulators.a
src/amgMes: /cm/shared/ex3-modules/0.5.0/pkgs/boost-1.73.0/lib/libboost_date_time.so
src/amgMes: /home/andreast/fork_opm/opm-grid/build-rpp-dune-amgcpr/lib/libopmgrid.a
src/amgMes: /cm/shared/ex3-modules/0.5.0/pkgs/openblas-0.3.12/lib/libopenblas.so
src/amgMes: /cm/shared/ex3-modules/0.5.0/pkgs/boost-1.73.0/lib/libboost_program_options.so
src/amgMes: /cm/shared/ex3-modules/0.5.0/pkgs/gcc-10.1.0/lib64/libgomp.so
src/amgMes: /lib/x86_64-linux-gnu/libpthread.a
src/amgMes: src/CMakeFiles/amgMes.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/erikhs/MatrixProblems/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable amgMes"
	cd /home/erikhs/MatrixProblems/build/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/amgMes.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/CMakeFiles/amgMes.dir/build: src/amgMes

.PHONY : src/CMakeFiles/amgMes.dir/build

src/CMakeFiles/amgMes.dir/clean:
	cd /home/erikhs/MatrixProblems/build/src && $(CMAKE_COMMAND) -P CMakeFiles/amgMes.dir/cmake_clean.cmake
.PHONY : src/CMakeFiles/amgMes.dir/clean

src/CMakeFiles/amgMes.dir/depend:
	cd /home/erikhs/MatrixProblems/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/erikhs/MatrixProblems /home/erikhs/MatrixProblems/src /home/erikhs/MatrixProblems/build /home/erikhs/MatrixProblems/build/src /home/erikhs/MatrixProblems/build/src/CMakeFiles/amgMes.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/CMakeFiles/amgMes.dir/depend

