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
include src/CMakeFiles/betterIOTest2.dir/depend.make

# Include the progress variables for this target.
include src/CMakeFiles/betterIOTest2.dir/progress.make

# Include the compile flags for this target's objects.
include src/CMakeFiles/betterIOTest2.dir/flags.make

src/CMakeFiles/betterIOTest2.dir/betterIOTest2.cpp.o: src/CMakeFiles/betterIOTest2.dir/flags.make
src/CMakeFiles/betterIOTest2.dir/betterIOTest2.cpp.o: ../src/betterIOTest2.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/erikhs/MatrixProblems/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/CMakeFiles/betterIOTest2.dir/betterIOTest2.cpp.o"
	cd /home/erikhs/MatrixProblems/build/src && /cm/shared/ex3-modules/0.5.0/pkgs/gcc-10.1.0/bin/g++-10.1  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/betterIOTest2.dir/betterIOTest2.cpp.o -c /home/erikhs/MatrixProblems/src/betterIOTest2.cpp

src/CMakeFiles/betterIOTest2.dir/betterIOTest2.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/betterIOTest2.dir/betterIOTest2.cpp.i"
	cd /home/erikhs/MatrixProblems/build/src && /cm/shared/ex3-modules/0.5.0/pkgs/gcc-10.1.0/bin/g++-10.1 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/erikhs/MatrixProblems/src/betterIOTest2.cpp > CMakeFiles/betterIOTest2.dir/betterIOTest2.cpp.i

src/CMakeFiles/betterIOTest2.dir/betterIOTest2.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/betterIOTest2.dir/betterIOTest2.cpp.s"
	cd /home/erikhs/MatrixProblems/build/src && /cm/shared/ex3-modules/0.5.0/pkgs/gcc-10.1.0/bin/g++-10.1 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/erikhs/MatrixProblems/src/betterIOTest2.cpp -o CMakeFiles/betterIOTest2.dir/betterIOTest2.cpp.s

# Object files for target betterIOTest2
betterIOTest2_OBJECTS = \
"CMakeFiles/betterIOTest2.dir/betterIOTest2.cpp.o"

# External object files for target betterIOTest2
betterIOTest2_EXTERNAL_OBJECTS =

src/betterIOTest2: src/CMakeFiles/betterIOTest2.dir/betterIOTest2.cpp.o
src/betterIOTest2: src/CMakeFiles/betterIOTest2.dir/build.make
src/betterIOTest2: /home/andreast/modules/dune/0.5.0/rpp-mod-2.7.1/lib/libdunegrid.a
src/betterIOTest2: /home/andreast/modules/dune/0.5.0/rpp-mod-2.7.1/lib/libdunegeometry.a
src/betterIOTest2: /home/andreast/modules/dune/0.5.0/rpp-mod-2.7.1/lib/libdunecommon.a
src/betterIOTest2: /cm/shared/ex3-modules/0.5.0/pkgs/suitesparse-5.7.2/lib/libldl.so
src/betterIOTest2: /cm/shared/ex3-modules/0.5.0/pkgs/suitesparse-5.7.2/lib/libspqr.so
src/betterIOTest2: /cm/shared/ex3-modules/0.5.0/pkgs/suitesparse-5.7.2/lib/libumfpack.so
src/betterIOTest2: /cm/shared/ex3-modules/0.5.0/pkgs/suitesparse-5.7.2/lib/libcholmod.so
src/betterIOTest2: /cm/shared/ex3-modules/0.5.0/pkgs/suitesparse-5.7.2/lib/libamd.so
src/betterIOTest2: /cm/shared/ex3-modules/0.5.0/pkgs/suitesparse-5.7.2/lib/libcamd.so
src/betterIOTest2: /cm/shared/ex3-modules/0.5.0/pkgs/suitesparse-5.7.2/lib/libcolamd.so
src/betterIOTest2: /cm/shared/ex3-modules/0.5.0/pkgs/suitesparse-5.7.2/lib/libccolamd.so
src/betterIOTest2: /cm/shared/ex3-modules/0.5.0/pkgs/suitesparse-5.7.2/lib/libsuitesparseconfig.so
src/betterIOTest2: /cm/shared/ex3-modules/0.5.0/pkgs/openblas-0.3.12/lib/libopenblas.so
src/betterIOTest2: /cm/shared/ex3-modules/0.5.0/pkgs/openmpi-4.0.5/lib/libmpi.so
src/betterIOTest2: /home/andreast/modules/flow/0.5.0/rpp-dune-amgcpr-branch-2/lib/libopmsimulators.a
src/betterIOTest2: /cm/shared/ex3-modules/0.5.0/pkgs/boost-1.73.0/lib/libboost_date_time.so
src/betterIOTest2: /cm/shared/ex3-modules/0.5.0/pkgs/suitesparse-5.7.2/lib/libldl.so
src/betterIOTest2: /cm/shared/ex3-modules/0.5.0/pkgs/suitesparse-5.7.2/lib/libspqr.so
src/betterIOTest2: /cm/shared/ex3-modules/0.5.0/pkgs/suitesparse-5.7.2/lib/libcholmod.so
src/betterIOTest2: /cm/shared/ex3-modules/0.5.0/pkgs/suitesparse-5.7.2/lib/libcamd.so
src/betterIOTest2: /cm/shared/ex3-modules/0.5.0/pkgs/suitesparse-5.7.2/lib/libcolamd.so
src/betterIOTest2: /cm/shared/ex3-modules/0.5.0/pkgs/suitesparse-5.7.2/lib/libccolamd.so
src/betterIOTest2: /home/andreast/fork_opm/opm-grid/build-rpp-dune-amgcpr/lib/libopmgrid.a
src/betterIOTest2: /home/andreast/modules/flow/0.5.0/rpp-dune-amgcpr-branch-2/lib/libopmgrid.a
src/betterIOTest2: /home/andreast/modules/dune/0.5.0/rpp-mod-2.7.1/lib/libdunegeometry.a
src/betterIOTest2: /cm/shared/ex3-modules/0.5.0/pkgs/suitesparse-5.7.2/lib/libumfpack.so
src/betterIOTest2: /cm/shared/ex3-modules/0.5.0/pkgs/suitesparse-5.7.2/lib/libamd.so
src/betterIOTest2: /cm/shared/ex3-modules/0.5.0/pkgs/openblas-0.3.12/lib/libopenblas.so
src/betterIOTest2: /usr/lib/x86_64-linux-gnu/libm.so
src/betterIOTest2: /cm/shared/ex3-modules/0.5.0/pkgs/suitesparse-5.7.2/lib/libsuitesparseconfig.so
src/betterIOTest2: /home/andreast/modules/parmetis/0.5.0/4.0.3/lib/libparmetis.a
src/betterIOTest2: /cm/shared/ex3-modules/0.5.0/pkgs/metis-5.1.0/lib/libmetis.so
src/betterIOTest2: /cm/shared/ex3-modules/0.5.0/pkgs/openmpi-4.0.5/lib/libmpi_cxx.so
src/betterIOTest2: /cm/shared/ex3-modules/0.5.0/pkgs/openmpi-4.0.5/lib/libmpi.so
src/betterIOTest2: /home/andreast/modules/dune/0.5.0/rpp-mod-2.7.1/lib/libdunecommon.a
src/betterIOTest2: /home/andreast/modules/dune/0.5.0/rpp-mod-2.7.1/lib/libdunegrid.a
src/betterIOTest2: /home/andreast/modules/flow/0.5.0/rpp-dune-amgcpr-branch-2/lib/libopmcommon.a
src/betterIOTest2: /cm/shared/ex3-modules/0.5.0/pkgs/boost-1.73.0/lib/libboost_system.so
src/betterIOTest2: /cm/shared/ex3-modules/0.5.0/pkgs/boost-1.73.0/lib/libboost_unit_test_framework.so
src/betterIOTest2: /home/andreast/modules/zoltan/0.5.0/12.16/lib/libzoltan.a
src/betterIOTest2: /home/andreast/modules/flow/0.5.0/rpp-dune-amgcpr-branch-2/lib/libopmgrid.a
src/betterIOTest2: /home/andreast/modules/dune/0.5.0/rpp-mod-2.7.1/lib/libdunegeometry.a
src/betterIOTest2: /cm/shared/ex3-modules/0.5.0/pkgs/suitesparse-5.7.2/lib/libumfpack.so
src/betterIOTest2: /cm/shared/ex3-modules/0.5.0/pkgs/suitesparse-5.7.2/lib/libamd.so
src/betterIOTest2: /cm/shared/ex3-modules/0.5.0/pkgs/openblas-0.3.12/lib/libopenblas.so
src/betterIOTest2: /usr/lib/x86_64-linux-gnu/libm.so
src/betterIOTest2: /cm/shared/ex3-modules/0.5.0/pkgs/suitesparse-5.7.2/lib/libsuitesparseconfig.so
src/betterIOTest2: /home/andreast/modules/parmetis/0.5.0/4.0.3/lib/libparmetis.a
src/betterIOTest2: /cm/shared/ex3-modules/0.5.0/pkgs/metis-5.1.0/lib/libmetis.so
src/betterIOTest2: /cm/shared/ex3-modules/0.5.0/pkgs/openmpi-4.0.5/lib/libmpi_cxx.so
src/betterIOTest2: /cm/shared/ex3-modules/0.5.0/pkgs/openmpi-4.0.5/lib/libmpi.so
src/betterIOTest2: /home/andreast/modules/dune/0.5.0/rpp-mod-2.7.1/lib/libdunecommon.a
src/betterIOTest2: /home/andreast/modules/dune/0.5.0/rpp-mod-2.7.1/lib/libdunegrid.a
src/betterIOTest2: /home/andreast/modules/flow/0.5.0/rpp-dune-amgcpr-branch-2/lib/libopmcommon.a
src/betterIOTest2: /cm/shared/ex3-modules/0.5.0/pkgs/boost-1.73.0/lib/libboost_system.so
src/betterIOTest2: /cm/shared/ex3-modules/0.5.0/pkgs/boost-1.73.0/lib/libboost_unit_test_framework.so
src/betterIOTest2: /home/andreast/modules/zoltan/0.5.0/12.16/lib/libzoltan.a
src/betterIOTest2: /cm/shared/ex3-modules/0.5.0/pkgs/openblas-0.3.12/lib/libopenblas.so
src/betterIOTest2: /cm/shared/ex3-modules/0.5.0/pkgs/openmpi-4.0.5/lib/libmpi_cxx.so
src/betterIOTest2: /cm/shared/ex3-modules/0.5.0/pkgs/openmpi-4.0.5/lib/libmpi.so
src/betterIOTest2: /home/andreast/modules/flow/0.5.0/rpp-dune-amgcpr-branch-2/lib/libopmcommon.a
src/betterIOTest2: /cm/shared/ex3-modules/0.5.0/pkgs/boost-1.73.0/lib/libboost_system.so
src/betterIOTest2: /cm/shared/ex3-modules/0.5.0/pkgs/boost-1.73.0/lib/libboost_unit_test_framework.so
src/betterIOTest2: /home/andreast/modules/dune/0.5.0/rpp-mod-2.7.1/lib/libdunecommon.a
src/betterIOTest2: /home/andreast/modules/flow/0.5.0/rpp-dune-amgcpr-branch-2/lib/libopmcommon.a
src/betterIOTest2: /cm/shared/ex3-modules/0.5.0/pkgs/boost-1.73.0/lib/libboost_system.so
src/betterIOTest2: /cm/shared/ex3-modules/0.5.0/pkgs/boost-1.73.0/lib/libboost_unit_test_framework.so
src/betterIOTest2: /cm/shared/ex3-modules/0.5.0/pkgs/boost-1.73.0/lib/libboost_program_options.so
src/betterIOTest2: /home/andreast/modules/dune/0.5.0/rpp-mod-2.7.1/lib/libdunegeometry.a
src/betterIOTest2: /home/andreast/modules/dune/0.5.0/rpp-mod-2.7.1/lib/libdunecommon.a
src/betterIOTest2: /cm/shared/ex3-modules/0.5.0/pkgs/suitesparse-5.7.2/lib/libldl.so
src/betterIOTest2: /cm/shared/ex3-modules/0.5.0/pkgs/suitesparse-5.7.2/lib/libspqr.so
src/betterIOTest2: /cm/shared/ex3-modules/0.5.0/pkgs/suitesparse-5.7.2/lib/libcholmod.so
src/betterIOTest2: /cm/shared/ex3-modules/0.5.0/pkgs/suitesparse-5.7.2/lib/libcamd.so
src/betterIOTest2: /cm/shared/ex3-modules/0.5.0/pkgs/suitesparse-5.7.2/lib/libcolamd.so
src/betterIOTest2: /cm/shared/ex3-modules/0.5.0/pkgs/suitesparse-5.7.2/lib/libccolamd.so
src/betterIOTest2: /home/andreast/modules/flow/0.5.0/rpp-dune-amgcpr-branch-2/lib/libopmsimulators.a
src/betterIOTest2: /cm/shared/ex3-modules/0.5.0/pkgs/boost-1.73.0/lib/libboost_date_time.so
src/betterIOTest2: /home/andreast/fork_opm/opm-grid/build-rpp-dune-amgcpr/lib/libopmgrid.a
src/betterIOTest2: /cm/shared/ex3-modules/0.5.0/pkgs/openblas-0.3.12/lib/libopenblas.so
src/betterIOTest2: /cm/shared/ex3-modules/0.5.0/pkgs/boost-1.73.0/lib/libboost_program_options.so
src/betterIOTest2: /cm/shared/ex3-modules/0.5.0/pkgs/gcc-10.1.0/lib64/libgomp.so
src/betterIOTest2: /lib/x86_64-linux-gnu/libpthread.a
src/betterIOTest2: src/CMakeFiles/betterIOTest2.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/erikhs/MatrixProblems/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable betterIOTest2"
	cd /home/erikhs/MatrixProblems/build/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/betterIOTest2.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/CMakeFiles/betterIOTest2.dir/build: src/betterIOTest2

.PHONY : src/CMakeFiles/betterIOTest2.dir/build

src/CMakeFiles/betterIOTest2.dir/clean:
	cd /home/erikhs/MatrixProblems/build/src && $(CMAKE_COMMAND) -P CMakeFiles/betterIOTest2.dir/cmake_clean.cmake
.PHONY : src/CMakeFiles/betterIOTest2.dir/clean

src/CMakeFiles/betterIOTest2.dir/depend:
	cd /home/erikhs/MatrixProblems/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/erikhs/MatrixProblems /home/erikhs/MatrixProblems/src /home/erikhs/MatrixProblems/build /home/erikhs/MatrixProblems/build/src /home/erikhs/MatrixProblems/build/src/CMakeFiles/betterIOTest2.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/CMakeFiles/betterIOTest2.dir/depend

