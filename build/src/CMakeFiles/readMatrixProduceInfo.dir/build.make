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
include src/CMakeFiles/readMatrixProduceInfo.dir/depend.make

# Include the progress variables for this target.
include src/CMakeFiles/readMatrixProduceInfo.dir/progress.make

# Include the compile flags for this target's objects.
include src/CMakeFiles/readMatrixProduceInfo.dir/flags.make

src/CMakeFiles/readMatrixProduceInfo.dir/readMatrixProduceInfo.cpp.o: src/CMakeFiles/readMatrixProduceInfo.dir/flags.make
src/CMakeFiles/readMatrixProduceInfo.dir/readMatrixProduceInfo.cpp.o: ../src/readMatrixProduceInfo.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/erikhs/MatrixProblems/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/CMakeFiles/readMatrixProduceInfo.dir/readMatrixProduceInfo.cpp.o"
	cd /home/erikhs/MatrixProblems/build/src && /cm/shared/ex3-modules/0.5.0/pkgs/gcc-10.1.0/bin/g++-10.1  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/readMatrixProduceInfo.dir/readMatrixProduceInfo.cpp.o -c /home/erikhs/MatrixProblems/src/readMatrixProduceInfo.cpp

src/CMakeFiles/readMatrixProduceInfo.dir/readMatrixProduceInfo.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/readMatrixProduceInfo.dir/readMatrixProduceInfo.cpp.i"
	cd /home/erikhs/MatrixProblems/build/src && /cm/shared/ex3-modules/0.5.0/pkgs/gcc-10.1.0/bin/g++-10.1 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/erikhs/MatrixProblems/src/readMatrixProduceInfo.cpp > CMakeFiles/readMatrixProduceInfo.dir/readMatrixProduceInfo.cpp.i

src/CMakeFiles/readMatrixProduceInfo.dir/readMatrixProduceInfo.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/readMatrixProduceInfo.dir/readMatrixProduceInfo.cpp.s"
	cd /home/erikhs/MatrixProblems/build/src && /cm/shared/ex3-modules/0.5.0/pkgs/gcc-10.1.0/bin/g++-10.1 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/erikhs/MatrixProblems/src/readMatrixProduceInfo.cpp -o CMakeFiles/readMatrixProduceInfo.dir/readMatrixProduceInfo.cpp.s

# Object files for target readMatrixProduceInfo
readMatrixProduceInfo_OBJECTS = \
"CMakeFiles/readMatrixProduceInfo.dir/readMatrixProduceInfo.cpp.o"

# External object files for target readMatrixProduceInfo
readMatrixProduceInfo_EXTERNAL_OBJECTS =

src/readMatrixProduceInfo: src/CMakeFiles/readMatrixProduceInfo.dir/readMatrixProduceInfo.cpp.o
src/readMatrixProduceInfo: src/CMakeFiles/readMatrixProduceInfo.dir/build.make
src/readMatrixProduceInfo: /home/andreast/modules/dune/0.5.0/rpp-mod-2.7.1/lib/libdunegrid.a
src/readMatrixProduceInfo: /home/andreast/modules/dune/0.5.0/rpp-mod-2.7.1/lib/libdunegeometry.a
src/readMatrixProduceInfo: /home/andreast/modules/dune/0.5.0/rpp-mod-2.7.1/lib/libdunecommon.a
src/readMatrixProduceInfo: /cm/shared/ex3-modules/0.5.0/pkgs/suitesparse-5.7.2/lib/libldl.so
src/readMatrixProduceInfo: /cm/shared/ex3-modules/0.5.0/pkgs/suitesparse-5.7.2/lib/libspqr.so
src/readMatrixProduceInfo: /cm/shared/ex3-modules/0.5.0/pkgs/suitesparse-5.7.2/lib/libumfpack.so
src/readMatrixProduceInfo: /cm/shared/ex3-modules/0.5.0/pkgs/suitesparse-5.7.2/lib/libcholmod.so
src/readMatrixProduceInfo: /cm/shared/ex3-modules/0.5.0/pkgs/suitesparse-5.7.2/lib/libamd.so
src/readMatrixProduceInfo: /cm/shared/ex3-modules/0.5.0/pkgs/suitesparse-5.7.2/lib/libcamd.so
src/readMatrixProduceInfo: /cm/shared/ex3-modules/0.5.0/pkgs/suitesparse-5.7.2/lib/libcolamd.so
src/readMatrixProduceInfo: /cm/shared/ex3-modules/0.5.0/pkgs/suitesparse-5.7.2/lib/libccolamd.so
src/readMatrixProduceInfo: /cm/shared/ex3-modules/0.5.0/pkgs/suitesparse-5.7.2/lib/libsuitesparseconfig.so
src/readMatrixProduceInfo: /cm/shared/ex3-modules/0.5.0/pkgs/openblas-0.3.12/lib/libopenblas.so
src/readMatrixProduceInfo: /cm/shared/ex3-modules/0.5.0/pkgs/openmpi-4.0.5/lib/libmpi.so
src/readMatrixProduceInfo: /home/andreast/modules/dune/0.5.0/rpp-mod-2.7.1/lib/libdunegrid.a
src/readMatrixProduceInfo: /home/andreast/modules/dune/0.5.0/rpp-mod-2.7.1/lib/libdunegeometry.a
src/readMatrixProduceInfo: /home/andreast/modules/dune/0.5.0/rpp-mod-2.7.1/lib/libdunecommon.a
src/readMatrixProduceInfo: /home/andreast/modules/zoltan/0.5.0/12.16/lib/libzoltan.a
src/readMatrixProduceInfo: /cm/shared/ex3-modules/0.5.0/pkgs/openmpi-4.0.5/lib/libmpi_cxx.so
src/readMatrixProduceInfo: /cm/shared/ex3-modules/0.5.0/pkgs/openmpi-4.0.5/lib/libmpi.so
src/readMatrixProduceInfo: /cm/shared/ex3-modules/0.5.0/pkgs/boost-1.73.0/lib/libboost_program_options.so
src/readMatrixProduceInfo: /cm/shared/ex3-modules/0.5.0/pkgs/openblas-0.3.12/lib/libopenblas.so
src/readMatrixProduceInfo: /cm/shared/ex3-modules/0.5.0/pkgs/openmpi-4.0.5/lib/libmpi_cxx.so
src/readMatrixProduceInfo: /cm/shared/ex3-modules/0.5.0/pkgs/boost-1.73.0/lib/libboost_program_options.so
src/readMatrixProduceInfo: src/CMakeFiles/readMatrixProduceInfo.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/erikhs/MatrixProblems/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable readMatrixProduceInfo"
	cd /home/erikhs/MatrixProblems/build/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/readMatrixProduceInfo.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/CMakeFiles/readMatrixProduceInfo.dir/build: src/readMatrixProduceInfo

.PHONY : src/CMakeFiles/readMatrixProduceInfo.dir/build

src/CMakeFiles/readMatrixProduceInfo.dir/clean:
	cd /home/erikhs/MatrixProblems/build/src && $(CMAKE_COMMAND) -P CMakeFiles/readMatrixProduceInfo.dir/cmake_clean.cmake
.PHONY : src/CMakeFiles/readMatrixProduceInfo.dir/clean

src/CMakeFiles/readMatrixProduceInfo.dir/depend:
	cd /home/erikhs/MatrixProblems/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/erikhs/MatrixProblems /home/erikhs/MatrixProblems/src /home/erikhs/MatrixProblems/build /home/erikhs/MatrixProblems/build/src /home/erikhs/MatrixProblems/build/src/CMakeFiles/readMatrixProduceInfo.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/CMakeFiles/readMatrixProduceInfo.dir/depend

