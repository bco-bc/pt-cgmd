# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.19

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

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /home/andre/clion-2021.1.3/bin/cmake/linux/bin/cmake

# The command to remove a file.
RM = /home/andre/clion-2021.1.3/bin/cmake/linux/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/andre/wrk3/simploce/pt-cgmd/cpputil

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/andre/wrk3/simploce/pt-cgmd/cpputil/cmake-build-release

# Include any dependencies generated for this target.
include tests/CMakeFiles/matrix_map_test.dir/depend.make

# Include the progress variables for this target.
include tests/CMakeFiles/matrix_map_test.dir/progress.make

# Include the compile flags for this target's objects.
include tests/CMakeFiles/matrix_map_test.dir/flags.make

tests/CMakeFiles/matrix_map_test.dir/matrix-map-test.cpp.o: tests/CMakeFiles/matrix_map_test.dir/flags.make
tests/CMakeFiles/matrix_map_test.dir/matrix-map-test.cpp.o: ../tests/matrix-map-test.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/andre/wrk3/simploce/pt-cgmd/cpputil/cmake-build-release/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object tests/CMakeFiles/matrix_map_test.dir/matrix-map-test.cpp.o"
	cd /home/andre/wrk3/simploce/pt-cgmd/cpputil/cmake-build-release/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/matrix_map_test.dir/matrix-map-test.cpp.o -c /home/andre/wrk3/simploce/pt-cgmd/cpputil/tests/matrix-map-test.cpp

tests/CMakeFiles/matrix_map_test.dir/matrix-map-test.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/matrix_map_test.dir/matrix-map-test.cpp.i"
	cd /home/andre/wrk3/simploce/pt-cgmd/cpputil/cmake-build-release/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/andre/wrk3/simploce/pt-cgmd/cpputil/tests/matrix-map-test.cpp > CMakeFiles/matrix_map_test.dir/matrix-map-test.cpp.i

tests/CMakeFiles/matrix_map_test.dir/matrix-map-test.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/matrix_map_test.dir/matrix-map-test.cpp.s"
	cd /home/andre/wrk3/simploce/pt-cgmd/cpputil/cmake-build-release/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/andre/wrk3/simploce/pt-cgmd/cpputil/tests/matrix-map-test.cpp -o CMakeFiles/matrix_map_test.dir/matrix-map-test.cpp.s

# Object files for target matrix_map_test
matrix_map_test_OBJECTS = \
"CMakeFiles/matrix_map_test.dir/matrix-map-test.cpp.o"

# External object files for target matrix_map_test
matrix_map_test_EXTERNAL_OBJECTS =

tests/matrix_map_test: tests/CMakeFiles/matrix_map_test.dir/matrix-map-test.cpp.o
tests/matrix_map_test: tests/CMakeFiles/matrix_map_test.dir/build.make
tests/matrix_map_test: libcpputil.so
tests/matrix_map_test: tests/CMakeFiles/matrix_map_test.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/andre/wrk3/simploce/pt-cgmd/cpputil/cmake-build-release/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable matrix_map_test"
	cd /home/andre/wrk3/simploce/pt-cgmd/cpputil/cmake-build-release/tests && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/matrix_map_test.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tests/CMakeFiles/matrix_map_test.dir/build: tests/matrix_map_test

.PHONY : tests/CMakeFiles/matrix_map_test.dir/build

tests/CMakeFiles/matrix_map_test.dir/clean:
	cd /home/andre/wrk3/simploce/pt-cgmd/cpputil/cmake-build-release/tests && $(CMAKE_COMMAND) -P CMakeFiles/matrix_map_test.dir/cmake_clean.cmake
.PHONY : tests/CMakeFiles/matrix_map_test.dir/clean

tests/CMakeFiles/matrix_map_test.dir/depend:
	cd /home/andre/wrk3/simploce/pt-cgmd/cpputil/cmake-build-release && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/andre/wrk3/simploce/pt-cgmd/cpputil /home/andre/wrk3/simploce/pt-cgmd/cpputil/tests /home/andre/wrk3/simploce/pt-cgmd/cpputil/cmake-build-release /home/andre/wrk3/simploce/pt-cgmd/cpputil/cmake-build-release/tests /home/andre/wrk3/simploce/pt-cgmd/cpputil/cmake-build-release/tests/CMakeFiles/matrix_map_test.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tests/CMakeFiles/matrix_map_test.dir/depend

