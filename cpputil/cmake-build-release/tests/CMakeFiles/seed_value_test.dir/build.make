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
include tests/CMakeFiles/seed_value_test.dir/depend.make

# Include the progress variables for this target.
include tests/CMakeFiles/seed_value_test.dir/progress.make

# Include the compile flags for this target's objects.
include tests/CMakeFiles/seed_value_test.dir/flags.make

tests/CMakeFiles/seed_value_test.dir/seed-value-test.cpp.o: tests/CMakeFiles/seed_value_test.dir/flags.make
tests/CMakeFiles/seed_value_test.dir/seed-value-test.cpp.o: ../tests/seed-value-test.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/andre/wrk3/simploce/pt-cgmd/cpputil/cmake-build-release/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object tests/CMakeFiles/seed_value_test.dir/seed-value-test.cpp.o"
	cd /home/andre/wrk3/simploce/pt-cgmd/cpputil/cmake-build-release/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/seed_value_test.dir/seed-value-test.cpp.o -c /home/andre/wrk3/simploce/pt-cgmd/cpputil/tests/seed-value-test.cpp

tests/CMakeFiles/seed_value_test.dir/seed-value-test.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/seed_value_test.dir/seed-value-test.cpp.i"
	cd /home/andre/wrk3/simploce/pt-cgmd/cpputil/cmake-build-release/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/andre/wrk3/simploce/pt-cgmd/cpputil/tests/seed-value-test.cpp > CMakeFiles/seed_value_test.dir/seed-value-test.cpp.i

tests/CMakeFiles/seed_value_test.dir/seed-value-test.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/seed_value_test.dir/seed-value-test.cpp.s"
	cd /home/andre/wrk3/simploce/pt-cgmd/cpputil/cmake-build-release/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/andre/wrk3/simploce/pt-cgmd/cpputil/tests/seed-value-test.cpp -o CMakeFiles/seed_value_test.dir/seed-value-test.cpp.s

# Object files for target seed_value_test
seed_value_test_OBJECTS = \
"CMakeFiles/seed_value_test.dir/seed-value-test.cpp.o"

# External object files for target seed_value_test
seed_value_test_EXTERNAL_OBJECTS =

tests/seed_value_test: tests/CMakeFiles/seed_value_test.dir/seed-value-test.cpp.o
tests/seed_value_test: tests/CMakeFiles/seed_value_test.dir/build.make
tests/seed_value_test: libcpputil.so
tests/seed_value_test: tests/CMakeFiles/seed_value_test.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/andre/wrk3/simploce/pt-cgmd/cpputil/cmake-build-release/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable seed_value_test"
	cd /home/andre/wrk3/simploce/pt-cgmd/cpputil/cmake-build-release/tests && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/seed_value_test.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tests/CMakeFiles/seed_value_test.dir/build: tests/seed_value_test

.PHONY : tests/CMakeFiles/seed_value_test.dir/build

tests/CMakeFiles/seed_value_test.dir/clean:
	cd /home/andre/wrk3/simploce/pt-cgmd/cpputil/cmake-build-release/tests && $(CMAKE_COMMAND) -P CMakeFiles/seed_value_test.dir/cmake_clean.cmake
.PHONY : tests/CMakeFiles/seed_value_test.dir/clean

tests/CMakeFiles/seed_value_test.dir/depend:
	cd /home/andre/wrk3/simploce/pt-cgmd/cpputil/cmake-build-release && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/andre/wrk3/simploce/pt-cgmd/cpputil /home/andre/wrk3/simploce/pt-cgmd/cpputil/tests /home/andre/wrk3/simploce/pt-cgmd/cpputil/cmake-build-release /home/andre/wrk3/simploce/pt-cgmd/cpputil/cmake-build-release/tests /home/andre/wrk3/simploce/pt-cgmd/cpputil/cmake-build-release/tests/CMakeFiles/seed_value_test.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tests/CMakeFiles/seed_value_test.dir/depend
