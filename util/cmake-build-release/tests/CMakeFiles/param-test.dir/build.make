# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.20

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
CMAKE_COMMAND = /localdisk/clion-2021.2.3/bin/cmake/linux/bin/cmake

# The command to remove a file.
RM = /localdisk/clion-2021.2.3/bin/cmake/linux/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /wrk3/simploce/pt-cgmd/util

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /wrk3/simploce/pt-cgmd/util/cmake-build-release

# Include any dependencies generated for this target.
include tests/CMakeFiles/param-test.dir/depend.make
# Include the progress variables for this target.
include tests/CMakeFiles/param-test.dir/progress.make

# Include the compile flags for this target's objects.
include tests/CMakeFiles/param-test.dir/flags.make

tests/CMakeFiles/param-test.dir/param-test.cpp.o: tests/CMakeFiles/param-test.dir/flags.make
tests/CMakeFiles/param-test.dir/param-test.cpp.o: ../tests/param-test.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/wrk3/simploce/pt-cgmd/util/cmake-build-release/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object tests/CMakeFiles/param-test.dir/param-test.cpp.o"
	cd /wrk3/simploce/pt-cgmd/util/cmake-build-release/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/param-test.dir/param-test.cpp.o -c /wrk3/simploce/pt-cgmd/util/tests/param-test.cpp

tests/CMakeFiles/param-test.dir/param-test.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/param-test.dir/param-test.cpp.i"
	cd /wrk3/simploce/pt-cgmd/util/cmake-build-release/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /wrk3/simploce/pt-cgmd/util/tests/param-test.cpp > CMakeFiles/param-test.dir/param-test.cpp.i

tests/CMakeFiles/param-test.dir/param-test.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/param-test.dir/param-test.cpp.s"
	cd /wrk3/simploce/pt-cgmd/util/cmake-build-release/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /wrk3/simploce/pt-cgmd/util/tests/param-test.cpp -o CMakeFiles/param-test.dir/param-test.cpp.s

# Object files for target param-test
param__test_OBJECTS = \
"CMakeFiles/param-test.dir/param-test.cpp.o"

# External object files for target param-test
param__test_EXTERNAL_OBJECTS =

tests/param-test: tests/CMakeFiles/param-test.dir/param-test.cpp.o
tests/param-test: tests/CMakeFiles/param-test.dir/build.make
tests/param-test: libutil.so
tests/param-test: tests/CMakeFiles/param-test.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/wrk3/simploce/pt-cgmd/util/cmake-build-release/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable param-test"
	cd /wrk3/simploce/pt-cgmd/util/cmake-build-release/tests && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/param-test.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tests/CMakeFiles/param-test.dir/build: tests/param-test
.PHONY : tests/CMakeFiles/param-test.dir/build

tests/CMakeFiles/param-test.dir/clean:
	cd /wrk3/simploce/pt-cgmd/util/cmake-build-release/tests && $(CMAKE_COMMAND) -P CMakeFiles/param-test.dir/cmake_clean.cmake
.PHONY : tests/CMakeFiles/param-test.dir/clean

tests/CMakeFiles/param-test.dir/depend:
	cd /wrk3/simploce/pt-cgmd/util/cmake-build-release && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /wrk3/simploce/pt-cgmd/util /wrk3/simploce/pt-cgmd/util/tests /wrk3/simploce/pt-cgmd/util/cmake-build-release /wrk3/simploce/pt-cgmd/util/cmake-build-release/tests /wrk3/simploce/pt-cgmd/util/cmake-build-release/tests/CMakeFiles/param-test.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tests/CMakeFiles/param-test.dir/depend

