# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

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
CMAKE_SOURCE_DIR = /wrk3/simploce/pt-cgmd/simulation

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /wrk3/simploce/pt-cgmd/simulation/cmake-build-debug

# Include any dependencies generated for this target.
include tests/CMakeFiles/simulation-test.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include tests/CMakeFiles/simulation-test.dir/compiler_depend.make

# Include the progress variables for this target.
include tests/CMakeFiles/simulation-test.dir/progress.make

# Include the compile flags for this target's objects.
include tests/CMakeFiles/simulation-test.dir/flags.make

tests/CMakeFiles/simulation-test.dir/simulation-test.cpp.o: tests/CMakeFiles/simulation-test.dir/flags.make
tests/CMakeFiles/simulation-test.dir/simulation-test.cpp.o: ../tests/simulation-test.cpp
tests/CMakeFiles/simulation-test.dir/simulation-test.cpp.o: tests/CMakeFiles/simulation-test.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/wrk3/simploce/pt-cgmd/simulation/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object tests/CMakeFiles/simulation-test.dir/simulation-test.cpp.o"
	cd /wrk3/simploce/pt-cgmd/simulation/cmake-build-debug/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT tests/CMakeFiles/simulation-test.dir/simulation-test.cpp.o -MF CMakeFiles/simulation-test.dir/simulation-test.cpp.o.d -o CMakeFiles/simulation-test.dir/simulation-test.cpp.o -c /wrk3/simploce/pt-cgmd/simulation/tests/simulation-test.cpp

tests/CMakeFiles/simulation-test.dir/simulation-test.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/simulation-test.dir/simulation-test.cpp.i"
	cd /wrk3/simploce/pt-cgmd/simulation/cmake-build-debug/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /wrk3/simploce/pt-cgmd/simulation/tests/simulation-test.cpp > CMakeFiles/simulation-test.dir/simulation-test.cpp.i

tests/CMakeFiles/simulation-test.dir/simulation-test.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/simulation-test.dir/simulation-test.cpp.s"
	cd /wrk3/simploce/pt-cgmd/simulation/cmake-build-debug/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /wrk3/simploce/pt-cgmd/simulation/tests/simulation-test.cpp -o CMakeFiles/simulation-test.dir/simulation-test.cpp.s

# Object files for target simulation-test
simulation__test_OBJECTS = \
"CMakeFiles/simulation-test.dir/simulation-test.cpp.o"

# External object files for target simulation-test
simulation__test_EXTERNAL_OBJECTS =

tests/simulation-test: tests/CMakeFiles/simulation-test.dir/simulation-test.cpp.o
tests/simulation-test: tests/CMakeFiles/simulation-test.dir/build.make
tests/simulation-test: libsimulation.so
tests/simulation-test: tests/CMakeFiles/simulation-test.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/wrk3/simploce/pt-cgmd/simulation/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable simulation-test"
	cd /wrk3/simploce/pt-cgmd/simulation/cmake-build-debug/tests && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/simulation-test.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tests/CMakeFiles/simulation-test.dir/build: tests/simulation-test
.PHONY : tests/CMakeFiles/simulation-test.dir/build

tests/CMakeFiles/simulation-test.dir/clean:
	cd /wrk3/simploce/pt-cgmd/simulation/cmake-build-debug/tests && $(CMAKE_COMMAND) -P CMakeFiles/simulation-test.dir/cmake_clean.cmake
.PHONY : tests/CMakeFiles/simulation-test.dir/clean

tests/CMakeFiles/simulation-test.dir/depend:
	cd /wrk3/simploce/pt-cgmd/simulation/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /wrk3/simploce/pt-cgmd/simulation /wrk3/simploce/pt-cgmd/simulation/tests /wrk3/simploce/pt-cgmd/simulation/cmake-build-debug /wrk3/simploce/pt-cgmd/simulation/cmake-build-debug/tests /wrk3/simploce/pt-cgmd/simulation/cmake-build-debug/tests/CMakeFiles/simulation-test.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tests/CMakeFiles/simulation-test.dir/depend

