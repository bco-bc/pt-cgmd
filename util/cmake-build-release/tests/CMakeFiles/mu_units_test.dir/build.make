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
include tests/CMakeFiles/mu_units_test.dir/depend.make
# Include the progress variables for this target.
include tests/CMakeFiles/mu_units_test.dir/progress.make

# Include the compile flags for this target's objects.
include tests/CMakeFiles/mu_units_test.dir/flags.make

tests/CMakeFiles/mu_units_test.dir/units-mu-test.cpp.o: tests/CMakeFiles/mu_units_test.dir/flags.make
tests/CMakeFiles/mu_units_test.dir/units-mu-test.cpp.o: ../tests/units-mu-test.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/wrk3/simploce/pt-cgmd/util/cmake-build-release/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object tests/CMakeFiles/mu_units_test.dir/units-mu-test.cpp.o"
	cd /wrk3/simploce/pt-cgmd/util/cmake-build-release/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mu_units_test.dir/units-mu-test.cpp.o -c /wrk3/simploce/pt-cgmd/util/tests/units-mu-test.cpp

tests/CMakeFiles/mu_units_test.dir/units-mu-test.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mu_units_test.dir/units-mu-test.cpp.i"
	cd /wrk3/simploce/pt-cgmd/util/cmake-build-release/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /wrk3/simploce/pt-cgmd/util/tests/units-mu-test.cpp > CMakeFiles/mu_units_test.dir/units-mu-test.cpp.i

tests/CMakeFiles/mu_units_test.dir/units-mu-test.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mu_units_test.dir/units-mu-test.cpp.s"
	cd /wrk3/simploce/pt-cgmd/util/cmake-build-release/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /wrk3/simploce/pt-cgmd/util/tests/units-mu-test.cpp -o CMakeFiles/mu_units_test.dir/units-mu-test.cpp.s

# Object files for target mu_units_test
mu_units_test_OBJECTS = \
"CMakeFiles/mu_units_test.dir/units-mu-test.cpp.o"

# External object files for target mu_units_test
mu_units_test_EXTERNAL_OBJECTS =

tests/mu_units_test: tests/CMakeFiles/mu_units_test.dir/units-mu-test.cpp.o
tests/mu_units_test: tests/CMakeFiles/mu_units_test.dir/build.make
tests/mu_units_test: libutil.so
tests/mu_units_test: tests/CMakeFiles/mu_units_test.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/wrk3/simploce/pt-cgmd/util/cmake-build-release/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable mu_units_test"
	cd /wrk3/simploce/pt-cgmd/util/cmake-build-release/tests && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/mu_units_test.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tests/CMakeFiles/mu_units_test.dir/build: tests/mu_units_test
.PHONY : tests/CMakeFiles/mu_units_test.dir/build

tests/CMakeFiles/mu_units_test.dir/clean:
	cd /wrk3/simploce/pt-cgmd/util/cmake-build-release/tests && $(CMAKE_COMMAND) -P CMakeFiles/mu_units_test.dir/cmake_clean.cmake
.PHONY : tests/CMakeFiles/mu_units_test.dir/clean

tests/CMakeFiles/mu_units_test.dir/depend:
	cd /wrk3/simploce/pt-cgmd/util/cmake-build-release && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /wrk3/simploce/pt-cgmd/util /wrk3/simploce/pt-cgmd/util/tests /wrk3/simploce/pt-cgmd/util/cmake-build-release /wrk3/simploce/pt-cgmd/util/cmake-build-release/tests /wrk3/simploce/pt-cgmd/util/cmake-build-release/tests/CMakeFiles/mu_units_test.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tests/CMakeFiles/mu_units_test.dir/depend

