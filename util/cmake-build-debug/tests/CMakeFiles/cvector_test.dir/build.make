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
CMAKE_BINARY_DIR = /wrk3/simploce/pt-cgmd/util/cmake-build-debug

# Include any dependencies generated for this target.
include tests/CMakeFiles/cvector_test.dir/depend.make
# Include the progress variables for this target.
include tests/CMakeFiles/cvector_test.dir/progress.make

# Include the compile flags for this target's objects.
include tests/CMakeFiles/cvector_test.dir/flags.make

tests/CMakeFiles/cvector_test.dir/cvector-test.cpp.o: tests/CMakeFiles/cvector_test.dir/flags.make
tests/CMakeFiles/cvector_test.dir/cvector-test.cpp.o: ../tests/cvector-test.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/wrk3/simploce/pt-cgmd/util/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object tests/CMakeFiles/cvector_test.dir/cvector-test.cpp.o"
	cd /wrk3/simploce/pt-cgmd/util/cmake-build-debug/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/cvector_test.dir/cvector-test.cpp.o -c /wrk3/simploce/pt-cgmd/util/tests/cvector-test.cpp

tests/CMakeFiles/cvector_test.dir/cvector-test.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/cvector_test.dir/cvector-test.cpp.i"
	cd /wrk3/simploce/pt-cgmd/util/cmake-build-debug/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /wrk3/simploce/pt-cgmd/util/tests/cvector-test.cpp > CMakeFiles/cvector_test.dir/cvector-test.cpp.i

tests/CMakeFiles/cvector_test.dir/cvector-test.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/cvector_test.dir/cvector-test.cpp.s"
	cd /wrk3/simploce/pt-cgmd/util/cmake-build-debug/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /wrk3/simploce/pt-cgmd/util/tests/cvector-test.cpp -o CMakeFiles/cvector_test.dir/cvector-test.cpp.s

# Object files for target cvector_test
cvector_test_OBJECTS = \
"CMakeFiles/cvector_test.dir/cvector-test.cpp.o"

# External object files for target cvector_test
cvector_test_EXTERNAL_OBJECTS =

tests/cvector_test: tests/CMakeFiles/cvector_test.dir/cvector-test.cpp.o
tests/cvector_test: tests/CMakeFiles/cvector_test.dir/build.make
tests/cvector_test: libutil.so
tests/cvector_test: tests/CMakeFiles/cvector_test.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/wrk3/simploce/pt-cgmd/util/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable cvector_test"
	cd /wrk3/simploce/pt-cgmd/util/cmake-build-debug/tests && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/cvector_test.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tests/CMakeFiles/cvector_test.dir/build: tests/cvector_test
.PHONY : tests/CMakeFiles/cvector_test.dir/build

tests/CMakeFiles/cvector_test.dir/clean:
	cd /wrk3/simploce/pt-cgmd/util/cmake-build-debug/tests && $(CMAKE_COMMAND) -P CMakeFiles/cvector_test.dir/cmake_clean.cmake
.PHONY : tests/CMakeFiles/cvector_test.dir/clean

tests/CMakeFiles/cvector_test.dir/depend:
	cd /wrk3/simploce/pt-cgmd/util/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /wrk3/simploce/pt-cgmd/util /wrk3/simploce/pt-cgmd/util/tests /wrk3/simploce/pt-cgmd/util/cmake-build-debug /wrk3/simploce/pt-cgmd/util/cmake-build-debug/tests /wrk3/simploce/pt-cgmd/util/cmake-build-debug/tests/CMakeFiles/cvector_test.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tests/CMakeFiles/cvector_test.dir/depend

