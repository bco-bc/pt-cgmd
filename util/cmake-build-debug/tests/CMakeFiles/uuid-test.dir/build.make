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
CMAKE_COMMAND = /home/andre/clion-2021.2.3/bin/cmake/linux/bin/cmake

# The command to remove a file.
RM = /home/andre/clion-2021.2.3/bin/cmake/linux/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/andre/wrk3/simploce/pt-cgmd/util

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/andre/wrk3/simploce/pt-cgmd/util/cmake-build-debug

# Include any dependencies generated for this target.
include tests/CMakeFiles/uuid-test.dir/depend.make
# Include the progress variables for this target.
include tests/CMakeFiles/uuid-test.dir/progress.make

# Include the compile flags for this target's objects.
include tests/CMakeFiles/uuid-test.dir/flags.make

tests/CMakeFiles/uuid-test.dir/id-test.cpp.o: tests/CMakeFiles/uuid-test.dir/flags.make
tests/CMakeFiles/uuid-test.dir/id-test.cpp.o: ../tests/id-test.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/andre/wrk3/simploce/pt-cgmd/util/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object tests/CMakeFiles/uuid-test.dir/id-test.cpp.o"
	cd /home/andre/wrk3/simploce/pt-cgmd/util/cmake-build-debug/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/uuid-test.dir/id-test.cpp.o -c /home/andre/wrk3/simploce/pt-cgmd/util/tests/id-test.cpp

tests/CMakeFiles/uuid-test.dir/id-test.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/uuid-test.dir/id-test.cpp.i"
	cd /home/andre/wrk3/simploce/pt-cgmd/util/cmake-build-debug/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/andre/wrk3/simploce/pt-cgmd/util/tests/id-test.cpp > CMakeFiles/uuid-test.dir/id-test.cpp.i

tests/CMakeFiles/uuid-test.dir/id-test.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/uuid-test.dir/id-test.cpp.s"
	cd /home/andre/wrk3/simploce/pt-cgmd/util/cmake-build-debug/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/andre/wrk3/simploce/pt-cgmd/util/tests/id-test.cpp -o CMakeFiles/uuid-test.dir/id-test.cpp.s

# Object files for target uuid-test
uuid__test_OBJECTS = \
"CMakeFiles/uuid-test.dir/id-test.cpp.o"

# External object files for target uuid-test
uuid__test_EXTERNAL_OBJECTS =

tests/uuid-test: tests/CMakeFiles/uuid-test.dir/id-test.cpp.o
tests/uuid-test: tests/CMakeFiles/uuid-test.dir/build.make
tests/uuid-test: libutil.so
tests/uuid-test: tests/CMakeFiles/uuid-test.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/andre/wrk3/simploce/pt-cgmd/util/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable uuid-test"
	cd /home/andre/wrk3/simploce/pt-cgmd/util/cmake-build-debug/tests && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/uuid-test.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tests/CMakeFiles/uuid-test.dir/build: tests/uuid-test
.PHONY : tests/CMakeFiles/uuid-test.dir/build

tests/CMakeFiles/uuid-test.dir/clean:
	cd /home/andre/wrk3/simploce/pt-cgmd/util/cmake-build-debug/tests && $(CMAKE_COMMAND) -P CMakeFiles/uuid-test.dir/cmake_clean.cmake
.PHONY : tests/CMakeFiles/uuid-test.dir/clean

tests/CMakeFiles/uuid-test.dir/depend:
	cd /home/andre/wrk3/simploce/pt-cgmd/util/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/andre/wrk3/simploce/pt-cgmd/util /home/andre/wrk3/simploce/pt-cgmd/util/tests /home/andre/wrk3/simploce/pt-cgmd/util/cmake-build-debug /home/andre/wrk3/simploce/pt-cgmd/util/cmake-build-debug/tests /home/andre/wrk3/simploce/pt-cgmd/util/cmake-build-debug/tests/CMakeFiles/uuid-test.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tests/CMakeFiles/uuid-test.dir/depend
