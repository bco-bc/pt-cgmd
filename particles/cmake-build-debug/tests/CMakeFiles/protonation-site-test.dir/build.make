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
CMAKE_SOURCE_DIR = /wrk3/simploce/pt-cgmd/particles

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /wrk3/simploce/pt-cgmd/particles/cmake-build-debug

# Include any dependencies generated for this target.
include tests/CMakeFiles/protonation-site-test.dir/depend.make
# Include the progress variables for this target.
include tests/CMakeFiles/protonation-site-test.dir/progress.make

# Include the compile flags for this target's objects.
include tests/CMakeFiles/protonation-site-test.dir/flags.make

tests/CMakeFiles/protonation-site-test.dir/protonation-site-test.cpp.o: tests/CMakeFiles/protonation-site-test.dir/flags.make
tests/CMakeFiles/protonation-site-test.dir/protonation-site-test.cpp.o: ../tests/protonation-site-test.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/wrk3/simploce/pt-cgmd/particles/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object tests/CMakeFiles/protonation-site-test.dir/protonation-site-test.cpp.o"
	cd /wrk3/simploce/pt-cgmd/particles/cmake-build-debug/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/protonation-site-test.dir/protonation-site-test.cpp.o -c /wrk3/simploce/pt-cgmd/particles/tests/protonation-site-test.cpp

tests/CMakeFiles/protonation-site-test.dir/protonation-site-test.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/protonation-site-test.dir/protonation-site-test.cpp.i"
	cd /wrk3/simploce/pt-cgmd/particles/cmake-build-debug/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /wrk3/simploce/pt-cgmd/particles/tests/protonation-site-test.cpp > CMakeFiles/protonation-site-test.dir/protonation-site-test.cpp.i

tests/CMakeFiles/protonation-site-test.dir/protonation-site-test.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/protonation-site-test.dir/protonation-site-test.cpp.s"
	cd /wrk3/simploce/pt-cgmd/particles/cmake-build-debug/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /wrk3/simploce/pt-cgmd/particles/tests/protonation-site-test.cpp -o CMakeFiles/protonation-site-test.dir/protonation-site-test.cpp.s

# Object files for target protonation-site-test
protonation__site__test_OBJECTS = \
"CMakeFiles/protonation-site-test.dir/protonation-site-test.cpp.o"

# External object files for target protonation-site-test
protonation__site__test_EXTERNAL_OBJECTS =

tests/protonation-site-test: tests/CMakeFiles/protonation-site-test.dir/protonation-site-test.cpp.o
tests/protonation-site-test: tests/CMakeFiles/protonation-site-test.dir/build.make
tests/protonation-site-test: libparticles.so
tests/protonation-site-test: tests/CMakeFiles/protonation-site-test.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/wrk3/simploce/pt-cgmd/particles/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable protonation-site-test"
	cd /wrk3/simploce/pt-cgmd/particles/cmake-build-debug/tests && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/protonation-site-test.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tests/CMakeFiles/protonation-site-test.dir/build: tests/protonation-site-test
.PHONY : tests/CMakeFiles/protonation-site-test.dir/build

tests/CMakeFiles/protonation-site-test.dir/clean:
	cd /wrk3/simploce/pt-cgmd/particles/cmake-build-debug/tests && $(CMAKE_COMMAND) -P CMakeFiles/protonation-site-test.dir/cmake_clean.cmake
.PHONY : tests/CMakeFiles/protonation-site-test.dir/clean

tests/CMakeFiles/protonation-site-test.dir/depend:
	cd /wrk3/simploce/pt-cgmd/particles/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /wrk3/simploce/pt-cgmd/particles /wrk3/simploce/pt-cgmd/particles/tests /wrk3/simploce/pt-cgmd/particles/cmake-build-debug /wrk3/simploce/pt-cgmd/particles/cmake-build-debug/tests /wrk3/simploce/pt-cgmd/particles/cmake-build-debug/tests/CMakeFiles/protonation-site-test.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tests/CMakeFiles/protonation-site-test.dir/depend

