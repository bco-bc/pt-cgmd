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
CMAKE_SOURCE_DIR = /home/andre/wrk3/simploce/pt-cgmd/apps

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/andre/wrk3/simploce/pt-cgmd/apps/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/s-distance.dir/depend.make
# Include the progress variables for this target.
include CMakeFiles/s-distance.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/s-distance.dir/flags.make

CMakeFiles/s-distance.dir/src/s-distance.cpp.o: CMakeFiles/s-distance.dir/flags.make
CMakeFiles/s-distance.dir/src/s-distance.cpp.o: ../src/s-distance.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/andre/wrk3/simploce/pt-cgmd/apps/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/s-distance.dir/src/s-distance.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/s-distance.dir/src/s-distance.cpp.o -c /home/andre/wrk3/simploce/pt-cgmd/apps/src/s-distance.cpp

CMakeFiles/s-distance.dir/src/s-distance.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/s-distance.dir/src/s-distance.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/andre/wrk3/simploce/pt-cgmd/apps/src/s-distance.cpp > CMakeFiles/s-distance.dir/src/s-distance.cpp.i

CMakeFiles/s-distance.dir/src/s-distance.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/s-distance.dir/src/s-distance.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/andre/wrk3/simploce/pt-cgmd/apps/src/s-distance.cpp -o CMakeFiles/s-distance.dir/src/s-distance.cpp.s

# Object files for target s-distance
s__distance_OBJECTS = \
"CMakeFiles/s-distance.dir/src/s-distance.cpp.o"

# External object files for target s-distance
s__distance_EXTERNAL_OBJECTS =

s-distance: CMakeFiles/s-distance.dir/src/s-distance.cpp.o
s-distance: CMakeFiles/s-distance.dir/build.make
s-distance: /usr/lib/x86_64-linux-gnu/libboost_program_options.so.1.71.0
s-distance: CMakeFiles/s-distance.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/andre/wrk3/simploce/pt-cgmd/apps/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable s-distance"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/s-distance.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/s-distance.dir/build: s-distance
.PHONY : CMakeFiles/s-distance.dir/build

CMakeFiles/s-distance.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/s-distance.dir/cmake_clean.cmake
.PHONY : CMakeFiles/s-distance.dir/clean

CMakeFiles/s-distance.dir/depend:
	cd /home/andre/wrk3/simploce/pt-cgmd/apps/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/andre/wrk3/simploce/pt-cgmd/apps /home/andre/wrk3/simploce/pt-cgmd/apps /home/andre/wrk3/simploce/pt-cgmd/apps/cmake-build-debug /home/andre/wrk3/simploce/pt-cgmd/apps/cmake-build-debug /home/andre/wrk3/simploce/pt-cgmd/apps/cmake-build-debug/CMakeFiles/s-distance.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/s-distance.dir/depend

