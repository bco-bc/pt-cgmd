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
CMAKE_SOURCE_DIR = /wrk3/simploce/pt-cgmd/apps

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /wrk3/simploce/pt-cgmd/apps/cmake-build-release

# Include any dependencies generated for this target.
include CMakeFiles/s-validate.dir/depend.make
# Include the progress variables for this target.
include CMakeFiles/s-validate.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/s-validate.dir/flags.make

CMakeFiles/s-validate.dir/src/s-validate.cpp.o: CMakeFiles/s-validate.dir/flags.make
CMakeFiles/s-validate.dir/src/s-validate.cpp.o: ../src/s-validate.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/wrk3/simploce/pt-cgmd/apps/cmake-build-release/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/s-validate.dir/src/s-validate.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/s-validate.dir/src/s-validate.cpp.o -c /wrk3/simploce/pt-cgmd/apps/src/s-validate.cpp

CMakeFiles/s-validate.dir/src/s-validate.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/s-validate.dir/src/s-validate.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /wrk3/simploce/pt-cgmd/apps/src/s-validate.cpp > CMakeFiles/s-validate.dir/src/s-validate.cpp.i

CMakeFiles/s-validate.dir/src/s-validate.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/s-validate.dir/src/s-validate.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /wrk3/simploce/pt-cgmd/apps/src/s-validate.cpp -o CMakeFiles/s-validate.dir/src/s-validate.cpp.s

# Object files for target s-validate
s__validate_OBJECTS = \
"CMakeFiles/s-validate.dir/src/s-validate.cpp.o"

# External object files for target s-validate
s__validate_EXTERNAL_OBJECTS =

s-validate: CMakeFiles/s-validate.dir/src/s-validate.cpp.o
s-validate: CMakeFiles/s-validate.dir/build.make
s-validate: /usr/lib/x86_64-linux-gnu/libboost_program_options.so.1.71.0
s-validate: CMakeFiles/s-validate.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/wrk3/simploce/pt-cgmd/apps/cmake-build-release/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable s-validate"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/s-validate.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/s-validate.dir/build: s-validate
.PHONY : CMakeFiles/s-validate.dir/build

CMakeFiles/s-validate.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/s-validate.dir/cmake_clean.cmake
.PHONY : CMakeFiles/s-validate.dir/clean

CMakeFiles/s-validate.dir/depend:
	cd /wrk3/simploce/pt-cgmd/apps/cmake-build-release && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /wrk3/simploce/pt-cgmd/apps /wrk3/simploce/pt-cgmd/apps /wrk3/simploce/pt-cgmd/apps/cmake-build-release /wrk3/simploce/pt-cgmd/apps/cmake-build-release /wrk3/simploce/pt-cgmd/apps/cmake-build-release/CMakeFiles/s-validate.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/s-validate.dir/depend

