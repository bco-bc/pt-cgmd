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
CMAKE_SOURCE_DIR = /wrk3/simploce/pt-cgmd/forcefields

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /wrk3/simploce/pt-cgmd/forcefields/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/forcefields.dir/depend.make
# Include the progress variables for this target.
include CMakeFiles/forcefields.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/forcefields.dir/flags.make

CMakeFiles/forcefields.dir/src/no_bc.cpp.o: CMakeFiles/forcefields.dir/flags.make
CMakeFiles/forcefields.dir/src/no_bc.cpp.o: ../src/no_bc.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/wrk3/simploce/pt-cgmd/forcefields/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/forcefields.dir/src/no_bc.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/forcefields.dir/src/no_bc.cpp.o -c /wrk3/simploce/pt-cgmd/forcefields/src/no_bc.cpp

CMakeFiles/forcefields.dir/src/no_bc.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/forcefields.dir/src/no_bc.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /wrk3/simploce/pt-cgmd/forcefields/src/no_bc.cpp > CMakeFiles/forcefields.dir/src/no_bc.cpp.i

CMakeFiles/forcefields.dir/src/no_bc.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/forcefields.dir/src/no_bc.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /wrk3/simploce/pt-cgmd/forcefields/src/no_bc.cpp -o CMakeFiles/forcefields.dir/src/no_bc.cpp.s

# Object files for target forcefields
forcefields_OBJECTS = \
"CMakeFiles/forcefields.dir/src/no_bc.cpp.o"

# External object files for target forcefields
forcefields_EXTERNAL_OBJECTS =

libforcefields.so: CMakeFiles/forcefields.dir/src/no_bc.cpp.o
libforcefields.so: CMakeFiles/forcefields.dir/build.make
libforcefields.so: CMakeFiles/forcefields.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/wrk3/simploce/pt-cgmd/forcefields/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX shared library libforcefields.so"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/forcefields.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/forcefields.dir/build: libforcefields.so
.PHONY : CMakeFiles/forcefields.dir/build

CMakeFiles/forcefields.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/forcefields.dir/cmake_clean.cmake
.PHONY : CMakeFiles/forcefields.dir/clean

CMakeFiles/forcefields.dir/depend:
	cd /wrk3/simploce/pt-cgmd/forcefields/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /wrk3/simploce/pt-cgmd/forcefields /wrk3/simploce/pt-cgmd/forcefields /wrk3/simploce/pt-cgmd/forcefields/cmake-build-debug /wrk3/simploce/pt-cgmd/forcefields/cmake-build-debug /wrk3/simploce/pt-cgmd/forcefields/cmake-build-debug/CMakeFiles/forcefields.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/forcefields.dir/depend
