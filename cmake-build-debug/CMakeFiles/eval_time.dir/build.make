# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /home/horacio/.local/share/JetBrains/Toolbox/apps/CLion/ch-0/201.7223.86/bin/cmake/linux/bin/cmake

# The command to remove a file.
RM = /home/horacio/.local/share/JetBrains/Toolbox/apps/CLion/ch-0/201.7223.86/bin/cmake/linux/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/horacio/CLionProjects/valgraphcore

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/horacio/CLionProjects/valgraphcore/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/eval_time.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/eval_time.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/eval_time.dir/flags.make

CMakeFiles/eval_time.dir/eval_time.c.o: CMakeFiles/eval_time.dir/flags.make
CMakeFiles/eval_time.dir/eval_time.c.o: ../eval_time.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/horacio/CLionProjects/valgraphcore/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/eval_time.dir/eval_time.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/eval_time.dir/eval_time.c.o   -c /home/horacio/CLionProjects/valgraphcore/eval_time.c

CMakeFiles/eval_time.dir/eval_time.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/eval_time.dir/eval_time.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/horacio/CLionProjects/valgraphcore/eval_time.c > CMakeFiles/eval_time.dir/eval_time.c.i

CMakeFiles/eval_time.dir/eval_time.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/eval_time.dir/eval_time.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/horacio/CLionProjects/valgraphcore/eval_time.c -o CMakeFiles/eval_time.dir/eval_time.c.s

CMakeFiles/eval_time.dir/vgraph.c.o: CMakeFiles/eval_time.dir/flags.make
CMakeFiles/eval_time.dir/vgraph.c.o: ../vgraph.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/horacio/CLionProjects/valgraphcore/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object CMakeFiles/eval_time.dir/vgraph.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/eval_time.dir/vgraph.c.o   -c /home/horacio/CLionProjects/valgraphcore/vgraph.c

CMakeFiles/eval_time.dir/vgraph.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/eval_time.dir/vgraph.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/horacio/CLionProjects/valgraphcore/vgraph.c > CMakeFiles/eval_time.dir/vgraph.c.i

CMakeFiles/eval_time.dir/vgraph.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/eval_time.dir/vgraph.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/horacio/CLionProjects/valgraphcore/vgraph.c -o CMakeFiles/eval_time.dir/vgraph.c.s

# Object files for target eval_time
eval_time_OBJECTS = \
"CMakeFiles/eval_time.dir/eval_time.c.o" \
"CMakeFiles/eval_time.dir/vgraph.c.o"

# External object files for target eval_time
eval_time_EXTERNAL_OBJECTS =

eval_time: CMakeFiles/eval_time.dir/eval_time.c.o
eval_time: CMakeFiles/eval_time.dir/vgraph.c.o
eval_time: CMakeFiles/eval_time.dir/build.make
eval_time: CMakeFiles/eval_time.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/horacio/CLionProjects/valgraphcore/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking C executable eval_time"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/eval_time.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/eval_time.dir/build: eval_time

.PHONY : CMakeFiles/eval_time.dir/build

CMakeFiles/eval_time.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/eval_time.dir/cmake_clean.cmake
.PHONY : CMakeFiles/eval_time.dir/clean

CMakeFiles/eval_time.dir/depend:
	cd /home/horacio/CLionProjects/valgraphcore/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/horacio/CLionProjects/valgraphcore /home/horacio/CLionProjects/valgraphcore /home/horacio/CLionProjects/valgraphcore/cmake-build-debug /home/horacio/CLionProjects/valgraphcore/cmake-build-debug /home/horacio/CLionProjects/valgraphcore/cmake-build-debug/CMakeFiles/eval_time.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/eval_time.dir/depend

