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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/user/projects/CompetitivePrograming/Rosalind

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/user/projects/CompetitivePrograming/Rosalind/build

# Include any dependencies generated for this target.
include CMakeFiles/Rosalind.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/Rosalind.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/Rosalind.dir/flags.make

CMakeFiles/Rosalind.dir/Source.cpp.o: CMakeFiles/Rosalind.dir/flags.make
CMakeFiles/Rosalind.dir/Source.cpp.o: ../Source.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/user/projects/CompetitivePrograming/Rosalind/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/Rosalind.dir/Source.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Rosalind.dir/Source.cpp.o -c /home/user/projects/CompetitivePrograming/Rosalind/Source.cpp

CMakeFiles/Rosalind.dir/Source.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Rosalind.dir/Source.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/user/projects/CompetitivePrograming/Rosalind/Source.cpp > CMakeFiles/Rosalind.dir/Source.cpp.i

CMakeFiles/Rosalind.dir/Source.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Rosalind.dir/Source.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/user/projects/CompetitivePrograming/Rosalind/Source.cpp -o CMakeFiles/Rosalind.dir/Source.cpp.s

CMakeFiles/Rosalind.dir/BigInteger.cpp.o: CMakeFiles/Rosalind.dir/flags.make
CMakeFiles/Rosalind.dir/BigInteger.cpp.o: ../BigInteger.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/user/projects/CompetitivePrograming/Rosalind/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/Rosalind.dir/BigInteger.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Rosalind.dir/BigInteger.cpp.o -c /home/user/projects/CompetitivePrograming/Rosalind/BigInteger.cpp

CMakeFiles/Rosalind.dir/BigInteger.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Rosalind.dir/BigInteger.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/user/projects/CompetitivePrograming/Rosalind/BigInteger.cpp > CMakeFiles/Rosalind.dir/BigInteger.cpp.i

CMakeFiles/Rosalind.dir/BigInteger.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Rosalind.dir/BigInteger.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/user/projects/CompetitivePrograming/Rosalind/BigInteger.cpp -o CMakeFiles/Rosalind.dir/BigInteger.cpp.s

# Object files for target Rosalind
Rosalind_OBJECTS = \
"CMakeFiles/Rosalind.dir/Source.cpp.o" \
"CMakeFiles/Rosalind.dir/BigInteger.cpp.o"

# External object files for target Rosalind
Rosalind_EXTERNAL_OBJECTS =

Rosalind: CMakeFiles/Rosalind.dir/Source.cpp.o
Rosalind: CMakeFiles/Rosalind.dir/BigInteger.cpp.o
Rosalind: CMakeFiles/Rosalind.dir/build.make
Rosalind: CMakeFiles/Rosalind.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/user/projects/CompetitivePrograming/Rosalind/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable Rosalind"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/Rosalind.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/Rosalind.dir/build: Rosalind

.PHONY : CMakeFiles/Rosalind.dir/build

CMakeFiles/Rosalind.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/Rosalind.dir/cmake_clean.cmake
.PHONY : CMakeFiles/Rosalind.dir/clean

CMakeFiles/Rosalind.dir/depend:
	cd /home/user/projects/CompetitivePrograming/Rosalind/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/user/projects/CompetitivePrograming/Rosalind /home/user/projects/CompetitivePrograming/Rosalind /home/user/projects/CompetitivePrograming/Rosalind/build /home/user/projects/CompetitivePrograming/Rosalind/build /home/user/projects/CompetitivePrograming/Rosalind/build/CMakeFiles/Rosalind.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/Rosalind.dir/depend

