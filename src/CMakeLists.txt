#Name project and specify source languages
project(mir_pp Fortran)

# Some compilers will work with earlier versions of CMake. (The Cray Compiler Environment requires CMake 2.8.5)
cmake_minimum_required(VERSION 2.8.5)

set (CMAKE_Fortran_COMPILER "ifort")
# set (CMAKE_Fortran_FLAGS_RELEASE "-f77rtl -O3 -xHost -ipo -fast")

#Print an error message on an attempt to build inside the source directory tree:
IF ("${CMAKE_CURRENT_SOURCE_DIR}" STREQUAL "${CMAKE_CURRENT_BINARY_DIR}")
  MESSAGE(FATAL_ERROR "ERROR! "
    "CMAKE_CURRENT_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}"
    " == CMAKE_CURRENT_BINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}"
    "\nThis archive does not support in-source builds:\n"
    "You must now delete the CMakeCache.txt file and the CMakeFiles/ directory under"
    "the source directory for the book code archive or you will not be able to configure correctly!"
    "\nYou must now run something like:\n"
    "  $ rm -r CMakeCache.txt CMakeFiles/"
    "\n"
    "Please create a different directory and configure ssdSource under that different directory such as:\n"
    "  $ mkdir MY_BUILD\n"
    "  $ cd MY_BUILD\n"
    "  $ cmake [OPTIONS] .."
    )
ENDIF()

#Specify all the subdirectories containing codes that produces executable files:
set(folders_to_build driver)
foreach(folder ${folders_to_build})
  add_subdirectory(${folder})
endforeach()

#Report untested Fortran compiler unless explicitly directed to build all examples.
#if (${ALL_EXAMPLES_ENABLED} OR "${CMAKE_Fortran_COMPILER_ID}" MATCHES     "VisualAge|XL|NAG|GNU|Intel|Cray|PGI" )
#else()
#  message("\n*** Untested Fortran compiler detected: ${CMAKE_Fortran_COMPILER_ID}.  ***\n")
#  message("\n*** Attempting to build all examples. Please report any failures to aditya90@stanford.edu***\n")
#endif()
