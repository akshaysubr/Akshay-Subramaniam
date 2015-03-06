# README #

This README would normally document whatever steps are necessary to get your application up and running.

### What is this repository for? ###

* A tool to postprocess Miranda visualization files

* Version 0.1

### How do I get set up? ###

* Summary of set up
    This project uses the CMake build system. You need to have CMake 2.8 or higher for the build process to work.
* Configuration
    Open the CMakeLists.txt file in the ./src directory and set the parameters for the MPI fortran compiler and compiler flags.
    Then, move to the ./build directory using 'mv ./build'
    Delete all the files here using 'rm -rf *'
    Run the cmake configuration using 'cmake ../src/'
    Then run 'make' to build the code
* Dependencies
    The only external dependency for this code is an mpif90 compiler
* Database configuration
* How to run tests
* Deployment instructions
    To run the code from the ./build/prob directory, run mpirun -np <# of procs> ./<program name> ../../input.txt
    Look at the sample input file and create a similar input file for your case

### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###

* Author: Akshay Subramaniam
* Email: akshays@stanford.edu
