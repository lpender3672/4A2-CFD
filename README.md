## 4A2 - Computational Fluid Dynamics
Example and complete code for the IIB module

## Getting started
This directory contains all Fortran source files of the skeleton code.
"solver.f90" is the main program, each of the other .f90 files contains one
subroutine. It also contains the Python files for creating the test cases and
analysing the output of the solver.

## Example Fortran program  
This repository also contains an example program to calculate the advection of a
scalar variable "phi". Study all of the lines in "advection.f90" to familiarise
yourself with Fortran. The comments contain instructions for compiling and
executing this program. The output can then be plotted with the Python script
"plot\_advection.py". If you want to refresh your Python knowledge this is also
a good script to study line by line, there are only 11 of them.

## Completing the Euler solver
You need supply the correct lines in the .f90 files where indicated by "INSERT" 
to complete a workable program. The detailed instructions on how to do this are
given in the module handout available on moodle. The lecture notes and guides to
programming Fortran may also be useful.

## Compiling your Fortran program
The compilation of the code is managed by "makefile". Type 
```
make
```
the operating system will compile all the .f90 files into .o files and link the
.o files into an executable "solver.x". 

"make" is a unix bash system command which uses the information in "makefile" to
manage the compilation of executables.
 
There are many benefits of using a makefile. One is that the system only
re-compiles the newly updated .f90 files compared to their compiled .o versions.
That is, if a .f90 file is older than its .o file, the system will do nothing.
This is extremely efficient when dealing with large quantities of code. You may 
like to have a look at the makefile as you'll need to make a couple of changes 
later.

To force a recompile all files, delete all .o files with
```
make clean
```

## Creating a test case
The Python script "generate\_cases.py" can be used to create many different
geometries to model using your CFD solver. For the basic version of the code the
"bend" and "bump" cases form a good starting point. Move to the directory for 
your cases alongside your code directory and then execute the Python script.
```
cd ~/4A2/Cases
python ../Code/generate_case.py bend
```
it will show you a plot of the geometry and generate two input files
"geom\_bend.txt" and "input\_bend.txt"

## Running your code
Move to the directory with your input files in a terminal tab and execute the
solver. To run the code, type
```
../Code/solver.x < input_bend.txt
```
it will begin solving and print the status of the calculation in your terminal.
If you want to run the solver in the background and save this status to a log
file you should run
```
../Code/solver.x < input_bend.txt > log_bend.txt &
```

## Maintaining your code
If you haven't already, now is a very good time to learn a version control 
system. "Git" is the single most popular method and the university hosted
[Gitlab](https://gitlab.developers.cam.ac.uk/) has many different 
[tutorials](https://docs.gitlab.com/ee/tutorials/learn_git.html) for beginners 
and experienced users.

Alternatively there are some commands written into the makefile. Using
``` 
make save
```
will save your source files into a single compressed file "SaveSrc.tar.gz"
which is easy to transfer to other locations. To partially reverse the process 
type
```
make extract
```
Which restores the original files to a new directory "SaveSrc".

## Analysing your solutions
Four Python scripts are provided to analyse the output of your CFD solver.
"plot\_coord.py" will be the first you can use to plot the mesh created in
Fortran. You can run it on the bend test case with:
```
python ../Code/plot_coord.py bend"
```
Other test cases can be plotted by changing the casename after the script is
called. "plot\_guess.py" is the next that you can use to visualise the initial
guess, it is called in a similar way. Later "plot\_conv.py" can be used to show 
the convergence history, while "plot\_contours.py" shows the solution after a
successful run. This script has some blank sections for you to complete before
it will work.


