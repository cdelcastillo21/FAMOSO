# FAMOSO
A Fast Monte-Carlo Solver

# FAMOSO - a FAst MOnte carlo SOlver - v1.0.0
developed by Carlos del-Castillo-Negrete
Yale University - Senior Project for Computer Science Department
Spring 2015


## INSTALLING 

To install the modules contained in FAMOSO:

1) Obtain a copy of the BLAS and LAPACK libraries. These can be obtained from http://www.netlib.org/blas/ and http://www.netlib.org/lapack/ respectively. Once these have been downloaded and installed. 

2) Edit the file modules/makefile to change the path of the LAPACK library that is set to the path of where you downloaded your LAPACK library.

3) Inside the modules directory, type 

$ make modules

to make all the modules in FAMOSO.


## TEST SCRIPTS 

Contained in the directory test/ you will find a couple of test scripts that demonstrate how to use the modules in FAMOSO. The makefile contained in this directory demonstrates how to compile and link a program with the FAMOSO and BLAS/LAPACK libraries to run the Monte Carlo solver routines in FAMOSO. 


## PLOTTING DATA

In the octave/ director you will find a series of scripts used to plot and visualize data outputed by the FAMOSO Monte Carlo solver routines. The test scripts in the test/ directory demonstrate how to ouptut Monte Carlo simulation data (which is done by default to the data/ directory). These octave scripts are used to read in and plot this data. 

## PARALLEL VERSION

A parallel implementation (using MPI) of the FAMOSO library's Monte-Carlo integration step can be found in the parallel folder.
Scaling tests were done on the Texas Advanced Computing Center.
The corresponding report (and a presentation) can be found in the reports folder.


## MORE INFO

This project was completed back in 2015 to satisfy the senior undergraduate thesis requirements for Carlos del-Castillo-Negrete. 
The parallel work was completed in 2019 as part of a class project.
For more info or questions, contact at cdelcastillo21@gmail.com.
