# Carlos del-Castillo-Negrete, March 23, 2015
# Makefile for glram and montecarlo and supporting modules

# Compiler Information and Flags
FC = gfortran
DEBUG = -Og -g -Wall -Wextra -pedantic -fimplicit-none -fbounds-check -fbacktrace

# LAPACK and BLAS library info
LIBDIR = /home/carlos/Documents/Yale/SPRING2015/SeniorProject/lapack-3.5.0
LAPACK = lapack
BLAS = refblas

# FAMOSO Dependencies and Modules
FAMOSO = projective_montecarlo.o glram.o montecarlo.o mc_grid.o mc_reps.o utilities.o test_functions.o 
MODULES = projective_montecarlo.mod glram.mod montecarlo.mod mc_grid.mod mc_reps.mod utilities.mod 

modules: $(MODULES)

%.mod: %.f90 %.o
	@true

projective_montecarlo.o: projective_montecarlo.f90 montecarlo.mod glram.mod
	$(FC) -c $<

montecarlo.o: montecarlo.f90 mc_grid.mod mc_reps.mod 
	$(FC) -c $<

glram.o: glram.f90 utilities.mod
	$(FC) glram.f90 -O -c -L${LIBDIR} -l${LAPACK} -l${BLAS}

mc_reps.o: mc_reps.f90 mc_grid.mod
	$(FC) -c $<

mc_grid.o:  mc_grid.f90 utilities.mod
	$(FC) -c $<

utilities.o: utilities.f90
	$(FC) -c $<

clean:
	rm -f *.o ${MODULES} 
