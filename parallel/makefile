# Carlos del-Castillo-Negrete, April 2019
# Makefile for montecarlo and supporting modules

# Compiler Information and Flags
MPI_COMP = mpif90 -no-wrap-margin 
MPI_DEBUG = -g -O0 -g -traceback -check all -check bounds -check uninit -ftrapuv -debug all 
COMP = gfortran
DEBUG = -Og -g -Wall -Wextra -pedantic -fimplicit-none -fbounds-check -fbacktrace
FC = $(MPI_COMP) $(MPI_DEBUG)

# FAMOSO Dependencies and Modules
MODULES = montecarlo.mod mc_grid.mod mc_reps.mod utilities.mod 

modules: $(MODULES)

%.mod: %.f90 %.o
	@true

montecarlo.o: montecarlo.f90 mc_grid.mod mc_reps.mod 
	$(FC) -c $<

mc_reps.o: mc_reps.f90 mc_grid.mod
	$(FC) -c $<

mc_grid.o:  mc_grid.f90 utilities.mod
	$(FC) -c $<

utilities.o: utilities.f90
	$(FC) -c $<

clean:
	rm -f *.o ${MODULES} 
