# Carlos del-Castillo-Negrete
# April, 2015
# Makefile example for building programs that use the functions in the FAMOSO library 

# Compiler Information and Flags
FC = gfortran
DEBUG = -Og -g -Wall -Wextra -pedantic -fimplicit-none -fbounds-check -fbacktrace

# LAPACK and BLAS library info. Replace with path to local copy of lapack library
LIBDIR = /home/carlos/Documents/Yale/SPRING2015/SeniorProject/lapack-3.5.0
LAPACK = lapack
BLAS = refblas

# Absolute Path of FAMOSO library. Replace with appropriate path
FAMOSO_PATH =  /home/carlos/Documents/Yale/SPRING2015/SeniorProject/FAMOSO/modules

# Names of test script executable files
EXECUTABLES = test_pjmc test_mc test_glram

all: test_pjmc test_mc test_glram
	

# Test scripts
test_pjmc: test_functions.o test_pjmc.o 
	$(FC) -o $@ $^ $(FAMOSO_PATH)/*.o -L${LIBDIR} -l${LAPACK} -l${BLAS}

test_mc: test_functions.o test_mc.o 
	$(FC) -o $@ $^ $(FAMOSO_PATH)/*.o -L${LIBDIR} -l${LAPACK} -l${BLAS}

test_glram: test_functions.o test_glram.o 
	$(FC) -o $@ $^ $(FAMOSO_PATH)/*.o -L${LIBDIR} -l${LAPACK} -l${BLAS}

%.o: %.f90
	$(FC) -c $< -I$(FAMOSO_PATH) 

clean:
	rm -f *.o ${EXECUTABLES} 

