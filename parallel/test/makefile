# Carlos del-Castillo-Negrete
# April, 2019
# Makefile example for building programs that use the functions in the monte-carlo library 

# Compiler Information and Flags
MPI_COMP = mpif90 -no-wrap-margin 
MPI_DEBUG = -g -O0 -g -traceback -check all -check bounds -check uninit -ftrapuv -debug all 
COMP = gfortran
DEBUG = -Og -g -Wall -Wextra -pedantic -fimplicit-none -fbounds-check -fbacktrace
FC = $(MPI_COMP) $(MPI_DEBUG)

# Absolute Path of monte-carlo library. Replace with appropriate path
MC_PATH =  ..

# Names of test script executable files
EXECUTABLES = test_mc 

all: test_mc 

test_mc: test_functions.o test_mc.o 
	$(FC) -o $@ $^ $(MC_PATH)/*.o 

%.o: %.f90
	$(FC) -c $< -I$(MC_PATH) 

clean:
	rm -f *.o ${EXECUTABLES} 

