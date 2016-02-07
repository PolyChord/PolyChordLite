# List of available example likelihoods
#EXAMPLES = gaussian pyramidal rastrigin twin_gaussian random_gaussian himmelblau rosenbrock eggbox half_gaussian fitting gaussian_shell gaussian_shells

# Your likelihood programs
#PROGRAMS = my_likelihood my_cpp_likelihood 

# Whether to use MPI
MPI=

# Whether to compile in debugging mode
DEBUG=1

# ============ MPI settings ================
ifdef MPI
# Using MPI
# ---------
# Set a define flag for MPI
FCFLAGS += -DMPI
CCFLAGS += -DMPI
# Choose the mpi wrapper
FC = mpif90
CC = mpicc
LINKER = mpif90
else
FC = gfortran
CC = g++
LINKER = gfortran
endif


# Archive tool
AR = ar r


# default flags
# --------------
# free-line-length-none : turn of line length limitation (why is this not a default??)
# cpp  					: perform preprocessing
# fPIC                  : for compiling a shared object library
FCFLAGS += -ffree-line-length-none -cpp -fPIC

ifdef DEBUG
# Debugging mode
# --------------
# g             : enable gnu debugger compatibility
# O0            : no optimisation
# Wall          : all warnings
# Wextra        : even more warnings
# pedantic      : check for language features not part of f95 standard
# implicit-none : specify no implicit typing
# backtrace     : produce backtrace of error
# fpe-trap      : search for floating point exceptions (dividing by zero etc)
# fbounds-check : check array indices
FCFLAGS += -g -O0 -Wall -Wextra -pedantic -fcheck=all -fimplicit-none -fbacktrace -ffpe-trap=zero,overflow 
#
CCFLAGS += -g -O0 -Wall -Wextra -Wshadow -Weffc++
else
# Optimised mode
# --------------
# Ofast : maximum optimisation
FCFLAGS += -Ofast
CCFLAGS += -Ofast
endif

# Where polychord is stored
POLYCHORD_DIR = ./src

# Where likelihoods are stored
LIKELIHOOD_DIR = ./likelihoods

# Where likelihood examples are stored
EXAMPLE_LIKELIHOOD_DIR = ./example_likelihoods

# Where binaries are stored
BIN_DIR = ./bin

# Add these libraries to the program
LIBCHORD = -L$(POLYCHORD_DIR) -lchord
LIBSTDCPP = -lstdc++

# likelihood libraries, this is created by changing X to libX.a
#EXAMPLE_LIKELIHOODS = $(patsubst %,lib%.a,$(EXAMPLES))
#PROGRAM_LIKELIHOODS = $(patsubst %,lib%.a,$(PROGRAMS))


# Export all of the necessary variables
export DEBUG COMPILER_TYPE FCFLAGS MPI AR FC CC CCFLAGS EXAMPLE_LIKELIHOODS  


# "make" builds all
all: cpp_polychord


# Rule for building polychord static library
libchord.a:
	$(MAKE) -C $(POLYCHORD_DIR) libchord.a

# Rule for example programs
demo_gaussian: libchord.a polychord.o 
	$(LINKER) polychord.o -o $(BIN_DIR)/$@ $(LIBCHORD) 

cpp_polychord: libchord.a cpp_polychord.o 
	$(LINKER) cpp_polychord.o -o $(BIN_DIR)/$@ $(LIBCHORD) $(LIBSTDCPP)

# Rule for building fortran files
%.o: %.F90 
	$(FC) $(FCFLAGS) -I$(POLYCHORD_DIR) -c $< 
%.o: %.f90 
	$(FC) $(FCFLAGS) -I$(POLYCHORD_DIR) -c $< 
%.o: %.cpp
	$(CC) $(CCFLAGS) -I$(POLYCHORD_DIR) -c $< 


clean:
	rm -f *.o *.mod *.MOD
	$(MAKE) -C $(POLYCHORD_DIR) clean
	
veryclean: clean
	rm -f *~ 
	$(MAKE) -C $(POLYCHORD_DIR) veryclean
