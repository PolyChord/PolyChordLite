# List of available example likelihoods
EXAMPLES = gaussian pyramidal rastrigin twin_gaussian random_gaussian himmelblau rosenbrock eggbox half_gaussian fitting gaussian_shell gaussian_shells
SIMPLE_EXAMPLES = polychord_simple polychord_simple_C

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
FFLAGS += -DMPI
CXXFLAGS += -DMPI
# Choose the mpi wrapper
FC = mpif90
CXX = mpicc
LD = mpif90
else
FC = gfortran
CXX = g++
LD = gfortran
endif


# Archive tool
AR = ar r


# default flags
# --------------
# free-line-length-none : turn of line length limitation (why is this not a default??)
# cpp  					: perform preprocessing
# fPIC                  : for compiling a shared object library
FFLAGS += -ffree-line-length-none -cpp -fPIC
CXXFLAGS += -std=c++11

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
FFLAGS += -g -O0 -Wall -Wextra -pedantic -fcheck=all -fimplicit-none -fbacktrace -ffpe-trap=zero,overflow 
#
CXXFLAGS += -g -O0 -Wall -Wextra -Wshadow -Weffc++
else
# Optimised mode
# --------------
# Ofast : maximum optimisation
FFLAGS += -Ofast
CXXFLAGS += -Ofast
endif

# Where polychord is stored
POLYCHORD_DIR = ./src
# Where likelihoods are stored
LIKELIHOOD_DIR = ./likelihoods
# Where likelihood examples are stored
EXAMPLES_DIR = ./example_likelihoods 
# Where binaries are stored
BIN_DIR = ./bin

# Library flags
LDFLAGS += -L$(POLYCHORD_DIR)
LDFLAGS += -L$(EXAMPLES_DIR)
LDLIBS += -lchord 
LDLIBS += -lstdc++

# likelihood libraries, this is created by changing X to libX.a
EXAMPLE_LIBRARIES = $(patsubst %,lib%.a,$(EXAMPLES))


# Export all of the necessary variables
export DEBUG MPI AR FC CXX CXXFLAGS FFLAGS EXAMPLE_LIBRARIES


# "make" builds all
all: $(SIMPLE_EXAMPLES)

examples: $(EXAMPLES)

# Rule for building polychord static library
libchord.a:
	$(MAKE) -C $(POLYCHORD_DIR) libchord.a

# Rule for building likelihood libraries
$(EXAMPLE_LIBRARIES): libchord.a
	$(MAKE) -C $(EXAMPLES_DIR) $@

# Rule for example programs
$(EXAMPLES): %: libchord.a lib%.a polychord_examples.o
	$(LD) polychord_examples.o -o $(BIN_DIR)/$@ $(LDFLAGS) $(LDLIBS) -l$@

# Rule for simple example programs
$(SIMPLE_EXAMPLES): %: libchord.a %.o 
	$(LD) $@.o -o $(BIN_DIR)/$@ $(LDFLAGS) $(LDLIBS)

# Rule for building fortran files
%.o: %.F90 
	$(FC) $(FFLAGS) -I$(POLYCHORD_DIR) -I$(EXAMPLES_DIR) -c $< 
%.o: %.f90 
	$(FC) $(FFLAGS) -I$(POLYCHORD_DIR) -I$(EXAMPLES_DIR) -c $< 
# Rule for building c++ files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -I$(POLYCHORD_DIR) -I$(EXAMPLES_DIR) -c $< 


clean:
	rm -f *.o *.mod *.MOD
	$(MAKE) -C $(POLYCHORD_DIR) clean
	$(MAKE) -C $(EXAMPLES_DIR) clean
	
veryclean: clean
	rm -f *~ 
	$(MAKE) -C $(POLYCHORD_DIR) veryclean
	$(MAKE) -C $(EXAMPLES_DIR) veryclean
