# List of available example likelihoods
EXAMPLES = gaussian rastrigin

PROGRAMS = my_likelihood my_cpp_likelihood 



# The compiler type:
# This sets up the correct compiler flags below 
# current options are:
#  * gnu
#  * intel
COMPILER_TYPE=intel


# Whether to use MPI or not (leave blank if running without MPI)
MPI=1

# Whether to run in debugging mode
DEBUG=


# ============ MPI settings ================
ifdef MPI
# Using MPI
# ---------
# Set a define flag for MPI
FCFLAGS = -DMPI
CCFLAGS = -DMPI
# Choose the mpi wrapper
FC = mpif90
CC = mpicc
else
# No parallelisation
# ------------------
# Set the compiler to be that defined in COMPILER_TYPE
ifeq ($(COMPILER_TYPE),intel)
FC = ifort
CC = icc
endif

ifeq ($(COMPILER_TYPE),gnu)
FC = gfortran
CC = g++
endif

endif



FCFLAGS += -lstdc++

# ============ COMPILER_TYPE settings ================
 
# Archive tool
AR = ar r

# ifort settings
# ==============
ifeq ($(COMPILER_TYPE),intel)
# default flags
# --------------
# fpp  : perform preprocessing
FCFLAGS += -fpp 
CCFLAGS += -fpp 

ifdef DEBUG
# Debugging mode
# --------------
# g              : enable gnu debugger compatibility
# O0             : no optimisation
# traceback      : create a backtrace on system failure
# fp-stack-check : check for floating point errors from stack overflow, etc
# check all      : all checks (whilst compiling)
# fpe0           : halts program if dividing by zero, etc
# warn all       : all warnings (whilst running)
FCFLAGS += -g -O0 -traceback -fp-stack-check -check all,noarg_temp_created -fpe0 -warn all
CCFLAGS += -g -O0 
else
# Optimised mode
# --------------
#   ipo          : interprocedural optimization (optimize entire program)
#   O3           : maximum optimisation
#   no-prec-div  : slightly less precise floating point divides, but speeds up
#   static       : link intel libraries statically
#   xHost        : maximise architecture usage
FCFLAGS += -ipo -O3 -no-prec-div -xHost
CCFLAGS += -ipo -O3 -no-prec-div -xHost
# Archive tool for compiling with ipo.
AR = xiar r
endif

endif



# gfortran settings
# =================
ifeq ($(COMPILER_TYPE),gfortran) 
# default flags
# --------------
# free-line-length-none : turn of line length limitation (why is this not a default??)
# cpp  					: perform preprocessing
FCFLAGS += -ffree-line-length-none -cpp 

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
else
# Optimised mode
# --------------
# O3 : maximum optimisation
FCFLAGS += -O3
endif
endif

# Where polychord is stored
POLYCHORD_DIR = ./src

# Where likelihoods are stored
LIKELIHOOD_DIR = ./likelihoods

# Where binaries are stored
BIN_DIR = ./bin

# Add these libraries to the program
LIBCHORD += -L$(POLYCHORD_DIR) -lchord
LIBLIKE += -L$(LIKELIHOOD_DIR)


INC += -I$(POLYCHORD_DIR) -I$(LIKELIHOOD_DIR)

# likelihood libraries, this is created by changing X to libX.a
EXAMPLE_LIKELIHOODS = $(patsubst %,lib%.a,$(EXAMPLES))
PROGRAM_LIKELIHOODS = $(patsubst %,lib%.a,$(PROGRAMS))


# Export all of the necessary variables
export DEBUG COMPILER_TYPE FCFLAGS MPI AR FC CC CCFLAGS EXAMPLE_LIKELIHOODS


# "make" builds all
all: $(EXAMPLES) $(PROGRAMS)


# Rule for building polychord library
libchord.a:
	cd $(POLYCHORD_DIR) && make libchord.a

# Rule for building likelihood libraries
$(EXAMPLE_LIKELIHOODS) $(PROGRAM_LIKELIHOODS): libchord.a
	cd $(LIKELIHOOD_DIR) && make $@

# Rule for example programs
$(EXAMPLES) $(PROGRAMS): %: libchord.a lib%.a  polychord.o
	$(FC) $(FCFLAGS) -o $(BIN_DIR)/$@ polychord.o $(LIBCHORD) $(LIBLIKE) -l$@

# Rule for building main file
polychord.o: polychord.F90
	$(FC) $(FCFLAGS) $(INC) -c $< 



clean:
	rm -f *.o *.mod *.MOD
	cd $(LIKELIHOOD_DIR) && make clean 
	cd $(POLYCHORD_DIR) && make clean
	
veryclean: clean
	rm -f *~ 
	cd $(POLYCHORD_DIR) && make veryclean
	cd $(LIKELIHOOD_DIR) && make veryclean
	cd $(BIN_DIR) && rm -f $(EXAMPLES) $(PROGRAMS)
