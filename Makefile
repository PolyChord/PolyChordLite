# The compiler type:
# This sets up the correct compiler flags below 
# current options are:
#  * gfortran
#  * ifort
COMPILER_TYPE=ifort


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
# Choose the mpi wrapper
FC = mpif90
else
# No parallelisation
# ------------------
# Set the compiler to be that defined in COMPILER_TYPE
FC = $(COMPILER_TYPE)
endif




# ============ COMPILER_TYPE settings ================

# ifort settings
# ==============
ifeq ($(COMPILER_TYPE),ifort)
# default flags
# --------------
# fpp  : perform preprocessing
FCFLAGS += -fpp 

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
else
# Optimised mode
# --------------
#   ipo          : interprocedural optimization (optimize entire program)
#   O3           : maximum optimisation
#   no-prec-div  : slightly less precise floating point divides, but speeds up
#   static       : link intel libraries statically
#   xHost        : maximise architecture usage
FCFLAGS += -ipo -O3 -no-prec-div -xHost
endif

# Archive tool for compiling with ipo.
AR = xiar r

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
 
# Archive tool
AR = ar r
endif

# Export all of the variables
export DEBUG COMPILER_TYPE FCFLAGS MPI AR FC

# List of available executables
PROGRAMS = gaussian rastrigin

# likelihood libraries, this is created by changing X from PROGRAMS into libX.a
LIKELIHOODS = $(patsubst %,lib%.a,$(PROGRAMS))



# "make" builds all
all: $(PROGRAMS)


# Rule for building polychord library
libchord.a: ./src/*90
	cd ./src && make libchord.a

# Rule for building likelihood libraries
$(LIKELIHOODS): libchord.a
	cd ./likelihoods && make $@

# Rule for building each program
$(PROGRAMS): %: libchord.a lib%.a  polychord.o
	$(FC) $(FCFLAGS) -o bin/$@ polychord.o -L./src/ -L./likelihoods/ -lchord -l$@

# Rule for building main file
polychord.o: polychord.F90
	$(FC) $(FCFLAGS) -I./src/ -I./likelihoods/ -c $< 



clean:
	rm -f *.o *.mod *.MOD
	cd ./src && make clean && cd - && cd ./likelihoods && make clean
	
veryclean: clean
	rm -f *~ 
	cd ./src && make veryclean && cd -
	cd ./likelihoods && make veryclean && cd -
	cd ./bin && rm -f $(PROGRAMS)
