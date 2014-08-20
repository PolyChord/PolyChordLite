# A simple hand-made makefile for a package including applications
# built from Fortran 90 sources, taking into account the usual
# dependency cases.

# This makefile works with the GNU make command, the one find on
# GNU/Linux systems and often called gmake on non-GNU systems, if you
# are using an old style make command, please see the file
# Makefile_oldstyle provided with the package.

# ======================================================================
# Let's start with the declarations
# ======================================================================

# whether its linear or parallel
PAR ?= 1

# The compiler
ifeq ($(PAR), 1)
	FC = mpif90
else
	FC = ifort
endif
# flags for debugging or for maximum performance, comment as necessary
DEBUG ?= 0
ifeq ($(DEBUG), 1)
FCFLAGS = -g -traceback -fp-stack-check -check all -fpe0 -warn all -mkl -fpp
else
FCFLAGS = -O3 -mkl -fpp -xHost
# O3  : maximum optimisation
# mkl : use the mkl library
# fpp : perform preprocessing
endif

ifeq ($(PAR), 1)
	FCFLAGS += -DMPI
endif


# flags forall (e.g. look for system .mod files, required in gfortran)
#FCFLAGS += -I/usr/include

# libraries needed for linking, unused in the examples
LDFLAGS = #-lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lmkl_lapack -lmkl_blas95 -lmkl_vml

# List of executables to be built within the package
PROGRAMS = main

# "make" builds all
all: $(PROGRAMS)

# ======================================================================
# Here comes the most interesting part: the rules for prog1, prog2,
# prog3 and prog4, modify to suit your needs
# ======================================================================
ifeq ($(PAR), 1)
endif

ifeq ($(PAR), 1)
NESTED_SAMPLING = nested_sampling_parallel.o nested_sampling_linear.o 
MPI_UTILS = mpi_utils.o
else
NESTED_SAMPLING = nested_sampling_linear.o
endif

MODULES = $(MPI_UTILS) utils.o model.o settings.o feedback.o read_write.o random_utils.o evidence.o likelihoods.o chordal_sampling.o $(NESTED_SAMPLING)




$(MPI_UTILS):
utils.o:
model.o:                    utils.o
random_utils.o:             $(MPI_UTILS) utils.o
settings.o:                 model.o 
likelihoods.o:              $(MPI_UTILS) model.o utils.o 
evidence.o:                 settings.o
feedback.o:                 model.o settings.o utils.o
read_write.o:               model.o settings.o utils.o
chordal_sampling.o:         model.o utils.o random_utils.o
$(NESTED_SAMPLING):   		$(MPI_UTILS) read_write.o feedback.o model.o utils.o random_utils.o settings.o chordal_sampling.o
main.o:                     $(MODULES)

main: $(MODULES) main.o 


# ======================================================================
# And now the general rules, these should not require modification
# ======================================================================

# General rule for building prog from prog.o; $^ (GNU extension) is
# used in order to list additional object files on which the
# executable depends
%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

# General rules for building prog.o from prog.f90 or prog.F90; $< is
# used in order to list only the first prerequisite (the source file)
# and not the additional prerequisites such as module or include files
%.o: %.f90
	$(FC) $(FCFLAGS) -c $< $(LDFLAGS)

%.o: %.F90
	$(FC) $(FCFLAGS) -c $<

# Utility targets
.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD

veryclean: clean
	rm -f *~ $(PROGRAMS)


