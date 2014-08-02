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

# The compiler
FC = mpif90
# flags for debugging or for maximum performance, comment as necessary
DEBUG ?= 0
ifeq ($(DEBUG), 1)
FCFLAGS = -g -traceback -fp-stack-check -check all -fpe0 -warn all -mkl
else
FCFLAGS = -O3 -mkl
endif
# flags forall (e.g. look for system .mod files, required in gfortran)
#FCFLAGS += -I/usr/include

# libraries needed for linking, unused in the examples
#LDFLAGS = -li_need_this_lib

# List of executables to be built within the package
PROGRAMS = main

# "make" builds all
all: $(PROGRAMS)

# ======================================================================
# Here comes the most interesting part: the rules for prog1, prog2,
# prog3 and prog4, modify to suit your needs
# ======================================================================

MODULES = nested_sampling.o settings.o model.o random_utils.o likelihoods.o galilean_sample.o test_sample.o evidence.o feedback.o

feedback.o:         model.o settings.o
settings.o:			model.o
evidence.o:         model.o settings.o 
galilean_sample.o:  model.o settings.o random_utils.o
likelihoods.o:      model.o settings.o 
nested_sampling.o:  model.o settings.o random_utils.o galilean_sample.o feedback.o
main.o:   		  	$(MODULES)

main: main.o $(MODULES)


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
	$(FC) $(FCFLAGS) -c $<

%.o: %.F90
	$(FC) $(FCFLAGS) -c $<

# Utility targets
.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD

veryclean: clean
	rm -f *~ $(PROGRAMS)


