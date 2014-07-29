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
FCFLAGS = -g -traceback -fp-stack-check -check all -fpe0 -warn all -mkl
#FCFLAGS = -O3 -mkl
# flags forall (e.g. look for system .mod files, required in gfortran)
#FCFLAGS += -I/usr/include

# libraries needed for linking, unused in the examples
#LDFLAGS = -li_need_this_lib

# List of executables to be built within the package
PROGRAMS = galileo

# "make" builds all
all: $(PROGRAMS)

# ======================================================================
# Here comes the most interesting part: the rules for prog1, prog2,
# prog3 and prog4, modify to suit your needs
# ======================================================================

MODULES = nested_sampling.o types.o random_utils.o likelihoods.o galilean_sample.o

galilean_sample.o:  random_utils.o types.o
likelihoods.o:      types.o
nested_sampling.o:  types.o random_utils.o galilean_sample.o
galileo.o: 			$(MODULES)

galileo: galileo.o $(MODULES)


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


