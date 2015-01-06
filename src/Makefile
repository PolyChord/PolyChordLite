# Whether to run in debugging mode
DEBUG=

# Whether to use MPI or not
MPI = 1

# The compiler:
#COMPILER_TYPE=ifort
COMPILER_TYPE=gfortran

ifdef MPI
# Using MPI
# ---------
# Set a define flag for MPI
FCFLAGS = -DMPI
# Choose the mpi wrapper
FC = mpifc -fc=$(COMPILER_TYPE)
# Define the files to be used 
OBJECTS =                    \
		  utils.o            \
		  mpi_utils.o        \
		  random_utils.o     \
		  priors.o grades.o  \
		  settings.o         \
		  calculate.o        \
		  evidence.o         \
		  feedback.o         \
		  read_write.o       \
		  clustering.o       \
		  likelihoods.o      \
		  chordal_sampling.o \
		  generate.o         \
		  nested_sampling.o
else

# No parallelisation
# ------------------
# Set the compiler to be that defined in COMPILER_TYPE
FC = $(COMPILER_TYPE)
# Define the files to be used (note no mpi_utils file)
OBJECTS =                    \
		  utils.o            \
		  random_utils.o     \
		  priors.o grades.o  \
		  settings.o         \
		  calculate.o        \
		  evidence.o         \
		  feedback.o         \
		  read_write.o       \
		  clustering.o       \
		  likelihoods.o      \
		  chordal_sampling.o \
		  generate.o         \
		  nested_sampling.o

endif






# intel settings
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
FCFLAGS += -g -O0 #-traceback -fp-stack-check -check all -fpe0 -warn all
else
# Optimised mode
# --------------
# O3         : maximum optimisation
# xHost      : maximise architecture usage
FCFLAGS += -O3 -xHost
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
FCFLAGS += -g -O0 -Wall -Wextra -pedantic -fimplicit-none -fbacktrace -ffpe-trap=zero,overflow -fbounds-check
else
# Optimised mode
# --------------
# O3 : maximum optimisation
FCFLAGS += -O3
endif
 
endif








# libraries needed for linking, unused in the examples
LDFLAGS =

# List of executables to be built within the package
PROGRAMS = libchord.a main

# "make" builds all
all: $(PROGRAMS)


AR = ar r



libchord.a: $(OBJECTS)
	$(AR) $@ $^

main: main.o
	$(FC) $(FCFLAGS) -o ../main main.o -L. -lchord

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


print-%  : ; @echo $* = $($*)



