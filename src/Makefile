# Whether to run in debugging mode
DEBUG ?= 0

# The compiler
FC = mpif90

# flags
# mkl  : use the mkl library
# fpp  : perform preprocessing
FCFLAGS = -mkl -fpp

ifeq ($(DEBUG), 1)
# Debugging mode
# --------------
# g              : enable gnu debugger compatibility
# traceback      : create a backtrace on system failure
# fp-stack-check : check for floating point errors from stack overflow, etc
# check all      : all checks (whilst compiling)
# fpe0           : halts program if dividing by zero, etc
# warn all       : all warnings (whilst running)
FCFLAGS += -g -traceback -fp-stack-check -check all -fpe0 -warn all
else
# Optimised mode
# --------------
# O3         : maximum optimisation
# xHost      : maximise architecture usage
FCFLAGS += -O3 -xHost
endif



# libraries needed for linking, unused in the examples
LDFLAGS =

# List of executables to be built within the package
PROGRAMS = libchord.a main planck batch

# "make" builds all
all: $(PROGRAMS)

OBJECTS = utils.o priors.o grades.o settings.o calculate.o feedback.o read_write.o random_utils.o clustering.o evidence.o likelihoods.o chordal_sampling.o nested_sampling.o

AR = ar r



libchord.a: $(OBJECTS)
	$(AR) $@ $^

main: main.o
	$(FC) $(FCFLAGS) -o ../main main.o -L. -lchord

batch: batch.o
	$(FC) $(FCFLAGS) -o ../batch batch.o -L. -lchord

planck: planck.o
	$(FC) $(FCFLAGS) -o ../planck planck.o -L. -lchord


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


