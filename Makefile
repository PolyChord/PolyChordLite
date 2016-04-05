# List of available example likelihoods
EXAMPLES = gaussian pyramidal rastrigin twin_gaussian random_gaussian himmelblau rosenbrock eggbox half_gaussian fitting gaussian_shell gaussian_shells
SIMPLE_EXAMPLES = polychord_simple polychord_simple_C

# Your likelihood programs
PROGRAMS = my_likelihood my_cpp_likelihood 

# Whether to use MPI
MPI=1

# Whether to compile in debugging mode
DEBUG=
export MPI DEBUG

# We can autodetect the compiler type on unix systems via the shell.
# if you want to override this then just run make with
# make COMPILER_TYPE=<your type>
# where <your time> is gnu or intel
ifeq "$(shell which ifort >/dev/null; echo $$?)" "0" 
COMPILER_TYPE=intel
else ifeq "$(shell which gfortran >/dev/null; echo $$?)" "0"
COMPILER_TYPE=gnu
endif

ifeq ($(COMPILER_TYPE),intel)
include Makefile_intel
else ifeq ($(COMPILER_TYPE),gnu) 
include Makefile_gnu
endif



ifdef MPI
FFLAGS += -DMPI
CXXFLAGS += -DMPI
endif

# Remove command
RM = rm -f

# Where polychord is stored
POLYCHORD_DIR = $(PWD)/src
# Where likelihoods are stored
LIKELIHOOD_DIR = $(PWD)/likelihoods
# Where likelihood examples are stored
EXAMPLES_DIR = $(PWD)/example_likelihoods 
# Where binaries are stored
BIN_DIR = $(PWD)/bin

# Library flags
LDFLAGS += -L$(POLYCHORD_DIR)
LDFLAGS += -L$(EXAMPLES_DIR)
LDFLAGS += -L$(LIKELIHOOD_DIR)
LDLIBS += -lchord 

INCFLAGS += -I$(POLYCHORD_DIR) 
INCFLAGS += -I$(EXAMPLES_DIR)
INCFLAGS += -I$(LIKELIHOOD_DIR)

# example likelihood libraries, this is created by changing X to libX.a
EXAMPLE_LIBRARIES = $(patsubst %,lib%.a,$(EXAMPLES))

# likelihood libraries, this is created by changing X to libX.a
PROGRAM_LIBRARIES = $(patsubst %,lib%.a,$(PROGRAMS))


# Export all of the necessary variables
export DEBUG MPI AR ARFLAGS FC CXX CXXFLAGS FFLAGS RM LD
export EXAMPLES EXAMPLE_LIBRARIES SIMPLE_EXAMPLES
export POLYCHORD_DIR


# "make" builds all
all: gaussian
#all: $(SIMPLE_EXAMPLES) $(EXAMPLES) $(PROGRAMS)

examples: $(EXAMPLES)

# Rule for building polychord static library
libchord.a:
	$(MAKE) -C $(POLYCHORD_DIR) libchord.a

libchord.so:
	$(MAKE) -C $(POLYCHORD_DIR) libchord.so

# Rule for building example likelihood libraries
$(EXAMPLE_LIBRARIES): libchord.a
	$(MAKE) -C $(EXAMPLES_DIR) $@

# Rule for building likelihood libraries
$(PROGRAM_LIBRARIES): libchord.a
	$(MAKE) -C $(LIKELIHOOD_DIR) $@

# Rule for example programs
$(EXAMPLES): %: libchord.a lib%.a polychord_examples.o
	$(LD) polychord_examples.o -o $(BIN_DIR)/$@ $(LDFLAGS) $(LDLIBS) -l$@

$(PROGRAMS): %: libchord.a lib%.a polychord.o 
	$(LD) polychord.o -o $(BIN_DIR)/$@ $(LDFLAGS) $(LDLIBS) -l$@

# Rule for simple example programs
$(SIMPLE_EXAMPLES): %: libchord.a %.o 
	$(LD) $@.o -o $(BIN_DIR)/$@ $(LDFLAGS) $(LDLIBS)

PyPolyChord: _PyPolyChord.o libchord.so
	$(LD) -shared _PyPolyChord.o -o _PyPolyChord.so  $(LDFLAGS) $(LDLIBS) 

_PyPolyChord.o: 
	$(CXX) $(CXXFLAGS) -Isrc/ -I/usr/include/python2.7 -c _PyPolyChord.cpp -o _PyPolyChord.o


# Rule for building fortran files
%.o: %.F90 
	$(FC) $(FFLAGS) $(INCFLAGS) -c $< 
%.o: %.f90 
	$(FC) $(FFLAGS) $(INCFLAGS) -c $< 
# Rule for building c++ files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCFLAGS) -c $< 

.PHONY: clean veryclean

clean:
	$(RM) *.o *.mod *.MOD
	$(MAKE) -C $(POLYCHORD_DIR) clean
	$(MAKE) -C $(EXAMPLES_DIR) clean
	$(MAKE) -C $(LIKELIHOOD_DIR) clean
	$(MAKE) -C $(BIN_DIR) clean
	
veryclean: clean
	$(RM) *~ 
	$(RM) _PyPolyChord.so
	$(MAKE) -C $(POLYCHORD_DIR) veryclean
	$(MAKE) -C $(EXAMPLES_DIR) veryclean
	$(MAKE) -C $(LIKELIHOOD_DIR) veryclean
	$(MAKE) -C $(BIN_DIR) veryclean
