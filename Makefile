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
