# List of available example likelihoods
EXAMPLES = gaussian pyramidal rastrigin twin_gaussian random_gaussian himmelblau rosenbrock eggbox half_gaussian fitting gaussian_shell gaussian_shells

# Your likelihood programs
PROGRAMS = fortran_likelihood CC_likelihood 

# Directories
SRC_DIR = $(PWD)/src
DRIVERS_DIR = $(SRC_DIR)/drivers
POLYCHORD_DIR = $(SRC_DIR)/polychord
C_INTERFACE_DIR = $(SRC_DIR)/C_interface
PYPOLYCHORD_DIR = $(PWD)/PyPolyChord
LIKELIHOOD_DIR = $(PWD)/likelihoods
EXAMPLES_DIR = $(LIKELIHOOD_DIR)/examples
BIN_DIR = $(PWD)/bin
LIB_DIR = $(PWD)/lib
export DRIVERS_DIR POLYCHORD_DIR PYPOLYCHORD_DIR C_INTERFACE_DIR LIKELIHOOD_DIR EXAMPLES_DIR BIN_DIR LIB_DIR 


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

# Library flags
LDFLAGS += -L$(LIB_DIR)
LDLIBS += -lchord 

# example likelihood libraries, this is created by changing X to libX.a
EXAMPLE_LIBRARIES = $(patsubst %,$(LIB_DIR)/lib%.a,$(EXAMPLES))

# likelihood libraries, this is created by changing X to libX.a
PROGRAM_LIBRARIES = $(patsubst %,$(LIB_DIR)/lib%.a,$(PROGRAMS))

# Export all of the necessary variables
export CC CXX FC LD RM AR
export CFLAGS CXXFLAGS FFLAGS
export EXAMPLES EXAMPLE_LIBRARIES SIMPLE_EXAMPLES PROGRAMS


# make shortcuts
all: gaussian
examples: $(EXAMPLES)

.PHONY: $(EXAMPLES)

# PolyChord
# ---------
# static library
$(LIB_DIR)/libchord.a:
	$(MAKE) -C $(POLYCHORD_DIR) $(LIB_DIR)/libchord.a
# shared library
$(LIB_DIR)/libchord.so:
	$(MAKE) -C $(POLYCHORD_DIR) $(LIB_DIR)/libchord.so

# Examples
# --------
$(EXAMPLES): %: $(LIB_DIR)/libchord.a $(LIB_DIR)/lib%.a $(DRIVERS_DIR)/polychord_examples.o
	$(LD) $(DRIVERS_DIR)/polychord_examples.o -o $(BIN_DIR)/$@ $(LDFLAGS) $(LDLIBS) -l$@

$(EXAMPLE_LIBRARIES): $(LIB_DIR)/libchord.a
	$(MAKE) -C $(EXAMPLES_DIR) $@

$(DRIVERS_DIR)/polychord_examples.o:
	$(MAKE) -C $(DRIVERS_DIR) $@

# User Likelihoods
# ----------------
fortran_likelihood CC_likelihood: %_likelihood:  $(LIB_DIR)/libchord.a $(LIB_DIR)/lib%_likelihood.a %_polychord
	$(LD) $(DRIVERS_DIR)/polychord.o -o $(BIN_DIR)/$@ $(LDFLAGS) $(LDLIBS) -l$@

$(LIB_DIR)/libfortran_likelihood.a: $(LIB_DIR)/libchord.a
	$(MAKE) -C $(LIKELIHOOD_DIR)/fortran $@

$(LIB_DIR)/libCC_likelihood.a: $(LIB_DIR)/libchord.a
	$(MAKE) -C $(LIKELIHOOD_DIR)/CC $@

.PHONY: fortran_polychord CC_polychord
fortran_polychord CC_polychord:
	$(MAKE) -C $(DRIVERS_DIR) $@


# PyPolyChord
# -----------
PyPolyChord: $(LIB_DIR)/libchord.so
	$(MAKE) -C $(PYPOLYCHORD_DIR) PyPolyChord


.PHONY: clean veryclean

clean:
	$(RM) *.o *.mod *.MOD
	$(MAKE) -C $(POLYCHORD_DIR) clean
	$(MAKE) -C $(PYPOLYCHORD_DIR) clean
	$(MAKE) -C $(EXAMPLES_DIR) clean
	$(MAKE) -C $(LIKELIHOOD_DIR) clean
	$(MAKE) -C $(BIN_DIR) clean
	$(MAKE) -C $(LIB_DIR) clean
	$(MAKE) -C $(DRIVERS_DIR) clean
	
veryclean: clean
	$(RM) *~ 
	$(MAKE) -C $(POLYCHORD_DIR) veryclean
	$(MAKE) -C $(PYPOLYCHORD_DIR) veryclean
	$(MAKE) -C $(EXAMPLES_DIR) veryclean
	$(MAKE) -C $(LIKELIHOOD_DIR) veryclean
	$(MAKE) -C $(BIN_DIR) veryclean
	$(MAKE) -C $(LIB_DIR) veryclean
	$(MAKE) -C $(DRIVERS_DIR) veryclean
