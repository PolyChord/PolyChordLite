# List of available example likelihoods
EXAMPLES = gaussian pyramidal rastrigin twin_gaussian random_gaussian himmelblau rosenbrock eggbox half_gaussian fitting gaussian_shell gaussian_shells object_detection

# Your likelihood programs
PROGRAMS = polychord_fortran polychord_CC

# Directories
SRC_DIR = $(PWD)/src
DRIVERS_DIR = $(SRC_DIR)/drivers
POLYCHORD_DIR = $(SRC_DIR)/polychord
LIKELIHOOD_DIR = $(PWD)/likelihoods
EXAMPLES_DIR = $(LIKELIHOOD_DIR)/examples
BIN_DIR = $(PWD)/bin
LIB_DIR = $(PWD)/lib
export DRIVERS_DIR POLYCHORD_DIR PYPOLYCHORD_DIR LIKELIHOOD_DIR EXAMPLES_DIR BIN_DIR LIB_DIR 

# Whether to use MPI
ifeq "$(shell uname)" "Linux"
MPI=1
else
MPI= 
endif

# Whether to compile in debugging mode
DEBUG=
export MPI DEBUG

# We can autodetect the compiler type on unix systems via the shell.
# if you want to override this then just run make with
# make COMPILER_TYPE=<your type>
# where <your type> is gnu or intel
ifeq "$(shell which ifort >/dev/null 2>&1; echo $$?)" "0" 
COMPILER_TYPE=intel
else ifeq "$(shell which gfortran >/dev/null 2>&1; echo $$?)" "0"
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
RM = rm -rf

# Library flags
LDFLAGS += -L$(LIB_DIR)
LDLIBS += -lchord 


# Export all of the necessary variables
export CC CXX FC LD LDSHARED RM AR 
export CFLAGS CXXFLAGS FFLAGS
export EXAMPLES PROGRAMS
export PYTHON


# make shortcuts
all: $(LIB_DIR)/libchord.a $(LIB_DIR)/libchord.so
examples: $(EXAMPLES)
$(EXAMPLES): % : $(BIN_DIR)/%
$(PROGRAMS): % : $(BIN_DIR)/%

# PolyChord
# ---------
# static library
$(LIB_DIR)/libchord.a:
	$(MAKE) -C $(POLYCHORD_DIR) $@
# shared library
$(LIB_DIR)/libchord.so:
	$(MAKE) -C $(POLYCHORD_DIR) $@

PyPolyChord: $(LIB_DIR)/libchord.so
	@echo '======================================================================================='
	@echo ' now run:                                                                              '
	@echo '                                                                                       '
	@echo '    python3 setup.py install --user                                                    '
	@echo '                                                                                       '
	@echo ' or                                                                                    '
	@echo '                                                                                       '
	@echo '    python setup.py install --user                                                     '
	@echo '                                                                                       '
	@echo ' (The shared object and python compilation steps                                       '
	@echo '  separated to avoid compiler clashes)                                                 '
	@echo '                                                                                       '
	@echo ' OSX users may need to specify their C compilers to not be clang, e.g:                 '
	@echo '                                                                                       '
	@echo '    CC=/usr/local/bin/gcc-6 CXX=/usr/local/bin/g++-6  python3 setup.py install --user  '
	@echo '======================================================================================='



# Examples
# --------
$(patsubst %,$(BIN_DIR)/%,$(EXAMPLES)): $(BIN_DIR)/%: $(LIB_DIR)/libchord.a $(LIB_DIR)/lib%.a $(DRIVERS_DIR)/polychord_examples.o
	$(LD) $(DRIVERS_DIR)/polychord_examples.o -o $@ $(LDFLAGS) $(LDLIBS) -l$*

$(patsubst %,$(LIB_DIR)/lib%.a,$(EXAMPLES)): $(LIB_DIR)/libchord.a
	$(MAKE) -C $(EXAMPLES_DIR) $@

$(DRIVERS_DIR)/polychord_examples.o:
	$(MAKE) -C $(DRIVERS_DIR) $@

# User Likelihoods
# ----------------
$(patsubst %,$(BIN_DIR)/%,$(PROGRAMS)): $(BIN_DIR)/polychord_% : $(LIB_DIR)/libchord.a $(LIB_DIR)/lib%_likelihood.a $(DRIVERS_DIR)/polychord_%.o 
	$(LD) $(DRIVERS_DIR)/polychord_$*.o  -o $@ $(LDFLAGS) -l$*_likelihood $(LDLIBS) 

$(patsubst polychord_%,$(LIB_DIR)/lib%_likelihood.a,$(PROGRAMS)): $(LIB_DIR)/lib%_likelihood.a: $(LIB_DIR)/libchord.a
	$(MAKE) -C $(LIKELIHOOD_DIR)/$* $@

$(patsubst %,$(DRIVERS_DIR)/%.o,$(PROGRAMS)):
	$(MAKE) -C $(DRIVERS_DIR) $@


CLEANDIRS = $(POLYCHORD_DIR) $(PYPOLYCHORD_DIR) $(LIKELIHOOD_DIR) $(BIN_DIR) $(LIB_DIR) $(DRIVERS_DIR) 
.PHONY: clean veryclean $(addsuffix clean,$(CLEANDIRS)) $(addsuffix veryclean,$(CLEANDIRS))

clean: $(addsuffix clean,$(CLEANDIRS))
	$(RM) *.o *.mod *.MOD
$(addsuffix clean,$(CLEANDIRS)): %clean: 
	$(MAKE) -C $* clean

veryclean: clean $(addsuffix veryclean,$(CLEANDIRS))  
	$(RM) *~ build dist PyPolyChord.egg-info
$(addsuffix veryclean,$(CLEANDIRS))  : %veryclean: 
	$(MAKE) -C $* veryclean
	
