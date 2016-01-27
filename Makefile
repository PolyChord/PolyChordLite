# List of available example likelihoods
EXAMPLES = gaussian pyramidal rastrigin twin_gaussian random_gaussian himmelblau rosenbrock eggbox half_gaussian fitting gaussian_shell gaussian_shells

# Your likelihood programs
PROGRAMS = my_likelihood my_cpp_likelihood 



# The compiler type:
# This sets up the correct compiler flags below 
# current options are:
#  * gnu
#  * intel
COMPILER_TYPE=

# We can autodetect the compiler type on unix systems via the shell.
# if you want to override this then just run make with
# make COMPILER_TYPE=<your type>
# where <your time> is gnu or intel
ifeq "$(shell which gfortran >/dev/null; echo $$?)" "0"
COMPILER_TYPE=gnu
endif

ifeq "$(shell which ifort >/dev/null; echo $$?)" "0"
COMPILER_TYPE=intel
endif

# Whether to use MPI or not (leave blank if running without MPI)
MPI=1

# Whether to run in debugging mode
DEBUG=

# Whether to perform inter-procedural optimisation
# -- this results in a faster program, but much slower compilation
IPO = -ipo

# Host architecture (for intel)
HOST = -xHOST


# ============ MPI settings ================
ifdef MPI
# Using MPI
# ---------
# Set a define flag for MPI
FCFLAGS += -DMPI
CCFLAGS += -DMPI
# Choose the mpi wrapper
FC = mpif90
CC = mpicc
else
# No parallelisation
# ------------------
# Set the compiler to be that defined in COMPILER_TYPE
ifeq ($(COMPILER_TYPE),intel)
FC = ifort
CC = icc
endif

ifeq ($(COMPILER_TYPE),gnu)
FC = gfortran
CC = g++
endif

endif


# ============ COMPILER_TYPE settings ================
 
# Archive tool
AR = ar r

# ifort settings
# ==============
ifeq ($(COMPILER_TYPE),intel)
# default flags
# --------------
# fpp  : perform preprocessing
# fpic : shared object libraries
FCFLAGS += -fpp -fpic
CCFLAGS += -fpp -fpic

ifdef DEBUG
# Debugging mode
# --------------
# g              : enable gnu debugger compatibility
# O0             : no optimisation
# traceback      : create a backtrace on system failure
# check all      : all checks (whilst compiling)
# warn all       : all warnings (whilst running)
# ftrapuv        : Traps uninitialized variables by setting them to very large values
# debug all      : additional debugging information
# gen-interfaces : generate an interface block for each routine
# warn-interfaces: warn on these interface blocks
# fpe0           : halts program if dividing by zero, etc
FCFLAGS += -g -O0 -traceback -check all,noarg_temp_created -warn all -ftrapuv -debug all -gen-interfaces -warn-interfaces
# check-uninit   : check for uninitialised variables
CCFLAGS += -g -O0 -traceback -ftrapuv -debug all 
else
# Optimised mode
# --------------
#   ipo          : interprocedural optimization (optimize entire program)
#   O3           : maximum optimisation
#   no-prec-div  : slightly less precise floating point divides, but speeds up
#   static       : link intel libraries statically
#   xHost        : maximise architecture usage
#   w            : turn off all warnings
#   vec-report0  : disables printing of vectorizer diagnostic information
#   opt-report0  : disables printing of optimization reports
FCFLAGS += $(IPO) -O3 -no-prec-div $(HOST) -w -vec-report0 -opt-report0 
CCFLAGS += $(IPO) -O3 -no-prec-div $(HOST) -w -vec-report0 -opt-report0
ifdef IPO
# Archive tool for compiling with ipo.
AR = xiar r
LINKLIB = ld -shared

endif
endif

endif



# gfortran settings
# =================
ifeq ($(COMPILER_TYPE),gnu) 
# default flags
# --------------
# free-line-length-none : turn of line length limitation (why is this not a default??)
# cpp  					: perform preprocessing
# fPIC                  : for compiling a shared object library
FCFLAGS += -ffree-line-length-none -cpp -fPIC

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
#
CCFLAGS += -g -O0 -Wall -Wextra -Wshadow -Weffc++
else
# Optimised mode
# --------------
# Ofast : maximum optimisation
FCFLAGS += -Ofast
CCFLAGS += -Ofast
endif
LINKLIB = gfortran -shared -fPIC
endif

# Where polychord is stored
POLYCHORD_DIR = ./src

# Where likelihoods are stored
LIKELIHOOD_DIR = ./likelihoods

# Where likelihood examples are stored
EXAMPLE_LIKELIHOOD_DIR = $(LIKELIHOOD_DIR)/examples

# Where binaries are stored
BIN_DIR = ./bin

# Add these libraries to the program
LIBCHORD += -L$(POLYCHORD_DIR) -lchord
LIBLIKE += -L$(LIKELIHOOD_DIR) -L$(EXAMPLE_LIKELIHOOD_DIR)
EXT_LIBS += -lstdc++


INC += -I$(POLYCHORD_DIR) -I$(LIKELIHOOD_DIR) -I$(EXAMPLE_LIKELIHOOD_DIR)

# likelihood libraries, this is created by changing X to libX.a
EXAMPLE_LIKELIHOODS = $(patsubst %,lib%.a,$(EXAMPLES))
PROGRAM_LIKELIHOODS = $(patsubst %,lib%.a,$(PROGRAMS))


# Export all of the necessary variables
export DEBUG COMPILER_TYPE FCFLAGS MPI AR FC CC CCFLAGS EXAMPLE_LIKELIHOODS IPO EXT_LIBS LINKLIB


# "make" builds all
all: $(EXAMPLES) $(PROGRAMS)


# Rule for building polychord static library
libchord.a:
	$(MAKE) -C $(POLYCHORD_DIR) libchord.a

# Rule for building polychord shared library
libchord.so:
	$(MAKE) -C $(POLYCHORD_DIR) libchord.so

pypolychord.pyf:
	f2py -m pypolychord -h pypolychord/pypolychord.pyf src/interfaces.F90

pypolychord.so: pypolychord.pyf
	f2py -c -L./src/ -lchord -I./src/  pypolychord/pypolychord.pyf src/interfaces.F90
	mv pypolychord.so pypolychord


# Rule for building likelihood libraries
$(EXAMPLE_LIKELIHOODS) $(PROGRAM_LIKELIHOODS): libchord.a
	$(MAKE) -C $(LIKELIHOOD_DIR) $@

# Rule for example programs
$(EXAMPLES) $(PROGRAMS): %: libchord.a lib%.a  polychord.o
	$(FC) $(FCFLAGS) -o $(BIN_DIR)/$@ polychord.o $(LIBCHORD) $(LIBLIKE) -l$@ $(EXT_LIBS) 

# Rule for building main file
polychord.o: polychord.F90
	$(FC) $(FCFLAGS) $(INC) -c $< $(EXT_LIBS)  



clean:
	rm -f *.o *.mod *.MOD
	$(MAKE) -C $(LIKELIHOOD_DIR) clean 
	$(MAKE) -C $(POLYCHORD_DIR) clean
	
veryclean: clean
	rm -f *~ 
	$(MAKE) -C $(POLYCHORD_DIR) veryclean
	$(MAKE) -C $(LIKELIHOOD_DIR) veryclean
	cd $(BIN_DIR) && rm -f $(EXAMPLES) $(PROGRAMS)
