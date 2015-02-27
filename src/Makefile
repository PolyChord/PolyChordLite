# Object files
# ============
ifdef MPI
# Using MPI
# ---------
# Define the files to be used 
OBJECTS =                    \
		  utils.o            \
		  abort.o            \
		  mpi_utils.o        \
		  array_utils.o      \
		  random_utils.o     \
		  params.o           \
		  settings.o         \
		  priors.o           \
		  calculate.o        \
		  run_time_info.o    \
		  feedback.o         \
		  read_write.o       \
		  clustering.o       \
		  chordal_sampling.o \
		  generate.o         \
		  nested_sampling.o  \
		  ini.o
else

# No parallelisation
# ------------------
# Define the files to be used (note no mpi_utils file)
OBJECTS =                    \
		  utils.o            \
		  abort.o            \
		  array_utils.o      \
		  random_utils.o     \
		  params.o           \
		  settings.o         \
		  priors.o           \
		  calculate.o        \
		  run_time_info.o    \
		  feedback.o         \
		  read_write.o       \
		  clustering.o       \
		  chordal_sampling.o \
		  generate.o         \
		  nested_sampling.o  \
		  ini.o

endif


# List of executables to be built within the package
PROGRAMS = libchord.a

# "make" builds all
all: $(PROGRAMS)

# The polychord library
libchord.a: $(OBJECTS)
	$(AR) $@ $^

# General rule for building prog from prog.o; $^ (GNU extension) is
# used in order to list additional object files on which the
# executable depends
%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ 

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
