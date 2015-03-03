# Object files
# ============
OBJECTS_PAR =                \
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
ifdef MPI
# MPI pararellised
# -----------------
OBJECTS = $(OBJECTS_PAR)
include makefile_mpi.dep
else
# No parallelisaton
# -----------------
#  strip out the mpi_utils file
OBJECTS = $(patsubst mpi_utils.o,,$(OBJECTS_PAR))
include makefile.dep
endif

# "make" builds all
all: libchord.a

# The polychord library
libchord.a: $(OBJECTS)
	$(AR) $@ $^

# General rule for building object file (.o)  from fortran files (.f90)
%.o: %.f90
	$(FC) $(FCFLAGS) -c $< 

# General rule for building object file (.o)  from fortran files (.F90)
%.o: %.F90
	$(FC) $(FCFLAGS) -c $<

# Utility targets
.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD

veryclean: clean
	rm -f *~ *.a
