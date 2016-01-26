# Object files:
# find all files ending in 90, and change .f90 and .F90 to .o, 
OBJECTS = $(patsubst %.F90,%.o,$(patsubst %.f90,%.o,$(wildcard *90)))

# include the dependencies
include makefile.dep

# The polychord library depends on the objects
# we should make all the objects and the archive them with AR
libchord.a: $(OBJECTS)
	$(AR) $@ $^
libchord.so: $(OBJECTS)
	$(LINKLIB) -o $@ $^

# General rule for building object file (.o)  from fortran files (.f90)
%.o: %.f90
	$(FC) $(FCFLAGS) -c $< $(EXT_LIBS)  

# General rule for building object file (.o)  from fortran files (.F90)
%.o: %.F90
	$(FC) $(FCFLAGS) -c $< $(EXT_LIBS)  

# Utility targets
.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD

veryclean: clean
	rm -f *~ *.a *.so
