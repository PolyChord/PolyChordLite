# Object files:
# find all files ending in 90, and change .f90 and .F90 to .o, 
FOBJECTS = $(patsubst %.F90,%.o,$(patsubst %.f90,%.o,$(wildcard *90)))

OBJECTS = $(FOBJECTS)

# include the dependencies
include makefile.dep

# The polychord library depends on the objects
# we should make all the objects and the archive them with AR
$(LIB_DIR)/libchord.a: $(OBJECTS)
	$(AR) $@ $^

$(LIB_DIR)/libchord.so: $(OBJECTS)
	$(LD) -shared  $^ -o $@

# General rule for building object file (.o)  from fortran files (.f90/.F90)
%.o: %.f90
	$(FC) $(FFLAGS) -c $<
%.o: %.F90
	$(FC) $(FFLAGS) -c $<

# Utility targets
.PHONY: clean veryclean

clean:
	$(RM) *.o *.mod *.MOD

veryclean: clean
	$(RM) *~ *.a *.so
