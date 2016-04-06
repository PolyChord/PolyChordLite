# Object files:
# find all files ending in 90, and change .f90 and .F90 to .o, 
CXXOBJECTS = $(patsubst %.cpp,%.o,$(wildcard *.cpp))
COBJECTS = $(patsubst %.c,%.o,$(wildcard *.c))

OBJECTS = $(COBJECTS) $(CXXOBJECTS)

# include the dependencies
include makefile.dep

# The polychord library depends on the objects
# we should make all the objects and the archive them with AR
c_interface.a: $(OBJECTS)
	$(AR) $(LIB_DIR)/$@ $^

c_interface.so: $(OBJECTS)
	$(LD) -shared  $^ -o $(LIB_DIR)/$@

# General rule for building object file (.o)  from C files (.c)
%.o: %.c
	$(CC) $(CFLAGS) -c $<
# General rule for building object file (.o)  from CC files (.cpp)
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<

# Utility targets
.PHONY: clean veryclean

clean:
	$(RM) *.o *.mod *.MOD

veryclean: clean
	$(RM) *~ *.a *.so
