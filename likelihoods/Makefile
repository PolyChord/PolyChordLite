# What to include
POLYCHORD += -I../src/

libcpp_likelihood.a: cpp_likelihood.o likelihood_cpp.o

# general rule for building the archive
lib%.a: %.o
	$(AR) $@ $^ 

# General rules for building prog.o from prog.f90 or prog.F90; $< is
# used in order to list only the first prerequisite (the source file)
# and not the additional prerequisites such as module or include files
%.o: %.f90
	$(FC) $(FCFLAGS) $(POLYCHORD) -c $< 

%.o: %.F90
	$(FC) $(FCFLAGS) $(POLYCHORD) -c $<

# For making cpp files
%.o: %.cpp
	$(CC) $(CCFLAGS) -c $<



# Utility targets
.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD

veryclean: clean
	rm -f *~ *.a
