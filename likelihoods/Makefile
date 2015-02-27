# List of executables to be built within the package
PROGRAMS = libgaussian.a librastrigin.a

# What to include
INCLUDE += -I../src/

# "make" builds all
all: $(PROGRAMS)



lib%.a: %.o
	$(AR) $@ $^ 


# General rules for building prog.o from prog.f90 or prog.F90; $< is
# used in order to list only the first prerequisite (the source file)
# and not the additional prerequisites such as module or include files
%.o: %.f90
	$(FC) $(FCFLAGS) $(INCLUDE) -c $< 

%.o: %.F90
	$(FC) $(FCFLAGS) $(INCLUDE) -c $<



# Utility targets
.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD

veryclean: clean
	rm -f *~ $(PROGRAMS)
