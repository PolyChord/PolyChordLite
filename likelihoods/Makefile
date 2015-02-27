# List of executables to be built within the package
PROGRAMS = libgaussian.a librastrigin.a

# "make" builds all
all: $(PROGRAMS)



lib%.a: %.o
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
	$(FC) $(FCFLAGS) -I../src/ -c $< 

%.o: %.F90
	$(FC) $(FCFLAGS) -I../src/ -c $<



# Utility targets
.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD

veryclean: clean
	rm -f *~ $(PROGRAMS)
