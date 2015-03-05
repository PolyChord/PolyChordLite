# Include the polychord directory
INC = -I../src/

# Where examples are found
EXAMPLES_DIR = ./examples


# rule for building the default likelihood
libmy_likelihood.a: my_likelihood.o
	$(AR) $@ $^ 

# rule for building the c++ wrapper
libmy_cpp_likelihood.a: my_cpp_likelihood.o fortran_cpp_wrapper.o
	$(AR) $@ $^ 

# general rule for building each of the examples
# build the library lib<example name>.a  by archiving examples/<example name>.o using AR
$(EXAMPLE_LIKELIHOODS): lib%.a: $(EXAMPLES_DIR)/%.o
	$(AR) $@ $^ 



# General rule for building object file (.o)  from fortran files (.f90)
%.o: %.f90
	$(FC) $(FCFLAGS) $(INC) -c $< -o $@ $(EXT_LIBS)

# General rule for building object file (.o)  from fortran files (.F90)
%.o: %.F90
	$(FC) $(FCFLAGS) $(INC) -c $< -o $@ $(EXT_LIBS)

# General rule for building object file (.o)  from c++ files (.cpp)
%.o: %.cpp
	$(CC) $(CCFLAGS) -c $< $(EXT_LIBS)



# Utility targets
.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD
	cd $(EXAMPLES_DIR) && rm -f *.o *.mod *.MOD

veryclean: clean
	rm -f *~ *.a
	cd $(EXAMPLES_DIR) && rm -f *~ *.a 
