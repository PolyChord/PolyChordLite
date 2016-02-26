# rule for building the default likelihood
libmy_likelihood.a: my_likelihood.o
	$(AR) $(ARFLAGS) $@ $^ 

# rule for building the c++ wrapper
libmy_cpp_likelihood.a: my_cpp_likelihood.o fortran_cpp_wrapper.o
	$(AR) $(ARFLAGS) $@ $^ 



# General rule for building object file (.o)  from fortran files (.f90)
%.o: %.f90
	$(FC) $(FCFLAGS) -I$(POLYCHORD_DIR) -c $< -o $@

# General rule for building object file (.o)  from fortran files (.F90)
%.o: %.F90
	$(FC) $(FCFLAGS) -I$(POLYCHORD_DIR) -c $< -o $@

# General rule for building object file (.o)  from c++ files (.cpp)
%.o: %.cpp
	$(CC) $(CCFLAGS) -I$(POLYCHORD_DIR) -c $<

# Utility targets
.PHONY: clean veryclean

clean:
	$(RM) *.o *.mod *.MOD

veryclean: clean
	$(RM) *~ *.a
