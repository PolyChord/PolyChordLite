# rule for building the default likelihood
$(LIB_DIR)/libfortran_likelihood.a: likelihood.o
	$(AR) $@ $^ 


# General rule for building object file (.o)  from fortran files (.f90, .F90)
%.o: %.f90
	$(FC) $(FFLAGS) -I$(POLYCHORD_DIR) -c $< -o $@
%.o: %.F90
	$(FC) $(FFLAGS) -I$(POLYCHORD_DIR) -c $< -o $@

# Utility targets
.PHONY: clean veryclean

clean:
	$(RM) *.o *.mod *.MOD

veryclean: clean
	$(RM) *~
