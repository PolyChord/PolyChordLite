# rule for building the c++ wrapper
$(LIB_DIR)/libCC_likelihood.a: CC_likelihood.o fortran_cpp_wrapper.o
	$(AR) $@ $^ 



# General rule for building object file (.o)  from fortran files (.f90, .F90)
%.o: %.f90
	$(FC) $(FFLAGS) -I$(POLYCHORD_DIR) -c $< -o $@
%.o: %.F90
	$(FC) $(FFLAGS) -I$(POLYCHORD_DIR) -c $< -o $@
# General rule for building object file (.o)  from c++ files (.cpp)
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -I$(POLYCHORD_DIR) -I$(C_INTERFACE_DIR) -c $<
# General rule for building object file (.o)  from c files (.c)
%.o: %.c
	$(CC) $(CFLAGS) -I$(POLYCHORD_DIR) -I$(C_INTERFACE_DIR) -c $<

# Utility targets
.PHONY: clean veryclean

clean:
	$(RM) *.o *.mod *.MOD

veryclean: clean
	$(RM) *~ *.a
