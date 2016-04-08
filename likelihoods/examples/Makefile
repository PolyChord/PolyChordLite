# general rule for building each of the examples
# build the library lib<example name>.a  by archiving examples/<example name>.o using AR
$(patsubst %,$(LIB_DIR)/lib%.a,$(EXAMPLES)): $(LIB_DIR)/lib%.a: %.o
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
