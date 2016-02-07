# Include the polychord directory
INC = -I../src/

# general rule for building each of the examples
# build the library lib<example name>.a  by archiving examples/<example name>.o using AR
$(EXAMPLE_LIBRARIES): lib%.a: %.o
	$(AR) $@ $^ 



# General rule for building object file (.o)  from fortran files (.f90, .F90)
%.o: %.f90
	$(FC) $(FFLAGS) $(INC) -c $< -o $@
%.o: %.F90
	$(FC) $(FFLAGS) $(INC) -c $< -o $@
# General rule for building object file (.o)  from c++ files (.cpp)
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<



# Utility targets
.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD

veryclean: clean
	rm -f *~ *.a
