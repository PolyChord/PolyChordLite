
$(DRIVERS_DIR)/polychord_examples.o: polychord_examples.f90
	$(FC) $(FFLAGS) -I$(POLYCHORD_DIR) -I$(EXAMPLES_DIR) -c $< 

$(DRIVERS_DIR)/polychord_CC.o: polychord_CC.cpp
	$(CXX) $(CXXFLAGS) -I$(POLYCHORD_DIR) -I$(LIKELIHOOD_DIR)/CC -c $<

$(DRIVERS_DIR)/polychord_fortran.o: polychord_fortran.f90
	$(FC) $(FFLAGS) -I$(POLYCHORD_DIR) -I$(LIKELIHOOD_DIR)/fortran -c $<

clean:
	$(RM) *.o *.mod *.MOD

veryclean:
