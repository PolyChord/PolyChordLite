$(DRIVERS_DIR)/polychord_examples.o: polychord_examples.F90
	$(FC) $(FFLAGS) -I$(POLYCHORD_DIR) -I$(EXAMPLES_DIR) -c $< 

.PHONY: fortran_polychord CC_polychord

fortran_polychord: polychord.F90
	$(FC) $(FFLAGS) -I$(POLYCHORD_DIR) -I$(LIKELIHOOD_DIR)/fortran -c $<

CC_polychord: polychord.F90
	$(FC) $(FFLAGS) -I$(POLYCHORD_DIR) -I$(LIKELIHOOD_DIR)/CC -c $<

clean:
	$(RM) *.o *.mod *.MOD

veryclean:
