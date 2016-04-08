
$(DRIVERS_DIR)/polychord_examples.o: polychord_examples.F90
	$(FC) $(FFLAGS) -I$(POLYCHORD_DIR) -I$(EXAMPLES_DIR) -c $< 

$(patsubst %_likelihood,$(DRIVERS_DIR)/polychord_%.o,$(PROGRAMS)): $(DRIVERS_DIR)/polychord_%.o: polychord_%.F90
	$(FC) $(FFLAGS) -I$(POLYCHORD_DIR) -I$(LIKELIHOOD_DIR)/$* -c $<

clean:
	$(RM) *.o *.mod *.MOD

veryclean:
