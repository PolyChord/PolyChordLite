polychord_examples.o: polychord_examples.F90
	$(FC) $(FFLAGS) -I$(POLYCHORD_DIR) -I$(LIKELIHOOD_DIR) -c $< 

polychord.o: polychord.F90
	$(FC) $(FFLAGS) -I$(POLYCHORD_DIR) -I$(LIKELIHOOD_DIR) -c $< 

clean:
	$(RM) *.o *.mod *.MOD
