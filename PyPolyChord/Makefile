SRC=src

$(PYPOLYCHORD_DIR)/_PyPolyChord.so: $(PYPOLYCHORD_DIR)/src/_PyPolyChord.c
	$(MAKE) -C $(SRC) $@

clean:
	$(MAKE) -C $(SRC) clean

veryclean:
	$(RM) *.so *.pyc __pycache__
	$(MAKE) -C $(SRC) veryclean
