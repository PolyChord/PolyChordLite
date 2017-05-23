$(PYPOLYCHORD_DIR)/_PyPolyChord.so: _PyPolyChord.c
	python2 setup.py build_ext -b $(PYPOLYCHORD_DIR) -t .
	python3 setup.py build_ext -b $(PYPOLYCHORD_DIR) -t .

clean:
	$(RM) *.o

veryclean:
