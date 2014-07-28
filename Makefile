default: all

galileo: ./src/*90
	cd ./src && make galileo && mv galileo ../

clean:
	cd ./src && make clean
	
veryclean:
	cd ./src && make veryclean

all: galileo
