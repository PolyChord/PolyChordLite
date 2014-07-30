default: all

main: ./src/*90
	cd ./src && make main && mv main ../

clean:
	cd ./src && make clean
	
veryclean:
	cd ./src && make veryclean

all: main
