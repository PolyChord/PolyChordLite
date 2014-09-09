DEBUG ?= 0
PAR   ?= 1
export DEBUG

default: all

libchord.a: ./src/*90
	cd ./src && make libchord.a

main: ./src/*90
	cd ./src && make libchord.a && make main

planck: ./src/*90
	cd ./src && make libchord.a && make planck

clean:
	cd ./src && make clean
	
veryclean:
	cd ./src && make veryclean

all: main planck
