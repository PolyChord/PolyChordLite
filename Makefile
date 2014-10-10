DEBUG ?= 0
export DEBUG

default: all

libchord.a: ./src/*90
	cd ./src && make libchord.a

main: ./src/*90
	cd ./src && make libchord.a && make main

batch: ./src/*90
	cd ./src && make libchord.a && make batch

planck: ./src/*90
	cd ./src && make libchord.a && make planck

clean:
	cd ./src && make clean
	
veryclean:
	cd ./src && make veryclean

all: batch main planck
