# CC=mpicc
# C_FLAGS= -Wall -pedantic -std=c99 -O2 -fopenmp
# LIBRARY= -lpthread -lm -fopenmp
CC=clang
C_FLAGS= -Wall -pedantic -std=c99 -O0 -g
LIBRARY= -lpthread -lm


.PHONY: all
all: ncg

ncg: IO.o Utilities.o Clifford.o Random.o Actions.o MonteCarlo.o main.o
	${CC} -o $@ $^ $(C_FLAGS) $(LIBRARY) #-pg

%.o: %.c
	${CC} -c $(C_FLAGS) $< -o $@

clean:
	rm -f *.o
	rm -f ncg
	rm -f benchmark.x
	rm -f evs
	rm -f dist
