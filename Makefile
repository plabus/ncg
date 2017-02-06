CC=mpicc
C_FLAGS= -Wall -std=c99 -O2 -fopenmp -I${MKLROOT}/include #-pg
LIBRARY= -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a -Wl,--end-group -lpthread -lm -openmp


.PHONY: all
all: sim.x
bench: benchmark.x


sim.x: Clifford.o Random.o Actions.o MonteCarlo.o main.o
	${CC} -o $@ $^ $(C_FLAGS) $(LIBRARY) #-pg

benchmark.x: Clifford.o Random.o Actions.o MonteCarlo.o burn_in_benchmark.o
	${CC} -o $@ $^ $(C_FLAGS) $(LIBRARY) -pg



%.o: %.c
	${CC} -c $(C_FLAGS) $< -o $@


clean:
	rm -f *.o
	rm -f sim.x
	rm -f benchmark.x
	rm -f evs
	rm -f dist
