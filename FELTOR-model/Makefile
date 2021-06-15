device=omp

#configure machine
include config/default.mk # default configuration
include config/*.mk #

include config/devices/devices.mk # recombines variables for GPU and MPI

all: convection convection_hpc

convection: convection.cpp convection.h
	$(CC) $(OPT) $(CFLAGS) $< -o $@ $(INCLUDE) $(GLFLAGS) $(JSONLIB) -DDG_BENCHMARK  -g

convection_hpc: convection_hpc.cpp convection.h
	$(CC) $(OPT) $(CFLAGS) $< -o $@ $(INCLUDE) $(LIBS) $(JSONLIB) -DDG_BENCHMARK  -g

convection_mpi: convection_hpc.cpp convection.h
	$(MPICC) $(OPT) $(MPICFLAGS) $< -o $@ $(INCLUDE) $(LIBS) $(JSONLIB) -DFELTOR_MPI -DDG_BENCHMARK -g

.PHONY: clean

clean:
	rm -rf convection convection_hpc convection_mpi
