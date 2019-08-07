F90   = mpiifort 
CC    = mpiicc 
OBJ   = libfwrapper2.o bistable1d_val.o lattice.o val_mpi_lmp.o mt19937.o potentials_lmp.o main.o
LIB   = -L/home/akashi/LAMMPS_lib/lammps-5Jun19/src
INC   = -I/home/akashi/LAMMPS_lib/lammps-5Jun19/src
LIBFLAG = -llammps -lstdc++
SRC   = libfwrapper2.c bistable1d_val.f90 lattice.f90 val_mpi_lmp.f90 mt19937.f90 potentials_lmp.f90 main.f90
DEBUG = -g -traceback  -CB #-fpe0
a.out:$(OBJ)
	$(F90) $(DEBUG) $(LIB) $(OBJ) $(LIBFLAG)  -o $@
libfwrapper2.o: libfwrapper2.c
	$(CC)  $(DEBUG) $(INC)  -c $<
bistable1d_val.o: bistable1d_val.f90
	$(F90) $(DEBUG) -c $< 
lattice.o: lattice.f90
	$(F90) $(DEBUG) -c $< 
val_mpi_lmp.o: val_mpi_lmp.f90
	$(F90) $(DEBUG) -c $< 
mt19937.o: mt19937.f90
	$(F90) $(DEBUG) -c $< 
potentials_lmp.o: potentials_lmp.f90
	$(F90) $(DEBUG) -c $< 
main.o: main.f90
	$(F90) $(DEBUG) -c $<
clean:
	rm a.out *.o *.mod log* out* dist* *.sh.* path.dat atomic* E_atoms* ENERGY/*  COORDINATES/* ; clear 
