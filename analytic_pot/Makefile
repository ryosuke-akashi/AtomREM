F90   = mpiifort 
CC    = mpiicc 
OBJ   = bistable1d_val.o   lattice.o    val_mpi_lmp.o   mt19937.o   analytical_potential.o   main.o
SRC   = bistable1d_val.f90 lattice.f90  val_mpi_lmp.f90 mt19937.f90 analytical_potential.f90 main.f90
##OBJ   = libfwrapper2.o bistable1d_val.o lattice.o val_mpi_lmp.o mt19937.o potentials_lmp.o main.o
##SRC   = libfwrapper2.c bistable1d_val.f90 lattice.f90 val_mpi_lmp.f90 mt19937.f90 potentials_lmp.f90 main.f90
##LIB   = -L/home/iurii/src/lammps-5Jun19/src
##INC   = -I/home/iurii/src/lammps-5Jun19/src
#LIB   = -L/home/akashi/LAMMPS_lib/lammps-16Feb16/src
#INC   = -I/home/akashi/LAMMPS_lib/lammps-16Feb16/src
#LIBFLAG = -llammps -lstdc++
#SRC   = libfwrapper2.c val_mpi_lmp.f90 mt19937.f90 potentials_lmp.f90 main.f90
#DEBUG = -xHOST -O3 -ipo -static-intel -static_mpi #-fpe0
#DEBUG = -fast -pie #-fpe0
DEBUG = -g -traceback  -CB #-fpe0
#DEBUG = -O3#-g -traceback  #-CB
a.out:$(OBJ)
	$(F90) $(DEBUG) $(OBJ) -o $@
#a.out:$(OBJ)
#	$(F90) $(DEBUG) $(LIB) $(OBJ) $(LIBFLAG)  -o $@
#a.out:$(OBJ)
#	$(CPP) $(DEBUG) $(LIB) $(OBJ) $(LIBFLAG) $(LIBFFTW)  -o $@ $(LINK)
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
analytical_potential.o: analytical_potential.f90
	$(F90) $(DEBUG) -c $< 
##potentials_lmp.o: potentials_lmp.f90
##	$(F90) $(DEBUG) -c $< 
main.o: main.f90
	$(F90) $(DEBUG) -c $<
clean:
	rm in.cas* data.cas* Ener*  Chosen* Struct* After* Numb* fort* *.o *.mod log* out* dist* *.sh.* path.dat atomic* E_atoms* ENERGY/*  COORDINATES/* coordinates_ra* ; clear 
cleand:
	rm in.cas* data.cas* Ener*  Chosen* Struct* After* Numb* fort*   log* out* dist* *.sh.* path.dat atomic* E_atoms*  coordinates_ra* ; clear  

