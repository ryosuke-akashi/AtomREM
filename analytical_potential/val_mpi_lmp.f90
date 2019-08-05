module val_mpi
  implicit none
  integer :: ierr, myrank, nprocs
end module val_mpi
module val_mpi_lmp
  implicit none
  integer :: lammps, comm_lammps, ptr
end module val_mpi_lmp
