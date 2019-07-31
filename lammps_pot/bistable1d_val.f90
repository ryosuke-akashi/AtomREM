module bistable1d_val
  implicit none

  integer, save :: &
  & Nstep,         & ! number of steps
  & Nwalker,       & ! number of walkers
  & Natoms,        & ! number of atoms 
  & Nregions,        & ! number of lammps calculations done at one time
  & Ncycles,        & ! number of lammps calculations
  & steptowrite      ! number of steps in order to write the path and other data


  real(8), save ::    &
  & temp,             & ! temperature [Kelvin]
  & tempFin,          & ! temperature [Kelvin] for final stage in "Reacton" and "Initialization" modes
  & dt,               & ! timestep
  & ratio,            & ! delta parameter in equation
  & ratio_ctrl          ! delta parameter in equation !! AKASHI

  character(256), save :: &
  & mode                ! initial mode  

  logical, save :: &
  & Vswitch , Dist_switch, Debag          ! on or off V, saving of distribution and debagging info

  integer, allocatable, save :: &
  & itype(:)          ! index of atom types

  real(8), allocatable, save :: &
  & x_pos(:,:), &
  & y_pos(:,:), &
  & z_pos(:,:)

  real(8), allocatable, save :: &
  & x_pos_mod(:,:), &
  & y_pos_mod(:,:), &
  & z_pos_mod(:,:)

  real(8), save :: &
  & a_orig(3)         ! origin position

  real(8), save :: &
  & a_vec(3,3),    &  ! lattice vector
  & b_vec(3,3)        ! reciprocal lattice vector (a_i*b_j=\delta_ij)

  real(8), save :: &
  & xlo, xhi,      &
  & ylo, yhi,      &
  & zlo, zhi,      & ! components defining the lattice (in lammps format)
  & xy, xz, yz

  character(1), save :: &
  & bounds(3)         ! finite(f), periodic(p)
end module bistable1d_val
