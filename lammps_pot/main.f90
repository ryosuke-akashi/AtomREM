!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                                                       !
!        COPY RIGHT:    GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007                                              !
!                                                                                                                       !
!        AUTHORS:       Yuri NAGORNOV  &  RYOSUKE AKSHI                                                                 !
!                       University of Tokyo, Department of Physics                                                      !
!                                                                                                                       !
!        This research was supported by MEXT as Exploratory Challenge on Post-K                                         !
!        computer (Frontiers of Basic Science: Challenging the Limits).                                                 !
!        Project ID:    hp160257, hp170244, hp180184                                                                    !
!        "POST-K PROJECT"                                                                                               !
!                                                                                                                       !
!                                                                                                                       !
!        There is a _THEORY_ of method:                                                                                 !
!                  J. Phys. Soc. Jpn. 87, 063801 (2018) - https://journals.jps.jp/doi/abs/10.7566/JPSJ.87.063801        !
!        There is an _ALGORITHM_ of method and examples of simulations:                                                 !
!                  Yu.S.Nagornov, R.Akashi Physica A, 528, 121481 (2019) - https://doi.org/10.1016/j.physa.2019.121481  !
!                                                                                                                       !
!        This code has three modes of work:                                                                             !
!                       1 - Initialization - seekeng the initial positions (initialization mode),                       !
!                       2 - Reaction       - simulation under biasing potential (reaction path mode),                   !
!                       3 - Langevin       - simulation of Langevin mechanics (relaxation mode)                         !
!                                                                                                                       !
!                                                                                                                       !
!        The log scale for calculation of P and Q functions is used to avoid the overflow problem                       !
!                                                                                                                       !      
!        The reaction path is drawn by saving of maxima of Q distribution, in other words                               ! 
!        the reaction path is the average atomic coordinates and atomic energies.                                       !
!                                                                                                                       !
!        The code utilized the LAMMPS  (https://lammps.sandia.gov) as a shared library for potential calculation using  !
!        the diagonalHESSIAN package (reduction from "compute_partialHESSIAN")                                          !
!                                                                                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module routines
  implicit none
  include 'mpif.h'
  contains

   subroutine initmpi()
   use val_mpi, only : ierr, myrank, nprocs

   call MPI_INIT(ierr)

   call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
   call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)


   end subroutine initmpi

   subroutine init_lmp(restart)
   use bistable1d_val, only : Nwalker, Natoms, Nregions, Ncycles,itype, x_pos, y_pos, z_pos, &
 &                            a_orig, a_vec, bounds
   use val_mpi, only : ierr, myrank, nprocs
   use val_mpi_lmp, only : lammps, comm_lammps, ptr

   logical, intent(in) :: restart

   character (len=1024) :: line
   integer, parameter :: fp = 20
   integer, parameter :: fpair = 21
   integer :: n
   character(len=256), parameter :: fin = "in.case"

   real(8), parameter :: eunit=  1.671d0 ! 10^-14 erg instead 1.0d0
   !real(8), parameter :: dunit= 1.0d0
   real(8), parameter :: dunit= 3.4d0
   !real(8), parameter :: dcut = 2.0d0
   real(8), parameter :: dcut = 8.0d0
   real(8), parameter :: dcut2 = dcut*2d0
   !integer, parameter :: MAX_REGION = 16384 !! 
   integer, parameter :: MAX_REGION = 1 !! 

   real(8), parameter :: xmin = -dcut
   real(8), parameter :: ymin = -dcut
   real(8), parameter :: zmin = -dcut

   integer :: iw, idx, ina, ity
   character(16), allocatable :: id_pe(:)
   character(16), allocatable :: id_re(:)
   character(16) :: ich
  

   !! determine Ngroup (Number of lammps calculation done at the same time)
   IF(.not.restart)THEN
    IF(Nwalker.le.MAX_REGION) THEN
     Nregions = Nwalker
    ELSE
     Nregions = MAX_REGION
    ENDIF
    Ncycles = ceiling(dble(Nwalker)/dble(MAX_REGION))
   ELSE ! restart == .true.
    Nregions = 1
    Ncycles = Nregions
   ENDIF

   IF(myrank == 0) THEN
    allocate(id_pe(Nregions), id_re(Nregions))
   ENDIF
   !!/determine Ngroup (Number of lammps calculation done at the same time)

   !! Generate data.case file
   IF(myrank==0)THEN
    open(unit=fp,file="data.case", status="unknown")

    write(fp,*)"Input for lammps"
    write(fp,*)""
    write(fp,'(i10, a6)')(Natoms*Nregions), "atoms"
    write(fp,'(i10, a11)') maxval(itype(:)), "atom types"
    write(fp,'(2f16.6, a8)') a_orig(1),a_orig(1)+a_vec(1,1), "xlo xhi"   
    write(fp,'(2f16.6, a8)') a_orig(2),a_orig(2)+a_vec(2,2), "ylo yhi"
    write(fp,'(2f16.6, a8)') a_orig(3),a_orig(3)+a_vec(3,3), "zlo zhi"
    write(fp,'(3f16.6, a10)') a_vec(2,1),a_vec(3,1),a_vec(3,2), "xy xz yz"
    write(fp,*)""

    write(fp,*)"Masses"
    write(fp,*)""
    DO ity = 1, maxval(itype(:))
     write(fp,'(i10, f20.10)') ity, 1.0
    ENDDO
    write(fp,*)""
    write(fp,*)"Atoms"
    write(fp,*)""
    idx = 0
    DO iw = 1, Nregions
     DO ina = 1, Natoms
      idx = idx + 1
      write(fp,'(2i10, 3f20.10)') idx,  itype(ina),  x_pos(iw, ina), y_pos(iw, ina),&
&                               z_pos(iw, ina) 
     ENDDO
    ENDDO
    close(fp)

   ENDIF
   !!/Generate data.case file

   !! Generate lammps input file "in.case"

   IF(myrank==0)THEN
    DO iw = 1, Nregions
     write(ich,'(i0)')iw
     id_pe(iw) = "pe" //  trim(ich)
     id_re(iw) = "re" //  trim(ich)
    ENDDO
    open(unit=fp,file=fin, status="unknown")

    write(fp,*)"units          lj"
    write(fp,'(a9,3x,a1,1x,a1,1x,a1)')"boundary", bounds(1),bounds(2),bounds(3)
!    write(fp,*)"boundary    f f f"
    write(fp,*)"atom_style  atomic"
    write(fp,*)"atom_modify  map array"
    write(fp,*)"atom_modify  sort 0 0.0"
    write(fp,*)""
    write(fp,*)"read_data data.case"
    write(fp,*)""
    open(unit=fpair, file="in.pair", status="old", iostat=ierr)
    DO 
     read(unit=fpair, fmt='(A)', iostat=ierr) line
     IF(ierr.ne.0)exit
     write(fp, *)trim(line)
    ENDDO
    close(fpair)
!    write(fp,'(a10, 2a6, 3f9.4)')"pair_coeff", "*", "*", eunit, dunit, dcut
!    write(fp,'(a10, 2a6, a15, a6)')"pair_coeff", "*", "*", "CH.airebo", "C H"
!    DO iw = 1, Nregions
!     write(fp,'(a7, a7, a7, 6f16.8)')"region",trim(id_re(iw)),&
!&                               "block",xmin+dcut2*(iw-1), xmin+dcut2*iw, &
!&                               ymin, ymin+dcut2, zmin, zmin+dcut2     
!    ENDDO
    
    write(fp,*)""
    write(fp,*)"neighbor      0.3  bin"
    write(fp,*)"compute       peatom all pe/atom"
    write(fp,*)""

    write(fp,*)"compute         hsmtrxD    all  diagonalHessian  5e-9"

    write(fp,*)""
    write(fp,*)"compute         frcX    all property/atom fx"
    write(fp,*)"compute         frcY    all property/atom fy"
    write(fp,*)"compute         frcZ    all property/atom fz"

!    write(fp,*)"thermo_style custom step temp etotal c_hsmtrxD[*] "
!    write(fp,*)"thermo_modify  line multi format float %20.15g  flush yes"


!    write(fp,*)"  thermo         1    "

! PLEASE add log non after debagging
    write(fp,*)"log none"
   
    close(fp)
    deallocate(id_pe, id_re)
   ENDIF
   call MPI_BARRIER(MPI_COMM_WORLD, ierr)
   !!/Generate lammps in file

   
   !! MPI split for lammps
   lammps = myrank
   call MPI_COMM_SPLIT(MPI_COMM_WORLD, lammps, 0, comm_lammps, ierr)    
  
   call lammps_open(comm_lammps,ptr)
   !!/MPI split for lammps

   !! Open LAMMPS initial files
   open(unit=fp, file=fin, action='read', status='old', iostat=ierr)
   n = 0
   DO
    IF (myrank == 0)THEN
     read(unit=fp, fmt='(A)', iostat=ierr) line
     n = 0
     IF (ierr == 0) THEN
      n = len(trim(line))
      IF (n==0) THEN
      line = ' '
       n = 1
      ENDIF
     ENDIF
    ENDIF



    call MPI_BCAST(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)



    IF(n==0) EXIT
    call MPI_BCAST(line, n, MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
!   IF (lammps == 1) call lammps_command(ptr, line, n)



    call lammps_command(ptr, line, n)
   ENDDO

   close(fp)
   !!/Open LAMMPS initial files

   end subroutine init_lmp

   subroutine restart_lmp()
   use val_mpi, only: ierr, myrank, nprocs

   call finalize_lmp()

   IF(myrank.eq.0)THEN
    call system('mv in.case in.case1')
    call system('mv data.case data.case1')
   ENDIF
   call MPI_barrier(MPI_COMM_WORLD,ierr)
   call init_lmp(.true.)

   end subroutine restart_lmp


   subroutine finalize_lmp()
   use val_mpi_lmp, only :  ptr
  
   call lammps_close(ptr)

   end subroutine finalize_lmp
  !-----------
   subroutine stdin()
   use bistable1d_val, only : Nstep, Nwalker, Natoms, temp, tempFin, steptowrite, dt, ratio, &
  &                           mode, fin, fout,  Dist_switch, Debag, Vswitch,        &
  &                           a_orig, a_vec, b_vec, bounds,                         &
  &                           xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz

   use val_mpi, only: ierr, myrank, nprocs

   real(8) :: x_orig, y_orig, z_orig
   real(8) :: a_vec1(3), a_vec2(3), a_vec3(3)   
   namelist /input/ Nstep, Nwalker, Natoms, temp, tempFin, steptowrite,                      &
  &                 dt, ratio, mode, fin, fout,                                     &
  &                 x_orig, y_orig, z_orig, a_vec1, a_vec2, a_vec3, bounds
  
   mode = trim(mode)

   IF(myrank.eq.0)THEN
    open(unit=10,file="param_3N.in",status="unknown")
    read(unit=10,nml=input,err=100)
    close(10)
    write(*,*)"Nwalker, nprocs=",Nwalker, nprocs
    IF(mod(Nwalker,nprocs).ne.0)THEN
     write(*,*)"Error: Number of walkers must be multiple of nprocs"
     stop
    ENDIF
    Nwalker = Nwalker/nprocs
    
    ! Define the logical variables for each mode

    Dist_switch = .FALSE.       ! logical variable to construct the distribution, 
                                ! it takes a lot of hard disk space and memory, 
                                ! and utilized only for debagging  
    Debag = .FALSE.             ! the same reason - it's only for debagging process


    mode = trim(mode)

    IF(mode=="Initialization")THEN
       Vswitch = .FALSE. 

      ELSE IF(mode=="Reaction")THEN
       Vswitch = .TRUE.

        ELSE IF(mode=="Langevin")THEN
          Vswitch = .FALSE.

             ELSE
               write(*,*)"error: mode = invalid"
               write(*,*)"mode could be Initialization or Reaction or Langevin, but mode = ",mode

    ENDIF

    ! origin of the space
    a_orig(1) = x_orig
    a_orig(2) = y_orig
    a_orig(3) = z_orig
    write(*,*)"a_orig",a_orig(:)
    ! calculate reciprocal lattice vectors so that a_i*b_j=\delta_ij
    a_vec(1,:) = a_vec1(:)
    a_vec(2,:) = a_vec2(:)
    a_vec(3,:) = a_vec3(:)
    write(*,*)"a_vec1",a_vec(1,:)
    write(*,*)"a_vec2",a_vec(2,:)
    write(*,*)"a_vec3",a_vec(3,:)
    call recvec(a_vec, b_vec)
    write(*,*)"b_vec1",b_vec(1,:)
    write(*,*)"b_vec2",b_vec(2,:)
    write(*,*)"b_vec3",b_vec(3,:)
    
    !/calculate reciprocal lattice vectors so that a_i*b_j=\delta_ij
    ! define xlo, etc.
    ! NOTE!! only a1=(xx, 0, 0), a2=(xy, yy, 0), a3=(xz, yz, zz) is accepted
    ! according to the LAMMPS format
    xlo = a_orig(1)
    xhi = a_orig(1) + a_vec(1,1)
    ylo = a_orig(2)
    yhi = a_orig(2) + a_vec(2,2)
    zlo = a_orig(3)
    zhi = a_orig(3) + a_vec(3,3)
     xy = a_vec(2,1)
     xz = a_vec(3,1)
     yz = a_vec(3,2)
    !/define xlo, etc.
    ! check boundary condition
    write(*,*)"bounds"," ", bounds(1)," ", bounds(2)," ", bounds(3)
    !/check boundary condition
    write(*,*)"mode of simulation = ",mode
   ENDIF

   call MPI_BCAST(Nstep   , 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ierr)
   call MPI_BCAST(Nwalker , 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ierr)
   call MPI_BCAST(Natoms  , 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ierr)
   call MPI_BCAST(temp    , 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
   call MPI_BCAST(tempFin , 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
   call MPI_BCAST(steptowrite   , 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ierr)
   call MPI_BCAST(dt      , 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
   call MPI_BCAST(ratio   , 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
   call MPI_BCAST(mode    , 256, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
   call MPI_BCAST(Vswitch , 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
   call MPI_BCAST(Dist_switch , 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
   call MPI_BCAST(Debag , 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
   call MPI_BCAST(xlo , 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
   call MPI_BCAST(xhi , 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
   call MPI_BCAST(ylo , 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
   call MPI_BCAST(yhi , 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
   call MPI_BCAST(zlo , 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
   call MPI_BCAST(zhi , 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
   call MPI_BCAST(xy , 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
   call MPI_BCAST(xz , 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
   call MPI_BCAST(yz , 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
   call MPI_BCAST(bounds , 3, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
   call MPI_BCAST(a_vec , 9, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
   call MPI_BCAST(b_vec , 9, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
   return
100 write(*,*) "ERROR: reading namelist file"
   stop
   end subroutine stdin
  !-----------

  !-----------
   subroutine gen_walker()

     use mt19937
     use bistable1d_val, only: Nwalker, Natoms,itype, x_pos, y_pos, z_pos, mode, fin, &
  &                            x_pos_mod, y_pos_mod, z_pos_mod
     use val_mpi, only: ierr, myrank, nprocs

     integer iw,ina
     real(8) :: ranx
     real(8) :: rany
     real(8) :: ranz

     ALLOCATE(itype(Natoms))
     ALLOCATE(x_pos(Nwalker,Natoms))
     ALLOCATE(y_pos(Nwalker,Natoms))
     ALLOCATE(z_pos(Nwalker,Natoms))

     !! array for subroutine lattice
     ALLOCATE(x_pos_mod(Nwalker,Natoms))
     ALLOCATE(y_pos_mod(Nwalker,Natoms))
     ALLOCATE(z_pos_mod(Nwalker,Natoms))
     !!/array for subroutine lattice

     IF(myrank.eq.0)THEN

       open(unit=10,file=fin, status="unknown")
       rewind(10)
       DO ina=1, Natoms
        read(10,*)itype(ina),x_pos(1,ina), y_pos(1,ina), z_pos(1,ina)
        x_pos(2:Nwalker,ina) = x_pos(1,ina)
        y_pos(2:Nwalker,ina) = y_pos(1,ina)
        z_pos(2:Nwalker,ina) = z_pos(1,ina)
       ENDDO
       close(10)

     ENDIF
     call MPI_BCAST(x_pos   , Nwalker*Natoms, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
     call MPI_BCAST(y_pos   , Nwalker*Natoms, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
     call MPI_BCAST(z_pos   , Nwalker*Natoms, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
     call MPI_BCAST(itype   , Natoms, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

   end subroutine gen_walker
  !-----------

   subroutine gen_noise(nw,na,dx,dy,dz)
     use mt19937
     integer, intent(in) :: nw,na
     real(8), intent(out) :: dx(nw,na)
     real(8), intent(out) :: dy(nw,na)
     real(8), intent(out) :: dz(nw,na)

     real(8), parameter :: pi=3.14159265359d0
     integer :: iw,ina
     real(8) :: xval1, xval2, yval1, yval2, zval1, zval2

     !! Box-Muller
    DO ina=1, na
     DO iw=1, nw

        xval1=genrand_real3()
        xval2=genrand_real3()
        dx(iw,ina)=sqrt(-2d0*log(xval1))*cos(2d0*pi*xval2)

        yval1=genrand_real3()
        yval2=genrand_real3()
        dy(iw,ina)=sqrt(-2d0*log(yval1))*cos(2d0*pi*yval2)

        zval1=genrand_real3()
        zval2=genrand_real3()
        dz(iw,ina)=sqrt(-2d0*log(zval1))*cos(2d0*pi*zval2)

     ENDDO
    ENDDO

   end subroutine gen_noise
  !-----------
   subroutine solve()

   use bistable1d_val, only: Nstep, Nwalker,  Natoms,Nregions,Ncycles, itype,x_pos, y_pos, z_pos, steptowrite, dt, mode, ratio,ratio_ctrl, fout, temp, tempFin, Vswitch, Dist_switch, Debag, xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz, bounds
   use val_mpi, only: ierr, myrank, nprocs
   use potentials, only: calc_Udiff
   use mt19937

     real(8) :: dWX(Nwalker,Natoms)
     real(8) :: dWY(Nwalker,Natoms)
     real(8) :: dWZ(Nwalker,Natoms)

     real(8) :: mom1(Nwalker)
     real(8) :: mom1Natoms(Nwalker,Natoms)

     real(8) :: UdiffX(Nwalker,Natoms)
     real(8) :: UdiffY(Nwalker,Natoms)
     real(8) :: UdiffZ(Nwalker,Natoms)

     real(8) :: Udiff2X(Nwalker,Natoms)
     real(8) :: Udiff2Y(Nwalker,Natoms)
     real(8) :: Udiff2Z(Nwalker,Natoms)

     real(8) :: Vdiff0(Nwalker)

     real(8) :: Vexp(Nwalker)

     real(8) :: VdiffX(Nwalker,Natoms)
     real(8) :: VdiffY(Nwalker,Natoms)
     real(8) :: VdiffZ(Nwalker,Natoms)

     real(8) :: Vdiff2X(Nwalker,Natoms)
     real(8) :: Vdiff2Y(Nwalker,Natoms)
     real(8) :: Vdiff2Z(Nwalker,Natoms)

     real(8) :: weight(Nwalker)
     real(8) :: F_func(Nwalker)

     real(8) :: A_coeffX(Nwalker,Natoms)
     real(8) :: A_coeffY(Nwalker,Natoms)
     real(8) :: A_coeffZ(Nwalker,Natoms)

     real(8) :: prob(Nwalker)
     integer :: id_walker(Nwalker)

     integer :: max_d(3,Natoms)
     integer :: max_d_corr(3,Natoms)

     real(8) :: x_pos_tmp(Nwalker*nprocs, Natoms)
     real(8) :: y_pos_tmp(Nwalker*nprocs, Natoms)
     real(8) :: z_pos_tmp(Nwalker*nprocs, Natoms)

     real(8) :: x_ave_q(Natoms)
     real(8) :: y_ave_q(Natoms)
     real(8) :: z_ave_q(Natoms)
     real(8) :: x_ave_p(Natoms)
     real(8) :: y_ave_p(Natoms)
     real(8) :: z_ave_p(Natoms)
     logical :: Voff

     real(8)  :: coorX_max(Natoms)
     real(8)  :: coorY_max(Natoms)
     real(8)  :: coorZ_max(Natoms)

     integer, parameter :: Nx=30
     integer, parameter :: Ny=30
     integer, parameter :: Nz=30

     real(8) :: stepXdist(Nx,Natoms)
     real(8) :: stepYdist(Ny,Natoms)
     real(8) :: stepZdist(Nz,Natoms)

     real(8) :: dist_Q(1:Nx,1:Ny,1:Nz,1:Natoms)
     real(8) :: dist_P(1:Nx,1:Ny,1:Nz,1:Natoms)

     real(8) :: Vexp_ave

     real(8) :: displacement(1:Natoms),x_ini_r(1:Natoms), y_ini_r(1:Natoms), z_ini_r(1:Natoms)

     real(8) :: lnweight_min
     real(8) :: B_coeff
    
     real(8) :: Energy_q
     real(8) :: Energy_p
     real(8) :: E_atoms(Natoms)
     real(8) :: sq_force
     real(8) :: min_abs_force
     real(8) :: forceX(Natoms)
     real(8) :: forceY(Natoms)
     real(8) :: forceZ(Natoms)
    
     real(8) :: V_min
     real(8) :: V_ini
     integer :: Nw_tmp
     integer :: incr

     real(8) :: F_ave
     real(8) :: pr
     
     real(8) :: pr_accum(Nwalker*nprocs)
     
     real(8) :: rand_p
     real(8) :: inv_G
     real(8) :: lc

     integer :: istep
     integer :: iw,jw
     integer :: ix
     integer :: iy
     integer :: iz
     integer :: c,ina, i_ratio 
     integer :: i_restart,N_restart,restart_step
     real(8) :: ratio_ini,E_start

     !! mpi process local values
     real(8) :: V_min_pr
     real(8) :: Vexp_ave_pr
     real(8) :: F_ave_pr
     real(8) :: min_abs_force_pr
     integer :: n1, n2, npr
     integer :: iproc
     integer :: myincr
     real(8) :: pr_accum_pr(Nwalker*nprocs)
     real(8) :: x_ave_pr_q(Natoms)
     real(8) :: y_ave_pr_q(Natoms)
     real(8) :: z_ave_pr_q(Natoms)
     real(8) :: x_ave_pr_p(Natoms)
     real(8) :: y_ave_pr_p(Natoms)
     real(8) :: z_ave_pr_p(Natoms)
     real(8) :: w_sum_pr
     real(8) :: w_sum
     real(8) :: x_min_pr
     real(8) :: y_min_pr
     real(8) :: z_min_pr
     real(8) :: x_max_pr
     real(8) :: y_max_pr
     real(8) :: z_max_pr
     logical :: Com_switch
     !!/mpi process local values
     character(6)  :: istring
     character(64) :: f_rank

     !! mpi shared values with pudding
     real(8) :: x_pos_tmp2(Nwalker*nprocs, Natoms)
     real(8) :: y_pos_tmp2(Nwalker*nprocs, Natoms)
     real(8) :: z_pos_tmp2(Nwalker*nprocs, Natoms)
     !!/mpi shared values with pudding
     integer :: ct0, ct1, count_rate, count_max
     integer :: ct2, ct3, ct4, ct5, ct6

     !!       The part to find initial points with a fix two atoms which have a highest
     !        energy deviations from equilibrium state 
     !        to find initial points
     integer :: istruct,N_struct,UniqueNumber
     integer :: time_second_stage 
     real(8) :: timer_s, timer_f

     real(8) :: Energy_min(Natoms),Energy(Natoms)
     real(8) :: E_struct(Nwalker*nprocs,Natoms),U_write(Nwalker*nprocs),E_pr_struct(Nwalker*nprocs,Natoms),U_write_pr(Nwalker*nprocs)
     real(8) :: x_pos_str(Nwalker*nprocs,Natoms),y_pos_str(Nwalker*nprocs,Natoms),z_pos_str(Nwalker*nprocs,Natoms)
     integer :: number_struct(Nwalker*nprocs),number_pr_struct(Nwalker*nprocs),time_pr(Nwalker*nprocs)

     logical :: check_crossing(Nwalker),check_str(Nwalker*nprocs),check_walker(Nwalker),check_pr_str(Nwalker*nprocs)

     real(8) :: dE  ! , parameter :: dE = 0.01d0
     real(8) :: E_threshold ! , parameter :: E_threshold = -3.92d0
     real(8) :: E_cut  ! , parameter :: E_cut = -3.8d0

     character(32)  :: f_walker
     
     !!/to find the initial points !!!!!!!

      ! Initialization for seeking of the start points

      n1 = myrank*Nwalker + 1
      n2 = n1 + Nwalker - 1



      IF(mode=="Initialization")THEN
       U_write_pr(:) = 0.0d0

       open(unit=48, file="Energies_of_structures.dat",status="unknown")

       !n1 = myrank*Nwalker + 1
       !n2 = n1 + Nwalker - 1
       U_write(:) = 0.0d0
       !time(:) = 0
       time_pr(:) = 0
       number_struct(:) = 0
       x_pos_str(:,:) = 0.0d0
       y_pos_str(:,:) = 0.0d0
       z_pos_str(:,:) = 0.0d0
       E_struct(:,:) = 0.0d0
       E_pr_struct(:,:) = 0.0d0
       Energy_min(:) = 0.0d0
       Energy(:) = 0.0d0
       time_second_stage = 0

       check_crossing(:) = .FALSE.
       check_walker(:) = .FALSE.
       check_str(:) = .FALSE.
       check_pr_str(:) = .FALSE.

       ! calc  calc_Udiff0_atom
       call calc_Udiff(Nwalker,Natoms,Nregions,Ncycles,x_pos(1:Nwalker,1:Natoms),&
       &y_pos(1:Nwalker,1:Natoms),z_pos(1:Nwalker,1:Natoms),UdiffX(1:Nwalker,1:Natoms))

       IF (myrank.eq.0) THEN
        DO ina=1, Natoms
         Energy_min(ina) = UdiffX(1,ina)
         write(48,*)"i_atom, Energy_min", ina, Energy_min(ina)
        ENDDO       
        write(48,*)"Energy_ave",  SUM(Energy_min(1:Natoms))/Natoms
       ENDIF

       call MPI_BCAST(Energy_min(1:Natoms), Natoms, MPI_REAL8,0,MPI_COMM_WORLD,ierr)

       dE = 1.0d-1 * temp 
       E_cut = temp + SUM(Energy_min(1:Natoms))/Natoms 
       E_threshold = 1.5d-1 * temp  +  SUM(Energy_min(1:Natoms))/Natoms

       IF (myrank.eq.0) THEN
        write(48,*)"Energy_cut = ", E_cut
        write(48,*)"delta_Energy = ",dE
        write(48,*)"Energy_threshold = ",E_threshold
       ENDIF


       !!!!!!!!!!!!  \   ! Initial ENERGIES FOR GROUND STATE :  find start points

      ENDIF   ! \ mode=="Initialization"
      !!!!!!!!!!!!  \   ! Initialization for seeking of start points
 

      IF(nprocs.gt.1)Com_switch = .true.

      x_ini_r(1:Natoms) = x_pos(1,1:Natoms)
      y_ini_r(1:Natoms) = y_pos(1,1:Natoms)
      z_ini_r(1:Natoms) = z_pos(1,1:Natoms)

      i_ratio = 0

      ! initial delta parameter 
      ratio_ctrl = ratio + 0.02d0             ! 0.42d0

 200  continue           ! the label to reboot the simulation witha new delta

      i_ratio = i_ratio +1
 
      ratio_ctrl = ratio_ctrl - 0.02d0

      IF (myrank.eq.0) THEN
       write(*,*)"Change the delta parameter for ... times - i_ratio = ",i_ratio
       write(*,*)"Change the delta parameter - delta = ",ratio_ctrl
      ENDIF

      DO ina=1, Natoms 
       x_pos(1:Nwalker,ina) = x_ini_r(ina)
       y_pos(1:Nwalker,ina) = y_ini_r(ina)
       z_pos(1:Nwalker,ina) = z_ini_r(ina)
      ENDDO

      IF (i_ratio.gt.20) THEN
       stop
      ENDIF

      !  \initial delta parameter
     
      inv_G=0.1d0    ! inverse Gamma parameter in Langevin equation / term

      ratio_ini = ratio_ctrl

      call system_clock(c)
      c = (c + myrank)*(myrank+1)  ! different seeds for each proc  

      call init_genrand(c)
      IF(myrank.eq.0) THEN
       open(unit=11, file=fout, status="unknown")
       open(unit=12, file="atomic_coord_q.lammpstrj", status="unknown")

       open(unit=17, file="path.dat", status="unknown")
       open(unit=18, file="E_atoms.dat", status="unknown")

       open(unit=21, file="log_const.dat",status="unknown")

       IF(mode=="Initialization")THEN
        ! to save initial points: find start points  
             open(unit=49, file="Structures.lammpstrj", status="unknown")
             open(unit=50, file="Number_struct.dat", status="unknown")
        !     open(unit=30, file="Before_struct.lammpstrj", status="unknown")
             open(unit=31, file="After_heating_energy.dat", status="unknown")
       ENDIF  !\ to save initial points: find start points
       

       write(*,*)"initial ratio is ",ratio_ini
      ENDIF ! myrank .eq. 0

      ! Define the parameters of simulation 

      IF(mode=="Initialization")THEN
       Vswitch = .FALSE.
       N_restart = 2

       ELSE IF(mode=="Reaction")THEN
        Vswitch = .TRUE.
        N_restart = 10

         ELSE IF(mode=="Langevin")THEN
           Vswitch = .FALSE.
           N_restart = 1

             ELSE
               write(*,*)"error: mode = invalid"
               write(*,*)"mode could be Initialization or Reaction or Langevin, but mode = ",mode
               stop 
      ENDIF




      ! ITERATION FOR SEVERAL simulation with PULL OUT walkers at first steps ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
      restart_step = 0   ! initial number of step
      Nstep = floor(dble(Nstep)/dble(N_restart))

      IF (myrank.eq.0) THEN
       write(*,*)"N_restart= ",N_restart
       write(*,*)"Nstep= ",Nstep
      ENDIF

      DO i_restart =1, N_restart           !!!!!!!!!!!!!    RESTART    !!!!!!!!!!!!!!!

       ! the function of delta parameter could be changed for a constant, for example
       ratio = 0.48d0 - (0.48d0 - ratio_ini) * (N_restart - i_restart+1) / N_restart 
       !ratio = ratio_ini       

       IF (myrank.eq.0) THEN 
        write(*,*)"i_restart = ",i_restart
        write(*,*)"ratio = ",ratio
       ENDIF

       ! FOR RELAXATION - LAST CIRCLE
       IF ((i_restart .eq. N_restart) .AND. (mode.NE."Initialization")) THEN

        call restart_lmp()

        temp = tempFin    ! to increase the temperture at final stage

        Vswitch = .false.
        Nwalker = 1                     ! This time is for RELAXATION !!! 
        Com_switch = .false.            ! This time is for RELAXATION !!! 

        ! To write the positions of the saddle point
        IF((myrank.eq.0) .AND. (mode.EQ."Reaction")) THEN
         open(unit=50, file="saddle.dat", status="unknown")

         DO ina=1, Natoms
          write(50,'(1i10,3f16.8)')itype(ina),x_ave_q(ina), y_ave_q(ina), z_ave_q(ina)
         ENDDO
         close(50)

        ENDIF 

        ! To write all walkers for each processor (not more than 50) (one walker for one processor)

         IF ((myrank.lt.50) .AND. (mode.NE."Initialization")) THEN

          iw = myrank+100000
          write (istring, '(i6)') iw

          f_rank = "coordinates_rank_"//istring//".lammpstrj"

          open(unit=iw, file=f_rank, status="unknown")

          iw = myrank+200000
          f_rank = "E_atoms_"//istring//".dat"

          open(unit=iw, file=f_rank, status="unknown")

          Nstep = Nstep * N_restart
         ENDIF

       ENDIF      ! \(i_restart .eq. N_restart) .AND. (mode.NE."Initialization")




       ! FOR RELAXATION - LAST CIRCLE: find start points  
 
       IF ((i_restart .eq. N_restart) .AND. (mode.eq."Initialization")) THEN

        Vswitch = .false.
        Com_switch = .false.            ! This time is for RELAXATION !!! 

        IF (myrank.eq.0) THEN
           write(48,*)"Point i_restart = 2"
        ENDIF

        call calc_Udiff(Nwalker,Natoms,Nregions,Ncycles,&
        &x_pos(1:Nwalker,1:Natoms),y_pos(1:Nwalker,1:Natoms),z_pos(1:Nwalker,1:Natoms),Vdiff0(1:Nwalker))

        IF (myrank.eq.0) THEN
            write(48,*)"Calculate potential for each walker"
        ENDIF

        temp =  tempFin    ! 0.1d-2 ! the temperature for relaxation 

            ! We must to mark the walkers with E_cut-dE < E < E_cut+dE before
            ! relaxation at
            ! the SECOND STAGE and ignore the remaining walkers

        Vdiff0(1:Nwalker) = Vdiff0(1:Nwalker)/Natoms

        DO iw=1, Nwalker

             IF ((Vdiff0(iw) < E_cut+0.5d-1).AND.(Vdiff0(iw) > E_cut - 0.5d-1)) THEN
             !  IF (Vdiff0(iw) > E_cut) THEN
                    check_walker(iw) = .TRUE.  ! we mark this walker as a potential candidate 
             ENDIF
        ENDDO

        U_write_pr(:) = 0.0d0
        U_write_pr(n1:n2) = Vdiff0(1:Nwalker)

        call MPI_BARRIER(MPI_COMM_WORLD, ierr)

        call MPI_ALLREDUCE(U_write_pr,U_write, Nwalker*nprocs,MPI_REAL8, &
        &                     MPI_SUM, MPI_COMM_WORLD,ierr)

        ! To save all data from each processor to the zero rank (combine the data)

        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        IF (myrank.eq.0) THEN
            ! write(48,*)"Calculate E_cut"

            ! save energies of all walkers after heating
            DO iw=1, Nwalker*nprocs
               write(31,'(1i6,1e20.8e4)')iw,U_write(iw)
            ENDDO

            ! save the energies of choosen walkers
            DO iw=1, Nwalker
              IF (check_walker(iw)) THEN
                    write(48,*)iw,check_walker(iw),Vdiff0(iw)
              ENDIF
            ENDDO

        ENDIF ! \myrank=0
        time_second_stage = istep

     ENDIF      ! \(i_restart .eq. N_restart)

     !    \FOR RELAXATION - LAST CIRCLE: find start points 





       !      PULL OUT of walkers for second or more stage of restarting iteration


       IF ((i_restart .gt. 1) .AND. (mode=="Reaction")) THEN

        DO ina=1, Natoms
         x_pos(1,ina) = x_ave_q(ina)
         y_pos(1,ina) = y_ave_q(ina)
         z_pos(1,ina) = z_ave_q(ina)

         x_pos(2:Nwalker,ina) = x_pos(1,ina)
         y_pos(2:Nwalker,ina) = y_pos(1,ina)
         z_pos(2:Nwalker,ina) = z_pos(1,ina)
        ENDDO

       ENDIF




      ! Initialization for simulation circle 

      


       id_walker(:)= 0  
       UdiffX(:,:)=0d0
       UdiffY(:,:)=0d0
       UdiffZ(:,:)=0d0

       
       Vdiff0(:)=0d0
       !Vdiff0X(:,:)=0d0
       !Vdiff0Y(:,:)=0d0
       !Vdiff0Z(:,:)=0d0

       Vexp(:)=0d0
       Vexp_ave=0d0
       !Vexp(:,:)=0d0
       !Vexp_ave(:)=0d0
       V_ini=0d0

       F_ave=0d0

       VdiffX(:,:)=0d0
       VdiffY(:,:)=0d0
       VdiffZ(:,:)=0d0

       Vdiff2X(:,:)=0d0
       Vdiff2Y(:,:)=0d0
       Vdiff2Z(:,:)=0d0

       
      
       weight(:)=1d0
       F_func(:)=0d0
       mom1(:)=0d0
       mom1Natoms(:,:)=0d0
       prob(:)=0d0
       lc=0d0

       
       x_ave_q(:) = 0d0
       y_ave_q(:) = 0d0
       z_ave_q(:) = 0d0
       x_ave_p(:) = 0d0
       y_ave_p(:) = 0d0
       z_ave_p(:) = 0d0
       sq_force = 0d0
       
       stepXdist(:,:) = 0d0
       stepYdist(:,:) = 0d0
       stepZdist(:,:) = 0d0

       DO iw =1, Nwalker
        id_walker(iw) = iw
       ENDDO

       !! initial normalization of p

     

       !call  calc_rad_m2(Nwalker,Natoms,x_pos(1:Nwalker,1:Natoms),y_pos(1:Nwalker,1:Natoms),z_pos(1:Nwalker,1:Natoms),rad_m2(1:Nwalker,1:Natoms,1:Natoms))



       call calc_Udiff(Nwalker,Natoms,Nregions,Ncycles,&
        &          x_pos(1:Nwalker,1:Natoms),y_pos(1:Nwalker,1:Natoms),z_pos(1:Nwalker,1:Natoms), Vdiff0(1:Nwalker))

       IF (i_restart.eq.1) THEN

        E_start = Vdiff0(1)/Natoms   ! Initial energy of system

        IF (myrank.eq.0)THEN
         write(*,*)"Initial Energy = ",E_start

         write(*,*)"x_pos = ",x_pos(1,1:Natoms)
         write(*,*)"y_pos = ",y_pos(1,1:Natoms)
         write(*,*)"z_pos = ",z_pos(1,1:Natoms)
         write(*,*)"Natoms = ",Natoms
         write(*,*)"Nwalker = ",Nwalker
         write(*,*)"Vdiff0 = ",Vdiff0(1)
         write(*,*)"Nregions,Ncycles = ",Nregions,Ncycles
         write(*,*)"  " 
         ENDIF
       ENDIF
 

       Vdiff0(:) = (1d0 - ratio)*Vdiff0(:)

       ! prevent overflow problem
       !V_min = minval(Vdiff0(1:Nwalker))         ! prevent overflow problem 
       V_min_pr = minval(Vdiff0(1:Nwalker))         ! prevent overflow problem 
       call MPI_ALLREDUCE(V_min_pr, V_min, 1, MPI_REAL8, & 
        &                     MPI_MIN, MPI_COMM_WORLD,ierr)
       Vdiff0(1:Nwalker)=Vdiff0(1:Nwalker) - V_min

       ! prevent overflow problem
       Vexp(1:Nwalker) = dexp(-1d0 * Vdiff0(1:Nwalker)/temp)



       IF (myrank.eq.0) THEN
        write(*,*)"INITIALIZATION!!!!!!!!!"
        write(*,*)"V_min = ",V_min
        write(*,*)"Temperature of simulation T = ",temp
        write(*,*)"Vdiff_0(1) = "
        write(*,*)Vdiff0(1)+V_min
        write(*,*)"V_exp(1) = "
        write(*,*)Vexp(1)
        write(*,*)"  "

       ENDIF 
     
       weight(:)=1d0

       call calc_moment1(Vexp(1:Nwalker),Nwalker,weight,mom1(1:Nwalker),  &
        &               Vexp_ave_pr)
       call MPI_ALLREDUCE(Vexp_ave_pr, Vexp_ave, 1, MPI_REAL8, & 
        &                     MPI_SUM, MPI_COMM_WORLD,ierr)
       Vexp_ave = Vexp_ave/dble(nprocs)
       V_ini=V_min-temp*dlog(Vexp_ave)

       IF (myrank.eq.0) THEN

        write(*,*)"V_ini = ",V_ini

        write(*,*)"average values of Vexp = "
        write(*,*)ina,Vexp_ave
        write(*,*)"  "
       ENDIF


     !/END of initialization


     ! START of ITERATIONS     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

       B_coeff=sqrt(2d0*temp*inv_G)

       Voff = .false.
       DO istep =1, Nstep

        call system_clock(ct0, count_rate, count_max)

        weight(:)=1d0

        IF(Vswitch .and. (istep .le. 30) .and. (i_restart .gt. 1)) THEN
         Voff = .true.
         Vswitch = .false.
        ENDIF

        IF(Voff .and. istep.gt.30)THEN !! first 30 steps V is off
         Voff = .false.
         Vswitch=.true.
        ENDIF

        ! birth-death process1
        !
        !
        !! calculate factor P

        IF(Vswitch) THEN

!         call cpu_time(timer_s)
         call calc_Udiff(Nwalker,Natoms,Nregions,Ncycles,&
&                       x_pos(1:Nwalker,1:Natoms),y_pos(1:Nwalker,1:Natoms),z_pos(1:Nwalker,1:Natoms),&
&                    UdiffX(1:Nwalker,1:Natoms),UdiffY(1:Nwalker,1:Natoms),UdiffZ(1:Nwalker,1:Natoms),                                                     &
&                    Vdiff2X(1:Nwalker,1:Natoms),Vdiff2Y(1:Nwalker,1:Natoms),Vdiff2Z(1:Nwalker,1:Natoms) )

!        IF (myrank.eq.0) THEN
!          call cpu_time(timer_f)
!          write(*,*)"time of calculation is  ", timer_f - timer_s
!        ENDIF

!        write(*,*)"ID_atom  UdiffX UdiffY UdiffZ  Vdiff2X  Vdiff2Y  Vdiff2Z"
!        DO iw=1, Nwalker
!         write(*,*)"Number_of_walker = ", iw
!         DO ina=1, Natoms
!         write(*,*)ina,UdiffX(iw,ina),UdiffY(iw,ina),UdiffZ(iw,ina),Vdiff2X(iw,ina),Vdiff2Y(iw,ina),Vdiff2Z(iw,ina)
!         ENDDO
!        ENDDO

         VdiffX(:,:) = (1d0 - ratio)*UdiffX(:,:)
         VdiffY(:,:) = (1d0 - ratio)*UdiffY(:,:)
         VdiffZ(:,:) = (1d0 - ratio)*UdiffZ(:,:)

         Vdiff2X(:,:) = (1d0 - ratio)*Vdiff2X(:,:)
         Vdiff2Y(:,:) = (1d0 - ratio)*Vdiff2Y(:,:)
         Vdiff2Z(:,:) = (1d0 - ratio)*Vdiff2Z(:,:)


         F_func(1:Nwalker) =0d0
         DO ina=1, Natoms
          F_func(1:Nwalker) = F_func(1:Nwalker) + & 
&                          Vdiff2X(1:Nwalker,ina)+1d0/temp*VdiffX(1:Nwalker,ina)*(VdiffX(1:Nwalker,ina)-UdiffX(1:Nwalker,ina)) + &
&                          Vdiff2Y(1:Nwalker,ina)+1d0/temp*VdiffY(1:Nwalker,ina)*(VdiffY(1:Nwalker,ina)-UdiffY(1:Nwalker,ina)) + &
&                          Vdiff2Z(1:Nwalker,ina)+1d0/temp*VdiffZ(1:Nwalker,ina)*(VdiffZ(1:Nwalker,ina)-UdiffZ(1:Nwalker,ina))
         ENDDO

         !call calc_moment1(F_func(1:Nwalker),Nwalker,weight,mom1(1:Nwalker),F_ave)
         F_ave_pr = maxval(F_func)
         call MPI_ALLREDUCE(F_ave_pr, F_ave, 1, MPI_REAL8, &
&                         MPI_MAX, MPI_COMM_WORLD,ierr)
         mom1(1:Nwalker) = F_func(1:Nwalker) - F_ave
         prob(1:Nwalker)=exp(0.5d0*dt*inv_G*mom1(1:Nwalker))
         !prob(1:Nwalker)=exp(1.0d0*dt*inv_G*mom1(1:Nwalker)) !Langevin-->BD-->Langevin

         !!  exactly number-conserving algorithm (2018/09/05)
         !!! calculated accumulated probability 
         n1 = myrank*Nwalker + 1      
         n2 = n1 + Nwalker - 1      

         pr_accum(:) =0d0
         pr_accum_pr(:) =0d0
         pr_accum_pr(n1) = prob(1)

         DO iw = 2, Nwalker
          pr_accum_pr(n1+iw-1) = pr_accum_pr(n1 + iw -2)+prob(iw)
         ENDDO

         call MPI_ALLREDUCE(pr_accum_pr,pr_accum, Nwalker*nprocs, MPI_REAL8, & 
   &                     MPI_SUM, MPI_COMM_WORLD,ierr)
         DO iproc = 2, nprocs
          npr = (iproc-1)*Nwalker
          pr_accum((npr+1):(npr+Nwalker)) =                   &
   &       pr_accum((npr+1):(npr+Nwalker)) + pr_accum(npr)
         ENDDO

         pr_accum(:) = pr_accum(:)/pr_accum(Nwalker*nprocs)


         rand_p=genrand_real3()

         call MPI_BCAST(rand_p, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
 
         iw = 1
         myincr = 0
         IF(myrank.gt.0)THEN
          DO while (iw .le. Nwalker*nprocs)
           IF( ((dble(myrank*Nwalker)+rand_p )/dble(Nwalker*nprocs)).gt.pr_accum(iw) )THEN
            iw = iw + 1
            cycle
           ELSE
            myincr = iw -1
            exit
           ENDIF
          ENDDO
         ENDIF


         iw   =  1
         incr = myincr + 1

         DO while (iw .le. Nwalker)
          IF( ((dble(iw-1 + myrank*Nwalker) + rand_p )/dble(Nwalker*nprocs)).le.pr_accum(incr) )THEN
           id_walker(iw) = incr
           iw = iw + 1
          ELSE
           incr = incr + 1
          ENDIF
         ENDDO


        !! Pay attention: ALLREDUCE  can be used like (instead) ALLGATHER

         x_pos_tmp(:,:) = 0d0
         y_pos_tmp(:,:) = 0d0
         z_pos_tmp(:,:) = 0d0
         x_pos_tmp2(:,:) = 0d0
         y_pos_tmp2(:,:) = 0d0
         z_pos_tmp2(:,:) = 0d0

         x_pos_tmp(n1:n2, 1:Natoms) = x_pos(1:Nwalker, 1:Natoms)
         y_pos_tmp(n1:n2, 1:Natoms) = y_pos(1:Nwalker, 1:Natoms)
         z_pos_tmp(n1:n2, 1:Natoms) = z_pos(1:Nwalker, 1:Natoms)

         call MPI_ALLREDUCE(x_pos_tmp,x_pos_tmp2, Nwalker*nprocs*Natoms, MPI_REAL8, & 
   &                     MPI_SUM, MPI_COMM_WORLD,ierr)
         call MPI_ALLREDUCE(y_pos_tmp,y_pos_tmp2, Nwalker*nprocs*Natoms, MPI_REAL8, & 
   &                     MPI_SUM, MPI_COMM_WORLD,ierr)
         call MPI_ALLREDUCE(z_pos_tmp,z_pos_tmp2, Nwalker*nprocs*Natoms, MPI_REAL8, & 
   &                     MPI_SUM, MPI_COMM_WORLD,ierr)
        
        !!/allreduce can be used like allgather

         DO iw = 1, Nwalker
          x_pos(iw, 1:Natoms) = x_pos_tmp2(id_walker(iw), 1:Natoms)
          y_pos(iw, 1:Natoms) = y_pos_tmp2(id_walker(iw), 1:Natoms)
          z_pos(iw, 1:Natoms) = z_pos_tmp2(id_walker(iw), 1:Natoms)
         ENDDO

        !! /exactly number-conserving algorithm (2018/09/05)
        ENDIF ! Vswitch
        !/birth-death process 1

        ! Langevin step
        !
        !
        ! generate random noise
        call gen_noise(Nwalker,Natoms,dWX(1:Nwalker,1:Natoms),dWY(1:Nwalker,1:Natoms),dWZ(1:Nwalker,1:Natoms))

        !/generate random noise

        ! calculate coefficients
        call calc_Udiff(Nwalker,Natoms,Nregions,Ncycles,x_pos(1:Nwalker,1:Natoms),y_pos(1:Nwalker,1:Natoms),z_pos(1:Nwalker,1:Natoms),&
&                     UdiffX(1:Nwalker,1:Natoms),UdiffY(1:Nwalker,1:Natoms),UdiffZ(1:Nwalker,1:Natoms))
        VdiffX(:,:) = 0d0
        VdiffY(:,:) = 0d0
        VdiffZ(:,:) = 0d0
        IF(Vswitch) THEN
         VdiffX(:,:) = (1d0 - ratio)*UdiffX(:,:)
         VdiffY(:,:) = (1d0 - ratio)*UdiffY(:,:)
         VdiffZ(:,:) = (1d0 - ratio)*UdiffZ(:,:)
        ENDIF

        A_coeffX(1:Nwalker,1:Natoms) = inv_G*(-UdiffX(1:Nwalker,1:Natoms)+2d0*VdiffX(1:Nwalker,1:Natoms))
        A_coeffY(1:Nwalker,1:Natoms) = inv_G*(-UdiffY(1:Nwalker,1:Natoms)+2d0*VdiffY(1:Nwalker,1:Natoms))
        A_coeffZ(1:Nwalker,1:Natoms) = inv_G*(-UdiffZ(1:Nwalker,1:Natoms)+2d0*VdiffZ(1:Nwalker,1:Natoms))
        !!/calculate coefficients
        !
        !! increment
        x_pos(1:Nwalker,1:Natoms) = x_pos(1:Nwalker,1:Natoms)+A_coeffX(1:Nwalker,1:Natoms)*1.0d0*dt+B_coeff*sqrt(1.0d0*dt)*dWX(1:Nwalker,1:Natoms)
        y_pos(1:Nwalker,1:Natoms) = y_pos(1:Nwalker,1:Natoms)+A_coeffY(1:Nwalker,1:Natoms)*1.0d0*dt+B_coeff*sqrt(1.0d0*dt)*dWY(1:Nwalker,1:Natoms)
        z_pos(1:Nwalker,1:Natoms) = z_pos(1:Nwalker,1:Natoms)+A_coeffZ(1:Nwalker,1:Natoms)*1.0d0*dt+B_coeff*sqrt(1.0d0*dt)*dWZ(1:Nwalker,1:Natoms)
        !!/increment
        !
        !
        !/Langevin step 


        ! birth-death process2
        !
        !
        !! calculate factor P

        IF(Vswitch) THEN
         call calc_Udiff(Nwalker,Natoms,Nregions,Ncycles,&
&                       x_pos(1:Nwalker,1:Natoms),y_pos(1:Nwalker,1:Natoms),z_pos(1:Nwalker,1:Natoms),&
&                    UdiffX(1:Nwalker,1:Natoms),UdiffY(1:Nwalker,1:Natoms),UdiffZ(1:Nwalker,1:Natoms),                                                     &
&                    Vdiff2X(1:Nwalker,1:Natoms),Vdiff2Y(1:Nwalker,1:Natoms),Vdiff2Z(1:Nwalker,1:Natoms) )


         VdiffX(:,:) = (1d0 - ratio)*UdiffX(:,:)
         VdiffY(:,:) = (1d0 - ratio)*UdiffY(:,:)
         VdiffZ(:,:) = (1d0 - ratio)*UdiffZ(:,:)

         Vdiff2X(:,:) = (1d0 - ratio)*Vdiff2X(:,:)
         Vdiff2Y(:,:) = (1d0 - ratio)*Vdiff2Y(:,:)
         Vdiff2Z(:,:) = (1d0 - ratio)*Vdiff2Z(:,:)


         F_func(1:Nwalker) =0d0
         DO ina=1, Natoms
          F_func(1:Nwalker) = F_func(1:Nwalker) + & 
&                          Vdiff2X(1:Nwalker,ina)+1d0/temp*VdiffX(1:Nwalker,ina)*(VdiffX(1:Nwalker,ina)-UdiffX(1:Nwalker,ina)) + &
&                          Vdiff2Y(1:Nwalker,ina)+1d0/temp*VdiffY(1:Nwalker,ina)*(VdiffY(1:Nwalker,ina)-UdiffY(1:Nwalker,ina)) + &
&                          Vdiff2Z(1:Nwalker,ina)+1d0/temp*VdiffZ(1:Nwalker,ina)*(VdiffZ(1:Nwalker,ina)-UdiffZ(1:Nwalker,ina))
         ENDDO

         !call calc_moment1(F_func(1:Nwalker),Nwalker,weight,mom1(1:Nwalker),F_ave)
         F_ave_pr = maxval(F_func)
         call MPI_ALLREDUCE(F_ave_pr, F_ave, 1, MPI_REAL8, &
   &                     MPI_MAX, MPI_COMM_WORLD,ierr)
         mom1(1:Nwalker) = F_func(1:Nwalker) - F_ave
         prob(1:Nwalker)=exp(0.5d0*dt*inv_G*mom1(1:Nwalker))
         !prob(1:Nwalker)=exp(1.0d0*dt*inv_G*mom1(1:Nwalker)) !Langevin-->BD-->Langevin

         !!  exactly number-conserving algorithm (2018/09/05)
         !!! calculated accumulated probability 
         n1 = myrank*Nwalker + 1      
         n2 = n1 + Nwalker - 1      

         pr_accum(:) =0d0
         pr_accum_pr(:) =0d0
         pr_accum_pr(n1) = prob(1)

         DO iw = 2, Nwalker
          pr_accum_pr(n1+iw-1) = pr_accum_pr(n1 + iw -2)+prob(iw)
         ENDDO

         call MPI_ALLREDUCE(pr_accum_pr,pr_accum, Nwalker*nprocs, MPI_REAL8, & 
   &                     MPI_SUM, MPI_COMM_WORLD,ierr)
         DO iproc = 2, nprocs
          npr = (iproc-1)*Nwalker
          pr_accum((npr+1):(npr+Nwalker)) =                   &
   &       pr_accum((npr+1):(npr+Nwalker)) + pr_accum(npr)
         ENDDO

         pr_accum(:) = pr_accum(:)/pr_accum(Nwalker*nprocs)

        !!!/calculated accumulated probability       

         rand_p=genrand_real3()

         call MPI_BCAST(rand_p, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
 
         iw = 1
         myincr = 0
         IF(myrank.gt.0)THEN
          DO while (iw .le. Nwalker*nprocs)
           IF( ((dble(myrank*Nwalker)+rand_p )/dble(Nwalker*nprocs)).gt.pr_accum(iw) )THEN
            iw = iw + 1
            cycle
           ELSE
            myincr = iw -1
            exit
           ENDIF
          ENDDO
         ENDIF

         iw   =  1
         incr = myincr + 1

         DO while (iw .le. Nwalker)
          IF( ((dble(iw-1 + myrank*Nwalker) + rand_p )/dble(Nwalker*nprocs)).le.pr_accum(incr) )THEN
           id_walker(iw) = incr
           iw = iw + 1
          ELSE
           incr = incr + 1
          ENDIF
         ENDDO


        !! Pay attention: ALLREDUCE  can be used like (instead) ALLGATHER

         x_pos_tmp(:,:) = 0d0
         y_pos_tmp(:,:) = 0d0
         z_pos_tmp(:,:) = 0d0
         x_pos_tmp2(:,:) = 0d0
         y_pos_tmp2(:,:) = 0d0
         z_pos_tmp2(:,:) = 0d0

         x_pos_tmp(n1:n2, 1:Natoms) = x_pos(1:Nwalker, 1:Natoms)
         y_pos_tmp(n1:n2, 1:Natoms) = y_pos(1:Nwalker, 1:Natoms)
         z_pos_tmp(n1:n2, 1:Natoms) = z_pos(1:Nwalker, 1:Natoms)

         call MPI_ALLREDUCE(x_pos_tmp,x_pos_tmp2, Nwalker*nprocs*Natoms, MPI_REAL8, & 
   &                     MPI_SUM, MPI_COMM_WORLD,ierr)
         call MPI_ALLREDUCE(y_pos_tmp,y_pos_tmp2, Nwalker*nprocs*Natoms, MPI_REAL8, & 
   &                     MPI_SUM, MPI_COMM_WORLD,ierr)
         call MPI_ALLREDUCE(z_pos_tmp,z_pos_tmp2, Nwalker*nprocs*Natoms, MPI_REAL8, & 
   &                     MPI_SUM, MPI_COMM_WORLD,ierr)
        
        !!/Pay attention: ALLREDUCE  can be used like (instead) ALLGATHER

         DO iw = 1, Nwalker
          x_pos(iw, 1:Natoms) = x_pos_tmp2(id_walker(iw), 1:Natoms)
          y_pos(iw, 1:Natoms) = y_pos_tmp2(id_walker(iw), 1:Natoms)
          z_pos(iw, 1:Natoms) = z_pos_tmp2(id_walker(iw), 1:Natoms)
         ENDDO

        !! /exactly number-conserving algorithm (2018/09/05)
       !/birth-death process 2

        ENDIF ! Vswitch

     
        ! This part in order to save data of simulation    
        IF((mod(istep,steptowrite).eq.0).or.(istep.eq.1))THEN
      
         weight(1:Nwalker)=1d0
         lnweight_min=0d0
 

         call calc_moment1Natoms(x_pos(1:Nwalker,1:Natoms),Nwalker,Natoms,weight(1:Nwalker),mom1Natoms(1:Nwalker,1:Natoms),x_ave_pr_q(1:Natoms))
         call calc_moment1Natoms(y_pos(1:Nwalker,1:Natoms),Nwalker,Natoms,weight(1:Nwalker),mom1Natoms(1:Nwalker,1:Natoms),y_ave_pr_q(1:Natoms))
         call calc_moment1Natoms(z_pos(1:Nwalker,1:Natoms),Nwalker,Natoms,weight(1:Nwalker),mom1Natoms(1:Nwalker,1:Natoms),z_ave_pr_q(1:Natoms))

         call MPI_ALLREDUCE(x_ave_pr_q, x_ave_q, Natoms, MPI_REAL8, & 
   &                    MPI_SUM, MPI_COMM_WORLD,ierr)
         call MPI_ALLREDUCE(y_ave_pr_q, y_ave_q, Natoms, MPI_REAL8, & 
   &                    MPI_SUM, MPI_COMM_WORLD,ierr)
         call MPI_ALLREDUCE(z_ave_pr_q, z_ave_q, Natoms, MPI_REAL8, & 
   &                    MPI_SUM, MPI_COMM_WORLD,ierr)

         x_ave_q(:) = x_ave_q(:)/dble(nprocs)
         y_ave_q(:) = y_ave_q(:)/dble(nprocs)
         z_ave_q(:) = z_ave_q(:)/dble(nprocs)

         !!call calc_dist(x_pos(1:Nwalker,1:Natoms),y_pos(1:Nwalker,1:Natoms),z_pos(1:Nwalker,1:Natoms),weight(1:Nwalker,1:Natoms),Nwalker,Natoms,stepXdist,stepYdist,stepZdist,dist_Q,Nx,Ny,Nz,lnweight_min)

         !! calculate potential correction and exp(-bV)       
         IF(Vswitch)THEN

          call calc_Udiff(Nwalker,Natoms,Nregions,Ncycles,x_pos(1:Nwalker,1:Natoms),y_pos(1:Nwalker,1:Natoms),z_pos(1:Nwalker,1:Natoms),&
&                        Vdiff0(1:Nwalker))

          Vdiff0(:) = (1d0 - ratio)*Vdiff0(:)


          DO iw=1, Nwalker
           Vdiff0(iw)=Vdiff0(iw) - V_ini
          ENDDO 


          V_min_pr = Vdiff0(1)
          DO iw = 1, Nwalker
           IF (V_min_pr.gt.Vdiff0(iw)) THEN 
            V_min_pr = Vdiff0(iw)
           ENDIF
          ENDDO
          call MPI_ALLREDUCE(V_min_pr, V_min, 1, MPI_REAL8, & 
   &                    MPI_MIN, MPI_COMM_WORLD,ierr)
       
          lnweight_min=lc-V_min/temp
       
      

          DO iw = 1, Nwalker
           weight(iw)=exp((V_min-Vdiff0(iw))/temp)
          ENDDO

         ENDIF !Vswitch

         call calc_moment1Natoms(x_pos(1:Nwalker,1:Natoms),Nwalker,Natoms,weight(1:Nwalker),mom1Natoms(1:Nwalker,1:Natoms),x_ave_pr_p(1:Natoms))
         call calc_moment1Natoms(y_pos(1:Nwalker,1:Natoms),Nwalker,Natoms,weight(1:Nwalker),mom1Natoms(1:Nwalker,1:Natoms),y_ave_pr_p(1:Natoms))
         call calc_moment1Natoms(z_pos(1:Nwalker,1:Natoms),Nwalker,Natoms,weight(1:Nwalker),mom1Natoms(1:Nwalker,1:Natoms),z_ave_pr_p(1:Natoms))
         w_sum_pr = sum(weight(1:Nwalker))
         IF(Com_switch)THEN
          x_ave_pr_p(:)=x_ave_pr_p(:)*w_sum_pr
          y_ave_pr_p(:)=y_ave_pr_p(:)*w_sum_pr
          z_ave_pr_p(:)=z_ave_pr_p(:)*w_sum_pr
          call MPI_ALLREDUCE(w_sum_pr, w_sum, 1, MPI_REAL8, & 
   &                    MPI_SUM, MPI_COMM_WORLD,ierr)
          call MPI_ALLREDUCE(x_ave_pr_p, x_ave_p, Natoms, MPI_REAL8, & 
   &                    MPI_SUM, MPI_COMM_WORLD,ierr)
          call MPI_ALLREDUCE(y_ave_pr_p, y_ave_p, Natoms, MPI_REAL8, & 
   &                    MPI_SUM, MPI_COMM_WORLD,ierr)
          call MPI_ALLREDUCE(z_ave_pr_p, z_ave_p, Natoms, MPI_REAL8, & 
   &                    MPI_SUM, MPI_COMM_WORLD,ierr)
          x_ave_p(:) = x_ave_p(:)/w_sum
          y_ave_p(:) = y_ave_p(:)/w_sum
          z_ave_p(:) = z_ave_p(:)/w_sum
         ELSE ! Com_switch ==.false.
          x_ave_p(:) = x_ave_pr_p(:)
          y_ave_p(:) = y_ave_pr_p(:)
          z_ave_p(:) = z_ave_pr_p(:)
         ENDIF !Com_switch

         IF(myrank.eq.0)THEN
          call calc_Energy(Nwalker,Natoms,Nregions,Ncycles,x_ave_p(1:Natoms),y_ave_p(1:Natoms),z_ave_p(1:Natoms),Energy_p,E_atoms(1:Natoms))
          call calc_Energy(Nwalker,Natoms,Nregions,Ncycles,x_ave_q(1:Natoms),y_ave_q(1:Natoms),z_ave_q(1:Natoms),Energy_q,E_atoms(1:Natoms))

!          call calc_Force(Nwalker,Natoms,Nregions, Ncycles,x_ave_q(1:Natoms),y_ave_q(1:Natoms),z_ave_q(1:Natoms),&
!      &                   forceX(1:Natoms), forceY(1:Natoms), forceZ(1:Natoms))

          call calc_Force(1,Natoms,Nregions, Ncycles,x_ave_q(1:Natoms),y_ave_q(1:Natoms),z_ave_q(1:Natoms),&
      &                   forceX(1:Natoms), forceY(1:Natoms), forceZ(1:Natoms))

          sq_force =0d0
          DO ina=1, Natoms
           sq_force = sq_force + forceX(ina)*forceX(ina) &
      &                        + forceY(ina)*forceY(ina) &
      &                        + forceZ(ina)*forceZ(ina)
          ENDDO
        
         ENDIF
         call MPI_BCAST(Energy_q, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
         call MPI_BCAST(Energy_p, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
         
         ! TO check the delta parameter for Rection mode and to change and to reboot the simulation if it is nessesary
         IF ((mode.eq."Reaction").and.(i_restart.lt.N_restart).and.(istep.gt.400)) THEN
          IF (E_start .gt. Energy_q .and. Vswitch ) THEN 
           Nstep = Nstep*N_restart
           GOTO 200
          ENDIF
         ENDIF
         !  \check the delta parameter 


         ! Write the reaction path with Q (and P) functions + Energy of system  
         !ina = 1
         !write(17,'(1i6,3e20.8e4)')istep,dist_Q(max_d(1,ina),max_d(2,ina),max_d(3,ina),ina),dist_P(max_d_corr(1,ina),max_d_corr(2,ina),max_d_corr(3,ina),ina),Energy
         IF(myrank.eq.0.and.(i_restart.lt.N_restart)) THEN 
          write(17,'(1i6,3e20.8e4,i10)')istep+restart_step,Energy_q, Energy_p,sqrt(sq_force),Nwalker*nprocs

          write(18,'(1i6,1e20.8e4)', advance='no')istep+restart_step,Energy_q

          ! To write the data only for first 20 atoms or less
          IF(Natoms.lt.21)THEN
           DO ina = 1, Natoms-1
            write(18,'(1e20.8e4)', advance='no')E_atoms(ina)
           ENDDO
           write(18,'(1e20.8e4)')E_atoms(Natoms)
            ELSE 
            DO ina = 1, 19
             write(18,'(1e20.8e4)', advance='no')E_atoms(ina)
            ENDDO
            write(18,'(1e20.8e4)')E_atoms(20)
          ENDIF 
         ENDIF !!/(myrank.eq.0.and.(i_restart.lt.N_restart))


         ! To write Energy and coordinates during relaxation for first 50 walkers

         IF((i_restart.EQ.N_restart) .AND. (mode.NE."Initialization") .AND. (myrank.LT.50))  THEN

          x_ave_q(1:Natoms) = x_pos(1,1:Natoms) 
          y_ave_q(1:Natoms) = y_pos(1,1:Natoms)
          z_ave_q(1:Natoms) = z_pos(1,1:Natoms)

          call calc_Energy(Nwalker,Natoms,Nregions,Ncycles,x_ave_q(1:Natoms),y_ave_q(1:Natoms),z_ave_q(1:Natoms),Energy_q,E_atoms(1:Natoms))

          iw=200000+myrank
          write(iw,'(1i6,1e20.8e4)', advance='no')istep+restart_step,Energy_q

          DO ina = 1, Natoms-1
           write(iw,'(1e20.8e4)', advance='no')E_atoms(ina)
          ENDDO
          write(iw,'(1e20.8e4)')E_atoms(Natoms)

          iw=100000+myrank
!          call write_positions(Natoms,x_min,x_max,y_min,y_max,z_min,z_max,iw,istep+restart_step,itype,x_ave_q(1:Natoms),y_ave_q(1:Natoms),z_ave_q(1:Natoms))
          call write_positions(Natoms,xlo,xhi,ylo,yhi,zlo,zhi,xy, xz, yz,bounds,  &
      &          iw,istep+restart_step,itype,x_ave_q(1:Natoms),y_ave_q(1:Natoms),z_ave_q(1:Natoms))

         ENDIF ! i_restart.eq.N_restart .AND. (mode.NE."Initialization") .AND. (myrank.LT.50)


        ! write average atomic coordinates for Q distribution in the LAMMPS format !!!!!!!
        IF(myrank.eq.0)THEN
         call write_positions(Natoms,xlo,xhi,ylo,yhi,zlo,zhi,xy, xz, yz,bounds, &
    &         12,istep+restart_step,itype,x_ave_q(1:Natoms),y_ave_q(1:Natoms),z_ave_q(1:Natoms))
!         call write_positions(Natoms,x_min,x_max,y_min,y_max,z_min,z_max,12,istep+restart_step,itype,x_ave_q(1:Natoms),y_ave_q(1:Natoms),z_ave_q(1:Natoms))
        ENDIF

       ENDIF !steptowrite                       


       ! Check the ENERGY at each time step during relaxation : find start points 

       IF ((mode.eq."Initialization").and.(i_restart.eq.2)) THEN

         ! Check and save the NEW STRUCTURES!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! calc  calc_Udiff0_atom 
         call calc_Udiff(Nwalker,Natoms,Nregions,Ncycles,    &
         &              x_pos(1:Nwalker,1:Natoms), y_pos(1:Nwalker,1:Natoms),z_pos(1:Nwalker,1:Natoms), UdiffX(1:Nwalker,1:Natoms))

         DO jw = 1, Nwalker

           IF ((SUM(UdiffX(jw,1:Natoms))/Natoms < E_threshold).AND.(check_walker(jw)).AND.(.NOT.check_crossing(jw))) THEN

             ! we check the structure only 1 time so 
             check_crossing(jw) = .TRUE.

             N_struct = n1 - 1 + jw
             E_struct(N_struct,1:Natoms) = UdiffX(jw,1:Natoms)

             time_pr(N_struct) = Nstep + istep - time_second_stage

             Energy(1:Natoms) = UdiffX(jw,1:Natoms)
             call get_unique(Natoms,Energy_min(1:Natoms),Energy(1:Natoms),UniqueNumber)

             number_struct(N_struct) = UniqueNumber       ! Unique Number of structure
             ! write(*,*)"myrank, N_struct, UniqueNumberi, time_pr",myrank, N_struct, UniqueNumber, time_pr(N_struct)
             x_pos_str(N_struct,1:Natoms) = x_pos(jw,1:Natoms)
             y_pos_str(N_struct,1:Natoms) = y_pos(jw,1:Natoms)
             z_pos_str(N_struct,1:Natoms) = z_pos(jw,1:Natoms)

           ENDIF   ! \ENDIF for check crossing

         ENDDO  ! \ENDDO for check all walker at each time step!!!

      ENDIF   ! \ (mode.eq."Initialization").and. (i_restart.EQ.2)

      !!! \Check the ENERGY at each time step during relaxation : find start points

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      ENDDO              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   END OF ITERATION     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      restart_step = restart_step + Nstep


     ENDDO              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   END OF RESTART LOOP  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     ! To save initial coordinates for mode Initialization
     IF(mode.eq."Initialization")THEN 
           call save_init_data(myrank,nprocs,n1,n2,check_crossing(1:Nwalker),check_walker(1:Nwalker),           &
      &                         Nwalker,Natoms,                                                                 &
      &                         time_pr(1:Nwalker*nprocs),                                                      &
      &                         x_pos_str(1:Nwalker*nprocs,1:Natoms),y_pos_str(1:Nwalker*nprocs,1:Natoms),      &
      &                         z_pos_str(1:Nwalker*nprocs,1:Natoms),E_struct(1:Nwalker*nprocs,1:Natoms),       &
      &                         number_struct(1:Nwalker*nprocs),xlo,xhi,ylo,yhi,zlo,zhi,xy,xz,yz,               &
      &                         bounds(1:3),itype(1:Natoms))
     ENDIF
    
!        IF (myrank.eq.0) THEN
!          call cpu_time(timer_f)
!          write(*,*)"time of all calculations is  ", timer_f - timer_s
!        ENDIF

 
     ! subroutine to close files for each rank   
     call Fclose(myrank,mode)

   end subroutine solve
  !-----------

   subroutine finalizempi()
   use val_mpi, only : ierr

   call MPI_FINALIZE(ierr)

   end subroutine finalizempi

  !---------------------------------------------
  ! To unite data from each rank to zero and save to files
   subroutine save_init_data(myrank,nprocs,n1,n2,check_crossing,check_walker,   &
      &                         Nwalker,Natoms,time_pr,                         &
      &                         x_pos_str,y_pos_str,z_pos_str,E_struct,         &
      &                         number_struct,xlo,xhi,ylo,yhi,zlo,zhi,xy,xz,yz, &
      &                         bounds,itype)
  !-----------------------
   
   use val_mpi, only : ierr
   integer, intent(in)  :: myrank,nprocs,n1,n2,Nwalker,Natoms,itype(Natoms)
   character(1), intent(in)  :: bounds(3)
   real(8), intent(in)  :: E_struct(Nwalker*nprocs,Natoms),xlo,xhi,ylo,yhi,zlo,zhi,xy,xz,yz
   real(8), intent(in)  :: x_pos_str(Nwalker*nprocs,Natoms),y_pos_str(Nwalker*nprocs,Natoms),z_pos_str(Nwalker*nprocs,Natoms)
   integer, intent(in)  :: number_struct(Nwalker*nprocs),time_pr(Nwalker*nprocs)

   logical, intent(in)  :: check_walker(Nwalker),check_crossing(Nwalker)

   logical              :: check_str(Nwalker*nprocs),check_pr_str(Nwalker*nprocs),check_walker_pr(Nwalker)
   integer              :: time(Nwalker*nprocs),iw,jw,istep,ina,istruct,N_struct,UniqueNumber,number_pr_struct(Nwalker*nprocs), number_struct_out(Nwalker*nprocs)
   real(8)              :: U_write(Nwalker*nprocs),E_pr_struct(Nwalker*nprocs,Natoms),U_write_pr(Nwalker*nprocs)
   real(8)              :: x_pos_tmp2(Nwalker*nprocs, Natoms),x_ave_q(Natoms)
   real(8)              :: y_pos_tmp2(Nwalker*nprocs, Natoms),y_ave_q(Natoms)
   real(8)              :: z_pos_tmp2(Nwalker*nprocs, Natoms),z_ave_q(Natoms)
   real(8)              :: buf_str(Nwalker*nprocs,Natoms)
   real(8)              :: x_pos_str_out(Nwalker*nprocs,Natoms),y_pos_str_out(Nwalker*nprocs,Natoms)
   real(8)              :: z_pos_str_out(Nwalker*nprocs,Natoms),E_struct_out(Nwalker*nprocs,Natoms)


     ! Save the data and choose the best walkers: find start points
     IF (myrank.eq.0) THEN
          write(48,*)"End of iteration and restart loop"
     ENDIF

     check_walker_pr(:) = .FALSE.
     check_str(:) = .FALSE.
     check_pr_str(:) = .FALSE. 
     DO jw = 1, Nwalker
        check_walker_pr(jw) = check_walker(jw) .AND. check_crossing(jw) 
!        check_pr_str(n1-1+jw) = check_walker(jw) .AND. check_crossing(jw)
     ENDDO

     ! To save all data from each processor to the zero rank (combine the data)
     call MPI_BARRIER(MPI_COMM_WORLD, ierr)

     ! check_walker_pr(1:Nwalker) = check_pr_str(n1:n2)
     call MPI_GATHER(check_walker_pr(1:Nwalker),Nwalker,MPI_LOGICAL,check_str,Nwalker,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)

     call MPI_BARRIER(MPI_COMM_WORLD, ierr)

     time(:) = 0
     call MPI_ALLREDUCE(time_pr,time, Nwalker*nprocs,MPI_INTEGER, &
     &                     MPI_SUM, MPI_COMM_WORLD,ierr)


     x_pos_tmp2(:,:) = 0.0d0
     y_pos_tmp2(:,:) = 0.0d0
     z_pos_tmp2(:,:) = 0.0d0

     x_pos_tmp2(n1:n2,1:Natoms) = x_pos_str(n1:n2,1:Natoms)
     x_pos_str_out(:,:) = 0.0d0

     call MPI_BARRIER(MPI_COMM_WORLD, ierr)

     call MPI_ALLREDUCE(x_pos_tmp2,x_pos_str_out,Nwalker*nprocs*Natoms,MPI_REAL8, &
     &                     MPI_SUM, MPI_COMM_WORLD,ierr)


     y_pos_tmp2(n1:n2,1:Natoms) = y_pos_str(n1:n2,1:Natoms)
     y_pos_str_out(:,:) = 0.0d0

     call MPI_BARRIER(MPI_COMM_WORLD, ierr)

     call MPI_ALLREDUCE(y_pos_tmp2,y_pos_str_out,Nwalker*nprocs*Natoms,MPI_REAL8,&
     &                     MPI_SUM, MPI_COMM_WORLD,ierr)

     z_pos_tmp2(n1:n2,1:Natoms) = z_pos_str(n1:n2,1:Natoms)
     z_pos_str_out(:,:) = 0.0d0

     call MPI_BARRIER(MPI_COMM_WORLD, ierr)

     call MPI_ALLREDUCE(z_pos_tmp2,z_pos_str_out,Nwalker*nprocs*Natoms,MPI_REAL8,&
     &                     MPI_SUM, MPI_COMM_WORLD,ierr)

     E_pr_struct(:,:) = 0.0d0
     E_pr_struct(n1:n2,1:Natoms) = E_struct(n1:n2,1:Natoms)
     E_struct_out(:,:) = 0.0d0
     call MPI_BARRIER(MPI_COMM_WORLD, ierr)

     call MPI_ALLREDUCE(E_pr_struct,E_struct_out,Nwalker*nprocs*Natoms,MPI_REAL8,&
     &                     MPI_SUM, MPI_COMM_WORLD,ierr)

     number_pr_struct(:) = 0
     number_pr_struct(n1:n2) = number_struct(n1:n2)
     number_struct_out(:) = 0

     call MPI_BARRIER(MPI_COMM_WORLD, ierr)

     call MPI_ALLREDUCE(number_pr_struct,number_struct_out,Nwalker*nprocs,MPI_INTEGER,&
     &                     MPI_SUM, MPI_COMM_WORLD,ierr)


     !IF (myrank.eq.0) THEN
     !  DO jw = 1, Nwalker*nprocs
     !   write(*,*)"jw, check_str, time, number_struct_out ",jw, check_str(jw), time(jw), number_struct_out(jw)
     !  ENDDO
     !ENDIF

     !write(*,*)myrank,"To save data"

     IF (myrank.eq.0) THEN

         write(48,*)"To save data"
         write(50,*)"To save data"

         ! Save ALL data from structures

         open(unit=90,file="ALL_STRUCTURES.dat",status="unknown")

         DO jw=1,Nwalker*nprocs

            ! to write all separated walkers
            ! write(*,*)"check_str", check_str(jw)
            IF (check_str(jw)) THEN
                  ! Save ALL data from structures
                  write(90,*)jw,time(jw),check_str(jw),number_struct_out(jw),x_pos_str_out(jw,1),y_pos_str_out(jw,1),z_pos_str_out(jw,1),SUM(E_struct_out(jw,1:Natoms))/Natoms

                  x_ave_q(1:Natoms) = x_pos_str_out(jw,1:Natoms)
                  y_ave_q(1:Natoms) = y_pos_str_out(jw,1:Natoms)
                  z_ave_q(1:Natoms) = z_pos_str_out(jw,1:Natoms)

                  call write_positions(Natoms,xlo,xhi,ylo,yhi,zlo,zhi,xy, xz,yz,bounds,&
                  &                     49,jw,itype,x_ave_q(1:Natoms),y_ave_q(1:Natoms),z_ave_q(1:Natoms))

                  write(48,'(1i10,1e20.8e4)',advance='no')jw,SUM(E_struct_out(jw,1:Natoms))/Natoms

                  DO ina=1,Natoms-1
                        write(48,'(1e20.8e4)', advance='no')E_struct_out(jw,ina)
                  ENDDO
                  write(48,'(1e20.8e4)')E_struct_out(jw,Natoms)
                  write(50,'(3i24)')jw,number_struct_out(jw),time(jw)
            ENDIF  !/ (check_str(jw)) 
         ENDDO !\jw=1,Nwalker*nprocs

     ENDIF  ! (myrank.eq.0)


     IF (myrank.eq.0) THEN

        open(unit=38,file="Chosen_walk_Energy.dat",status="unknown")
        open(unit=39,file="Chosen_walk_Structures.lammpstrj",status="unknown")
        open(unit=40,file="Chosen_walk_Number",status="unknown")

        iw = 0
        istep = 1  ! To use istep as a variable for choosen the slowest walker

        write(*,*)"To save the Choosen walkers"

        DO jw=1,Nwalker*nprocs

                ! to write ONLY choosen walkers, for example only with unique numbers  102XX,
                ! which means that 1-st and 2-nd atoms have a largest deviations 

                IF (check_str(jw).AND.(floor(number_struct_out(jw)/1.0d2).EQ.102))THEN     !  1st and 2nd atoms have largest deviations

                      x_ave_q(1:Natoms) = x_pos_str_out(jw,1:Natoms)
                      y_ave_q(1:Natoms) = y_pos_str_out(jw,1:Natoms)
                      z_ave_q(1:Natoms) = z_pos_str_out(jw,1:Natoms)

                      call write_positions(Natoms,xlo,xhi,ylo,yhi,zlo,zhi,xy,xz,yz,bounds,&
                      &                         39,jw,itype,x_ave_q(1:Natoms),y_ave_q(1:Natoms),z_ave_q(1:Natoms))

                      write(38,'(1i10,1e20.8e4)',advance='no')jw,SUM(E_struct_out(jw,1:Natoms))/Natoms
                      DO ina=1,Natoms-1
                              write(38,'(1e20.8e4)',advance='no')E_struct_out(jw,ina)
                      ENDDO

                      write(38,'(1e20.8e4)')E_struct_out(jw,Natoms)
                      write(40,'(3i12)')jw,number_struct_out(jw),time(jw)

                      IF (istep .LT. time(jw)) THEN
                          iw = jw
                          istep = time(jw)
                      ENDIF
                ENDIF  !/ (check_str(jw)) 
        ENDDO !\jw=1,Nwalker*nprocs


        close(38)
        close(39)
        close(40)
        close(89)
        close(90)

        open(unit=68,file="Chosen_BEST_walk.dat",status="unknown")

        write(68,*)"Energy"

        IF (iw.NE.0) THEN
           
           write(68,'(1i10,1e20.8e4)',advance='no')iw,SUM(E_struct_out(iw,1:Natoms))/Natoms

           DO ina=1,Natoms-1
                 write(68,'(1e20.8e4)', advance='no')E_struct_out(iw,ina)
           ENDDO

           write(68,'(1e20.8e4)')E_struct_out(iw,Natoms)

           write(68,*)"Number and time"

           write(68,'(3i6)')iw,number_struct_out(iw),time(iw)
           write(68,*)"Coordinates"

           DO ina=1,Natoms
                write(68,'(1i10,3e20.8e4)')itype(ina),x_pos_str_out(iw,ina),y_pos_str_out(iw,ina),z_pos_str_out(iw,ina)
           ENDDO
        ENDIF    ! \(iw.NE.0)

        close(68)

      ENDIF  ! \(myrank.eq.0)

      ! \Save the data and choose the best walkers: find start points

   end subroutine save_init_data



  !---------------------------------------------
  ! To close Files for each rank
   subroutine Fclose(myrank,mode)
  !-----------------------

   implicit none
   integer, intent(in)   :: myrank
   character, intent(in) :: mode
   integer :: iw

    IF(myrank.eq.0) THEN
     close(11)
     close(12)
     close(17)
     close(18)
     close(21)
    ENDIF

     IF((myrank.LT.50) .AND. (mode.NE."Initialization")) THEN
       iw = 100000 + myrank
       close(iw)
       iw = 200000 + myrank
       close(iw)
     ENDIF

     IF (myrank.eq.0) THEN
        write(*,*)"End of simulation."
     ENDIF

   end subroutine Fclose


   !!! Unique number : find start points
   !------------------------
   subroutine get_unique(Natoms,Energy_min,Energy,UniqueNumber)
   !-----------------------

   implicit none
   integer, intent(in) :: Natoms
   real(8), intent(in) :: Energy(Natoms),Energy_min(Natoms)
   integer, intent(out) :: UniqueNumber
   real(8) :: delta_E(Natoms)
   integer :: ina

   UniqueNumber = 0
   delta_E(1:Natoms) = Energy(1:Natoms) - Energy_min(1:Natoms)

   ! First number 

   ina = MAXLOC(delta_E(1:Natoms),1)

   UniqueNumber = UniqueNumber + 10000*ina

   delta_E(ina) = minval(delta_E(1:Natoms))

   ! Second number
   ina = MAXLOC(delta_E(1:Natoms),1)

   UniqueNumber = UniqueNumber + 100*ina

   delta_E(ina) = minval(delta_E(1:Natoms))

   ! Third number
   ina = MAXLOC(delta_E(1:Natoms),1)

   UniqueNumber = UniqueNumber + ina

   delta_E(ina) = minval(delta_E(1:Natoms))

    end subroutine get_unique
   !------------------------  

   !------------------------
   subroutine calc_rad_m2(Nwalker,Natoms,x_pos,y_pos,z_pos,rad_minus2) 
   !------------------------
   !This routine calculates the distances between atoms for each walker
   !and keep the values as 1/r^2  
   ! 
   implicit none
   integer, intent(in) :: Nwalker
   integer, intent(in) :: Natoms

   real(8), intent(in) :: x_pos(Nwalker,Natoms)
   real(8), intent(in) :: y_pos(Nwalker,Natoms)
   real(8), intent(in) :: z_pos(Nwalker,Natoms)

   real(8), intent(out) :: rad_minus2(Nwalker,Natoms,Natoms)
   
   integer :: iw,ina_i,ina_j

!  rad_minus2(:,:,:) = 0d0
   DO iw=1, Nwalker
      
      DO ina_i=1, Natoms-1
       DO ina_j=ina_i+1, Natoms
 
            rad_minus2(iw,ina_i,ina_j) = 1d0 / ((x_pos(iw,ina_i)-x_pos(iw,ina_j))**2d0 + (y_pos(iw,ina_i)-y_pos(iw,ina_j))**2d0 + (z_pos(iw,ina_i)-z_pos(iw,ina_j))**2d0)
       !     print *,rad_minus2(iw,ina_i,ina_j),iw,ina_i,ina_j 

       ENDDO
      ENDDO
   ENDDO

  ! print *,rad_minus2(1,1,2)

    end subroutine calc_rad_m2
   !------------------------  

   !------------------------
   subroutine calc_Energy(Nwalker,Natoms,Nregions,Ncycles,x_pos,y_pos,z_pos,Energy,E_atoms)
   !------------------------
   !
   !This routine calculates the potential Energy for ALL atoms or for atomic system
   !
   use val_mpi_lmp, only :  ptr
   use potentials, only: calc_Udiff
   implicit none
   
   integer, intent(in) :: Nwalker  
   integer, intent(in) :: Nregions
   integer, intent(in) :: Ncycles
   !! Nwalker, Nregions, Ncycles are needed for pudding of arrays, so that lammps_scatter_atoms
   !! works properly.
   integer, intent(in) :: Natoms 

   real(8), intent(in) :: x_pos(Natoms)
   real(8), intent(in) :: y_pos(Natoms)
   real(8), intent(in) :: z_pos(Natoms)

   real(8), intent(out) :: Energy
   real(8), intent(out) :: E_atoms(1:Natoms)

   real(8) :: x_pos_tmp(1,Natoms)
   real(8) :: y_pos_tmp(1,Natoms)
   real(8) :: z_pos_tmp(1,Natoms)
   real(8) :: Energy_tmp(1)
   real(8) :: pe_tmp(1,Natoms)

   integer :: ina



   Energy = 0d0
   !! lammps ver. not yet implemented
   E_atoms(1:Natoms) = 0d0
   !!/lammps ver. not yet implemented

   DO ina = 1, Natoms
    x_pos_tmp(1,ina) = x_pos(ina)
    y_pos_tmp(1,ina) = y_pos(ina)
    z_pos_tmp(1,ina) = z_pos(ina)
   ENDDO

   call calc_Udiff(1, Natoms,Nregions,1, x_pos_tmp, y_pos_tmp, z_pos_tmp, Energy_tmp)
   Energy = Energy_tmp(1)
   call calc_Udiff(1, Natoms,Nregions,1, x_pos_tmp, y_pos_tmp, z_pos_tmp, pe_tmp)
   E_atoms(1:Natoms) = pe_tmp(1,1:Natoms)

   !write(*,*)"E_atoms=", E_atoms(:)
   !stop

   !DO ina = 1, Natoms
   ! pos(1,3*(ina-1)+1) = x_pos(ina)
   ! pos(1,3*(ina-1)+2) = y_pos(ina)
   ! pos(1,3*(ina-1)+3) = z_pos(ina)
   !ENDDO
   !call lammps_scatter_atoms(ptr, 'x', 1, pos(1,:))
   !call lammps_command(ptr, 'run 0', 5)
   !call lammps_extract_compute(ptr, 'thermo_pe', 9, Energy)

   Energy = Energy/dble(Natoms)

   end subroutine calc_Energy
   !------------------------
   !

   !------------------------
   subroutine calc_Force(Nwalker,Natoms,Nregions,Ncycles,x_pos, y_pos, z_pos, forceX, forceY, forceZ)
   !------------------------
   !
   !This routine calculates the potential Energy for ALL atoms or for atomic system
   !
   use val_mpi_lmp, only :  ptr
   use potentials, only: calc_Udiff
   implicit none
   
   integer, intent(in) :: Nwalker  
   integer, intent(in) :: Nregions
   integer, intent(in) :: Ncycles
   !! Nwalker, Nregions, Ncycles are needed for pudding of arrays, so that lammps_scatter_atoms
   !! works properly.
   integer, intent(in) :: Natoms 

   real(8), intent(in) :: x_pos(Natoms)
   real(8), intent(in) :: y_pos(Natoms)
   real(8), intent(in) :: z_pos(Natoms)

   real(8), intent(out) :: forceX(Natoms)
   real(8), intent(out) :: forceY(Natoms)
   real(8), intent(out) :: forceZ(Natoms)

   real(8) :: x_pos_tmp(Nwalker,Natoms)
   real(8) :: y_pos_tmp(Nwalker,Natoms)
   real(8) :: z_pos_tmp(Nwalker,Natoms)
   real(8) :: forceX_tmp(Nwalker,Natoms)
   real(8) :: forceY_tmp(Nwalker,Natoms)
   real(8) :: forceZ_tmp(Nwalker,Natoms)

   integer :: ina




   DO ina = 1, Natoms
    x_pos_tmp(1:Nwalker,ina) = x_pos(ina)
    y_pos_tmp(1:Nwalker,ina) = y_pos(ina)
    z_pos_tmp(1:Nwalker,ina) = z_pos(ina)
   ENDDO

   call calc_Udiff(Nwalker, Natoms,Nregions,Ncycles, x_pos_tmp, y_pos_tmp, z_pos_tmp, forceX_tmp, forceY_tmp, forceZ_tmp)
   forceX(1:Natoms) = forceX_tmp(1,1:Natoms)
   forceY(1:Natoms) = forceY_tmp(1,1:Natoms)
   forceZ(1:Natoms) = forceZ_tmp(1,1:Natoms)

   end subroutine calc_Force

   !------------------------
   subroutine calc_dist(x_pos,y_pos,z_pos,weight,Nwalker,Natoms,stepXdist,stepYdist,stepZdist,dist_Q,Nx,Ny,Nz,lnweight_min)
   !------------------------
   !
   !This routine calculates the distribution of Q or P 
   !
   !
   implicit none
   real(8), intent(in) :: x_pos(Nwalker,Natoms)
   real(8), intent(in) :: y_pos(Nwalker,Natoms)
   real(8), intent(in) :: z_pos(Nwalker,Natoms)
   real(8), intent(in) :: weight(Nwalker,Natoms)
   real(8), intent(in) :: lnweight_min
   integer, intent(in) :: Nwalker
   integer, intent(in) :: Natoms
   real(8), intent(out) :: stepXdist(Nx,Natoms)
   real(8), intent(out) :: stepYdist(Ny,Natoms)
   real(8), intent(out) :: stepZdist(Nz,Natoms)

   real(8), intent(out) :: dist_Q(1:Nx,1:Ny,1:Nz,1:Natoms)
   integer, intent(in) :: Nx
   integer, intent(in) :: Ny
   integer, intent(in) :: Nz
   integer :: iw, ix, iy, iz, ina, iat

   real(8) :: xmax(1:Natoms), ymax(1:Natoms), zmax(1:Natoms), xmin(1:Natoms), ymin(1:Natoms), zmin(1:Natoms), xdelta(1:Natoms), ydelta(1:Natoms), zdelta(1:Natoms) 
   real(8) :: diffx, diffy, diffz

   real(8) :: degree
   real(8), parameter :: delta = 0.01d0
   real(8), parameter :: pi = 3.14159265358979d0


   stepXdist(:,:) = 0d0
   stepYdist(:,:) = 0d0
   stepZdist(:,:) = 0d0

  DO ina=1,Natoms
   xmax(ina)=maxval(x_pos(:,ina))+delta
   ymax(ina)=maxval(y_pos(:,ina))+delta
   zmax(ina)=maxval(z_pos(:,ina))+delta

   xmin(ina)=minval(x_pos(:,ina))-delta
   ymin(ina)=minval(y_pos(:,ina))-delta
   zmin(ina)=minval(z_pos(:,ina))-delta


   xdelta(ina) = (xmax(ina)-xmin(ina))/dble(Nx-1)
   ydelta(ina) = (ymax(ina)-ymin(ina))/dble(Ny-1)
   zdelta(ina) = (zmax(ina)-zmin(ina))/dble(Nz-1)

   stepXdist(1,ina)=xmin(ina)
   stepYdist(1,ina)=ymin(ina)
   stepZdist(1,ina)=zmin(ina)
  ENDDO

  DO ina=1,Natoms
   DO ix = 2, Nx
    stepXdist(ix,ina) = xmin(ina) + xdelta(ina)*(ix-1)
   ENDDO
   DO iy = 2, Ny
    stepYdist(iy,ina) = ymin(ina) + ydelta(ina)*(iy-1)
   ENDDO
   DO iz = 2, Nz
    stepZdist(iz,ina) = zmin(ina) + zdelta(ina)*(iz-1)
   ENDDO
  ENDDO

  ! Calculation for 3D case
  ! During calculation we need to get log(Q), that is why it is impossible that Q=0
  dist_Q(:,:,:,:)=1d-300
  
  DO iat=1,Natoms 
   DO iz=1, Nz
    DO iy=1, Ny
     DO ix=1, Nx

      DO iw=1, Nwalker

      degree = 0d0

        diffx = stepXdist(ix,iat)-x_pos(iw,iat)
        diffy = stepYdist(iy,iat)-y_pos(iw,iat)
        diffz = stepZdist(iz,iat)-z_pos(iw,iat)
        degree = degree + (diffx)**2d0 + (diffy)**2d0 + (diffz)**2d0 


!      dist_Q(ix,iy,iz) =     dist_Q(ix,iy,iz) +  weight(iw)*exp(-(diffx/delta)**2d0)*exp(-(diffy/delta)**2d0)*exp(-(diffz/delta)**2d0)
!      dist_Q(ix,iy,iz,ina) = dist_Q(ix,iy,iz,ina) +    delta*delta*log(weight(iw,ina)) - (diffx(ina)**2d0 + diffy(ina)**2d0 +   diffz(ina)**2d0)


      dist_Q(ix,iy,iz,iat) =   dist_Q(ix,iy,iz,iat) + weight(iw,iat) * exp(-degree/delta/delta)

      ENDDO
     ENDDO
    ENDDO
   ENDDO
  ENDDO

!  dist_Q(:,:,:,:)= dist_Q(:,:,:,:)/delta/delta - log(dble(Nwalker)) - log(delta * sqrt(pi))*3d0*dble(Natoms)


! dist_Q(:,:,:,:)= dist_Q(:,:,:,:)/((delta * sqrt(pi))**3d0)

! dist_Q(:,:,:,:)= log(dist_Q(:,:,:,:)) + lnweight_min

dist_Q(:,:,:,:)= log(dist_Q(:,:,:,:)) - 3d0 * log(delta * sqrt(pi))+ lnweight_min

   end subroutine calc_dist

   !------------------------
   subroutine calc_moment1(func,Nwalker,weight,mom1,ave)
   !------------------------
   !
   !This routine calculates the first moment of x 
   ! x - <x>
   !
   implicit none
   integer, intent(in) :: Nwalker
   real(8), intent(in) :: func(Nwalker)
   real(8), intent(in) :: weight(Nwalker)  
   real(8), intent(out) :: mom1(Nwalker)
   real(8), intent(out) :: ave

   real(8) :: mean
   integer :: iw
   real(8) :: w_sum

   mean = 0d0
   w_sum = 0d0

   DO iw=1, Nwalker
    mean = mean + func(iw)*weight(iw)
   ENDDO

   DO iw=1, Nwalker
    w_sum = w_sum + weight(iw)
   ENDDO
   !mean = mean/dble(Nwalker)
   mean = mean/w_sum
   mom1(:)= func(:)-mean
   ave = mean
   end subroutine calc_moment1


   !------------------------
   subroutine calc_moment1Natoms(func,Nwalker,Natoms,weight,mom1Natoms,aveNatoms)
   !------------------------
   !
   !This routine calculates the first moment of x for N atoms 
   ! x - <x>
   
   implicit none
   integer, intent(in) :: Nwalker,Natoms
   real(8), intent(in) :: func(Nwalker,Natoms)
   real(8), intent(in) :: weight(Nwalker)  
   real(8), intent(out) :: mom1Natoms(Nwalker,Natoms)
   real(8), intent(out) :: aveNatoms(Natoms)

   real(8) :: mean(Natoms)
   integer :: iw,ina

   real(8) :: w_sum

   mean(:) = 0d0
   w_sum = 0d0

   DO ina=1, Natoms
    DO iw=1, Nwalker
     mean(ina) = mean(ina) + func(iw,ina)*weight(iw)
    ENDDO
   ENDDO

   DO iw=1, Nwalker
    w_sum = w_sum + weight(iw)
   ENDDO
   DO ina=1, Natoms
    !mean(ina) = mean(ina)/dble(Nwalker)
    mean(ina) = mean(ina)/w_sum
    mom1Natoms(:,ina)= func(:,ina)-mean(ina)
    aveNatoms(ina) = mean(ina)
   ENDDO

   end subroutine calc_moment1Natoms


end module routines





   !------------------------
   !subroutine write_positions(Natoms,x_min,x_max,y_min,y_max,z_min,z_max,n_file,istep,itype,coorX_max,coorY_max,coorZ_max)
   subroutine write_positions(Natoms,xlo,xhi,ylo,yhi,zlo,zhi,xy, xz, yz,bounds, n_file,istep,itype,coorX_max,coorY_max,coorZ_max)
   !------------------------
   !This subroutine write to the file the coordinates of atomic positions in the
   !format of
   !LAMMPS or OVITO  
   ! 
   use bistable1d_val, only: b_vec
   use val_mpi, only:  myrank
   implicit none
   integer, intent(in) :: Natoms
   integer, intent(in) :: itype(Natoms)
   integer, intent(in) :: n_file,istep

   real(8), intent(in) :: xlo, xhi, ylo, yhi, zlo, zhi
   real(8), intent(in) :: xy, xz, yz
   character(1), intent(in) :: bounds(3)

   real(8), intent(in) :: coorX_max(Natoms)
   real(8), intent(in) :: coorY_max(Natoms)
   real(8), intent(in) :: coorZ_max(Natoms)

   real(8) :: coorX(Natoms)
   real(8) :: coorY(Natoms)
   real(8) :: coorZ(Natoms)

   real(8) :: coorX_max_tmp(1,Natoms)
   real(8) :: coorY_max_tmp(1,Natoms)
   real(8) :: coorZ_max_tmp(1,Natoms)
   real(8) :: coorX_max_mod(1,Natoms)
   real(8) :: coorY_max_mod(1,Natoms)
   real(8) :: coorZ_max_mod(1,Natoms)

   real(8) :: x_min,x_max,y_min,y_max,z_min,z_max
!   real(8), intent(out) :: rad_minus2(Nwalker,Natoms,Natoms)

   integer :: ina
   character(2) :: bds(3)

    !! 
   !! Call periodic treatment
   coorX_max_tmp(1,:) = coorX_max(:)
   coorY_max_tmp(1,:) = coorY_max(:)
   coorZ_max_tmp(1,:) = coorZ_max(:)
   call lattice(1, Natoms, coorX_max_tmp, coorY_max_tmp, coorZ_max_tmp, &
     &                coorX_max_mod, coorY_max_mod, coorZ_max_mod)
   !!/Call periodic treatment

! write atomic coordinates ONLY for FIST  WALKER !!! in the LAMMPS format
! !!!!!!!

    write(n_file,'(A14)')"ITEM: TIMESTEP"
    write(n_file,'(i0)')istep

    write(n_file,'(A21)')"ITEM: NUMBER OF ATOMS"

    write(n_file,'(i0)')Natoms

!    write(n_file,'(A25)')"ITEM: BOX BOUNDS ss ss ss"
    bds(:)="ff"
    IF(bounds(1).eq."p") bds(1)="pp"
    IF(bounds(2).eq."p") bds(2)="pp"
    IF(bounds(3).eq."p") bds(3)="pp"
    write(n_file,'(A25,3a4)')"ITEM: BOX BOUNDS xy xz yz", bds(1:3)

    !! see LAMMPS manual 8.3.2 "Triclinic simulation boxes"
    x_min = xlo + dmin1(0d0, xy, xz, xy+xz)
    x_max = xhi + dmax1(0d0, xy, xz, xy+xz)
    y_min = ylo + dmin1(0d0, yz)
    y_max = yhi + dmax1(0d0, yz)
    z_min = zlo
    z_max = zhi
    write(n_file,'(3f14.6)')x_min, x_max, xy
    write(n_file,'(3f14.6)')y_min, y_max, xz
    write(n_file,'(3f14.6)')z_min, z_max, yz
!    write(n_file,'(f0.4," ",f0.4)')x_min,x_max
!    write(n_file,'(f0.4," ",f0.4)')y_min,y_max
!    write(n_file,'(f0.4," ",f0.4)')z_min,z_max

!!!! Recalculation of coordinates - for box coordinates  !!!



    DO ina = 1, Natoms

      coorX(ina) = (coorX_max_mod(1,ina) - xlo)*b_vec(1,1)  &
 &                +(coorY_max_mod(1,ina) - ylo)*b_vec(1,2)  &
 &                +(coorZ_max_mod(1,ina) - zlo)*b_vec(1,3)

      coorY(ina) = (coorX_max_mod(1,ina) - xlo)*b_vec(2,1)  &
 &                +(coorY_max_mod(1,ina) - ylo)*b_vec(2,2)  &
 &                +(coorZ_max_mod(1,ina) - zlo)*b_vec(2,3)

      coorZ(ina) = (coorX_max_mod(1,ina) - xlo)*b_vec(3,1)  &
 &                +(coorY_max_mod(1,ina) - ylo)*b_vec(3,2)  &
 &                +(coorZ_max_mod(1,ina) - zlo)*b_vec(3,3)

!     coorX(ina) = (coorX_max_mod(1,ina) - x_min) / (x_max - x_min)

!     coorY(ina) = (coorY_max_mod(1,ina) - y_min) / (y_max - y_min)

!     coorZ(ina) = (coorZ_max_mod(1,ina) - z_min) / (z_max - z_min)

    ENDDO

    write(n_file,'(A28)')"ITEM: ATOMS id type xs ys zs"
    DO ina = 1, Natoms

      !write(n_file,'(i0," ",i0," ",f0.6," ",f0.6," ",f0.6)')ina,itype(ina),coorX(ina),coorY(ina),coorZ(ina)
      write(n_file,'(i6,i6,3f14.6)')ina,itype(ina),coorX(ina),coorY(ina),coorZ(ina)
      !write(n_file,'(i6,i6,3f14.6)')ina,itype(ina),coorX_max_mod(1,ina),coorY_max_mod(1,ina),coorZ_max_mod(1,ina)
      !write(n_file,'(i6,i6,3f14.6)')ina,itype(ina),coorX_max(ina),coorY_max(ina),coorZ_max(ina)

    ENDDO

   end subroutine write_positions

   subroutine recvec(a_vec, b_vec)
   ! calculate the reciprocal vector
   ! in the normalization a_i *b_j = \delta_ij
   implicit none

   real(8), intent(in)  :: a_vec(3,3)
   real(8), intent(out) :: b_vec(3,3)

   real(8) vol
    vol = a_vec(1,1)*(a_vec(2,2)*a_vec(3,3)-a_vec(2,3)*a_vec(3,2))  &
     &       +a_vec(1,2)*(a_vec(2,3)*a_vec(3,1)-a_vec(2,1)*a_vec(3,3))  &
     &       +a_vec(1,3)*(a_vec(2,1)*a_vec(3,2)-a_vec(2,2)*a_vec(3,1))
    vol = dabs(vol)
    b_vec(1,1) = a_vec(2,2)*a_vec(3,3)-a_vec(2,3)*a_vec(3,2)
    b_vec(1,2) = a_vec(2,3)*a_vec(3,1)-a_vec(2,1)*a_vec(3,3)
    b_vec(1,3) = a_vec(2,1)*a_vec(3,2)-a_vec(2,2)*a_vec(3,1)
    b_vec(2,1) = a_vec(3,2)*a_vec(1,3)-a_vec(3,3)*a_vec(1,2)
    b_vec(2,2) = a_vec(3,3)*a_vec(1,1)-a_vec(3,1)*a_vec(1,3)
    b_vec(2,3) = a_vec(3,1)*a_vec(1,2)-a_vec(3,2)*a_vec(1,1)
    b_vec(3,1) = a_vec(1,2)*a_vec(2,3)-a_vec(1,3)*a_vec(2,2)
    b_vec(3,2) = a_vec(1,3)*a_vec(2,1)-a_vec(1,1)*a_vec(2,3)
    b_vec(3,3) = a_vec(1,1)*a_vec(2,2)-a_vec(1,2)*a_vec(2,1)
    b_vec(:,:) = b_vec(:,:)/vol
       
   end subroutine

program main
  use routines
  implicit none

  call initmpi()
  !write(*,*)"initmpi is finished"
  call stdin()
  !write(*,*)"stdin is finished"  
  call gen_walker()
  !write(*,*)"gen_walker is finished"
  call init_lmp(.false.)
  !write(*,*)"init_lmp is finished"
  call solve()
  !write(*,*)"solve is finished"
  call finalize_lmp()
  !write(*,*)"finilize_lmp is finished"
  call finalizempi()
  !write(*,*)"Simulation is finished"

  !call write_data()
end program


