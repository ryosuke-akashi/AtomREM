!
!        THIS IS 3N dimensional extension for N ARGON atoms !!!!!!!
!
! This program solves the Ito SDE for the Smoluchowski equation with
! the artificially heated walker method
! 3D poential U(x,y,z)= 4*epsilon * [sigma^12/r^12 - sigma^6/r^6]            [where r^2=(x^2) + (y^2) + (z^2)] 
! The log scale for calculation of P and Q functions is used !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The reaction path is drawn by saving of maxima of Q and P distributions
! With reaction path the atomic coordinates also are saved  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module bistable1d_val
  implicit none

  integer, save :: &
  & Nstep,         & ! number of steps
  & Nwalker,       & ! number of walkers
  & Natoms,        & ! number of atoms 
  & steptowrite      ! number of steps in order to write the path and other data


  real(8), save :: &
  & temp,          & ! temperature [Kelvin]
  & factor,        & ! heating factor T_eff = temp*factor
  & x_ini(2),         & ! initial x1 and x2 positions for mode=initial
  & y_ini(2),         & ! initial y1 and y2 positions for mode=initial
  & z_ini(2),         & ! initial z1 and z2 positions for mode=initial
  & dt,               & ! timestep
  & ratio,            & ! delta parameter in equation
  & ratio_ctrl          ! delta parameter in equation !! AKASHI

  character(256), save :: &
  & mode,          &  ! initial mode 
  & &                 ! 'initial':delta function
  & &                 ! 'load'   :read from file
  & fin,           &  ! walker input filename
  & fout              ! walker output filename

  logical, save :: &
  & Vswitch , Dist_switch, Debag          ! on or off V, saving of distribution and debagging info

  real(8), allocatable, save :: &
  & x_pos(:,:), &
  & y_pos(:,:), &
  & z_pos(:,:)

end module bistable1d_val

module val_mpi
  implicit none
  integer :: ierr, myrank, nprocs
end module val_mpi
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

  !-----------
   subroutine stdin()
   use bistable1d_val, only : Nstep, Nwalker, Natoms, temp, factor,steptowrite,x_ini, y_ini, z_ini, dt, ratio, &
  &                           mode, fin, fout,Vswitch, Dist_switch,Debag 

   use val_mpi, only: ierr, myrank, nprocs
   namelist /input/ Nstep, Nwalker, Natoms, temp, factor,steptowrite,x_ini, y_ini, z_ini, &
  &                 dt, ratio, mode, fin, fout, Vswitch, Dist_switch,Debag 

   IF(myrank.eq.0)THEN
    open(unit=10,file="param_3N.in",status="unknown")
    read(unit=10,nml=input,err=100)
    close(10)
    IF(mod(Nwalker,nprocs).ne.0)THEN
     write(*,*)"Error: Nwalker must be multiple of nprocs"
     stop
    ENDIF
    Nwalker = Nwalker/nprocs
!   read(*,input)

!   write(*,*)"x_ini=",x_ini(:,Natoms)
!   write(*,*)"y_ini=",y_ini(:,Natoms)
!   write(*,*)"z_ini=",z_ini(:,Natoms)a
   ENDIF
   call MPI_BCAST(Nstep   , 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ierr)
   call MPI_BCAST(Nwalker , 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ierr)
   call MPI_BCAST(Natoms  , 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ierr)
   call MPI_BCAST(temp    , 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
   call MPI_BCAST(factor  , 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
   call MPI_BCAST(steptowrite   , 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ierr)
   call MPI_BCAST(x_ini   , 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
   call MPI_BCAST(y_ini   , 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
   call MPI_BCAST(z_ini   , 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
   call MPI_BCAST(dt      , 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
   call MPI_BCAST(ratio   , 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
   call MPI_BCAST(mode    , 1, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
   call MPI_BCAST(Vswitch , 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
   call MPI_BCAST(Dist_switch , 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
   call MPI_BCAST(Debag , 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
   return
100 write(*,*) "ERROR: reading namelist file"
   stop
   end subroutine stdin
  !-----------

  !-----------
   subroutine gen_walker()

     use mt19937
     use bistable1d_val, only: Nwalker, Natoms, x_pos, y_pos, z_pos, mode, x_ini, y_ini, z_ini, fin
     use val_mpi, only: ierr, myrank, nprocs

     integer iw,ina
     real(8) :: ranx
     real(8) :: rany
     real(8) :: ranz

     ALLOCATE(x_pos(Nwalker,Natoms))
     ALLOCATE(y_pos(Nwalker,Natoms))
     ALLOCATE(z_pos(Nwalker,Natoms))

     IF(myrank.eq.0)THEN

      IF(mode=="initial")THEN
       x_pos(:,:)= x_ini(1)
       y_pos(:,:)= y_ini(1)
       z_pos(:,:)= z_ini(1)

      ELSE IF(mode=="random")THEN
      DO ina=1, Natoms
       DO iw=1, Nwalker
        ranx=genrand_real1()
        rany=genrand_real1()
        ranz=genrand_real1()
        x_pos(iw,ina)= ranx*(x_ini(2)-x_ini(1))+x_ini(1)
        y_pos(iw,ina)= rany*(y_ini(2)-y_ini(1))+y_ini(1)
        z_pos(iw,ina)= ranz*(z_ini(2)-z_ini(1))+z_ini(1)
       ENDDO
      ENDDO     

      ELSE IF(mode=="load")THEN

       open(unit=10,file=fin, status="unknown")
       rewind(10)
       DO ina=1, Natoms
        read(10,*)x_pos(1,ina), y_pos(1,ina), z_pos(1,ina)
        x_pos(2:Nwalker,ina) = x_pos(1,ina)
        y_pos(2:Nwalker,ina) = y_pos(1,ina)
        z_pos(2:Nwalker,ina) = z_pos(1,ina)
       ENDDO
       close(10)

      ELSE
       write(*,*)"error: mode = invalid"
      ENDIF

     ENDIF
     call MPI_BCAST(x_pos   , Nwalker*Natoms, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
     call MPI_BCAST(y_pos   , Nwalker*Natoms, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
     call MPI_BCAST(z_pos   , Nwalker*Natoms, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

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

   use bistable1d_val, only: Nstep, Nwalker,  Natoms, x_pos, y_pos, z_pos, factor, steptowrite, dt, ratio,ratio_ctrl, fout, temp, x_ini, y_ini, z_ini, Vswitch, Dist_switch, Debag
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

     !! Akashi (2018/09/04)
     real(8) :: Vdiff0(Nwalker)

     !! Akashi (2018/09/04)
     real(8) :: Vexp(Nwalker)
     !real(8) :: Vexp(Nwalker,Natoms)

     real(8) :: VdiffX(Nwalker,Natoms)
     real(8) :: VdiffY(Nwalker,Natoms)
     real(8) :: VdiffZ(Nwalker,Natoms)

     real(8) :: Vdiff2X(Nwalker,Natoms)
     real(8) :: Vdiff2Y(Nwalker,Natoms)
     real(8) :: Vdiff2Z(Nwalker,Natoms)
     !! Akashi modified the shape
     !real(8) :: weight(Nwalker,Natoms)
     real(8) :: weight(Nwalker)
     !!/Akashi modified the shape
     real(8) :: F_func(Nwalker)

     real(8) :: A_coeffX(Nwalker,Natoms)
     real(8) :: A_coeffY(Nwalker,Natoms)
     real(8) :: A_coeffZ(Nwalker,Natoms)

     real(8) :: prob(Nwalker)
     integer :: id_walker(Nwalker)

     integer :: max_d(3,Natoms)
     integer :: max_d_corr(3,Natoms)

     !! Akashi added
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
     !!/Akashi added
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
     !! Akashi (2018/09/04)
     real(8) :: Vexp_ave
     !real(8) :: Vexp_ave(1:Natoms)

     real(8) :: displacement(1:Natoms),x_ini_r(1:Natoms), y_ini_r(1:Natoms), z_ini_r(1:Natoms)

     real(8) :: lnweight_min
     real(8) :: B_coeff
     !! Akashi added
     real(8) :: Energy_q
     real(8) :: Energy_p
     real(8) :: E_atoms(Natoms)
     real(8) :: sq_force(Nwalker)
     real(8) :: min_abs_force
     !!/Akashi added
     real(8) :: V_min
     real(8) :: V_ini
     integer :: Nw_tmp
     integer :: incr
     real(8) :: inv_f
     real(8) :: F_ave
     real(8) :: pr
     !! Akashi added
     real(8) :: pr_accum(Nwalker*nprocs)
     !!/Akashi added
     real(8) :: rand_p
     real(8) :: temp0
     real(8) :: inv_G
     real(8) :: lc
     real(8) :: x_min, x_max
     real(8) :: y_min, y_max
     real(8) :: z_min, z_max
     real(8) :: maxtimestep
     integer :: istep
     integer :: iw,jw
     integer :: ix
     integer :: iy
     integer :: iz
     integer :: c,ina, i_ratio 
     integer :: i_restart,N_restart,restart_step
     real(8) :: ratio_ini,E_start

     integer :: atomictype

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
     character(3)  :: istring
     character(32) :: f_rank

     !! mpi shared values with pudding
     real(8) :: x_pos_tmp2(Nwalker*nprocs, Natoms)
     real(8) :: y_pos_tmp2(Nwalker*nprocs, Natoms)
     real(8) :: z_pos_tmp2(Nwalker*nprocs, Natoms)
     !!/mpi shared values with pudding
     integer :: ct0, ct1, count_rate, count_max
     integer :: ct2, ct3, ct4

      x_min = minval(x_pos(1,1:Natoms))
      x_max = maxval(x_pos(1,1:Natoms))

      y_min = minval(y_pos(1,1:Natoms))
      y_max = maxval(y_pos(1,1:Natoms))

      z_min = minval(z_pos(1,1:Natoms))
      z_max = maxval(z_pos(1,1:Natoms))

      x_min = x_min - (x_max - x_min)
      x_max = x_max + (x_max - x_min) / 2.0d0


      y_min = y_min - (y_max - y_min)
      y_max = y_max + (y_max - y_min) / 2.0d0

      z_min = z_min - (z_max - z_min)
      z_max = z_max + (z_max - z_min) / 2.0d0



     IF(nprocs.gt.1)Com_switch = .true.

     x_ini_r(1:Natoms) = x_pos(1,1:Natoms)
     y_ini_r(1:Natoms) = y_pos(1,1:Natoms)
     z_ini_r(1:Natoms) = z_pos(1,1:Natoms)

! label to restart with a new delta = ratio 


     i_ratio = 0

!     ratio = 0.42d0   ! initial delta parameter
!!    AKASHI
      ratio_ctrl = 0.42d0
!!   /AKASHI
200  continue

     i_ratio = i_ratio +1
!!    AKASHI
     ratio_ctrl = ratio_ctrl - 0.02d0
!!   /AKASHI

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




     maxtimestep = dt * 5.0d3
     temp0 = temp*1000d0
     inv_G=0.1d0
     inv_f=1d0/factor

!!  AKASHI
!     ratio_ini = ratio
     ratio_ini = ratio_ctrl
!!  AKASHI

     call system_clock(c)
     c = (c + myrank)*(myrank+1)  ! different seeds for each proc  
!     write(*,*)"proc, rand seed=",myrank,c
     call init_genrand(c)
     IF(myrank.eq.0) THEN
      open(unit=11, file=fout, status="unknown")
      open(unit=12, file="atomic_coord_q.lammpstrj", status="unknown")
!      open(unit=13, file="atomic_coord_1st_walker.lammpstrj", status="unknown")
      !! Akashi added
!      open(unit=15, file="atomic_coord_p.lammpstrj", status="unknown")
      open(unit=17, file="path.dat", status="unknown")
      open(unit=18, file="E_atoms.dat", status="unknown")

!     open(unit=13, file="out_max.dat", status="unknown")
      open(unit=21, file="log_const.dat",status="unknown")

      write(*,*)"x_pos(0) y_pos(0) z_pos(0)", x_pos(1,:), y_pos(1,:), z_pos(1,:)

      write(*,*)"initial ratio is ",ratio_ini
     ENDIF ! myrank .eq. 0

     ! ITERATION FOR SEVERAL simulation with PULL OUT walkers at first steps ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     N_restart = 10
     restart_step = 0   ! initial number of step
!     B_coeff = Nstep / N_restart

!     Nstep = floor(B_coeff)
     Nstep = floor(dble(Nstep)/dble(N_restart))

IF (myrank.eq.0) THEN
     write(*,*)"N_restart= ",N_restart
     write(*,*)"Nstep= ",Nstep
ENDIF

     DO i_restart =1, N_restart           !!!!!!!!!!!!!    RESTART    !!!!!!!!!!!!!!!


!!    AKASHI
!      ratio = 0.4d0 - (0.4d0 - ratio_ini) * (N_restart - i_restart) / N_restart 
      ratio = 0.4d0 - (0.4d0 - ratio_ini) * (N_restart - i_restart+1) / N_restart 
!!/   AKASHI

IF (myrank.eq.0) THEN 
      write(*,*)"i_restart = ",i_restart
      write(*,*)"ratio = ",ratio
ENDIF

      IF (i_restart .eq. (N_restart-1)) THEN

       Nstep = Nstep * 2
      ENDIF



! FOR RELAXATION - LAST CIRCLE

      IF (i_restart .eq. N_restart) THEN

       Vswitch = .false.
       !Nwalker = nprocs              ! This time is for RELAXATION !!! 
       Nwalker = 1              ! This time is for RELAXATION !!! 
       Com_switch = .false.            ! This time is for RELAXATION !!! 

! To write the positions of the saddle point
     IF(myrank.eq.0) THEN
      open(unit=50, file="saddle.dat", status="unknown")

       DO ina=1, Natoms
        write(50,'(3f16.8)')x_ave_q(ina), y_ave_q(ina), z_ave_q(ina)
       ENDDO
      close(50)

     ENDIF 

! To write all walkers for each processor (one walker for one processor)

     iw = myrank+100
     write (istring, '(i3)') iw

     f_rank = "coordinates_rank_"//istring//".lammpstrj"

     !open(unit=iw, file="COORDINATES/"//f_rank, status="unknown")
     open(unit=iw, file=f_rank, status="unknown")

     iw = myrank+200
     f_rank = "E_atoms_"//istring//".dat"

     !open(unit=iw, file="ENERGY/"//f_rank, status="unknown")
     open(unit=iw, file=f_rank, status="unknown")

       Nstep = Nstep * N_restart

      ENDIF      ! \(i_restart .eq. N_restart)

!      PULL OUT of walkers for second or more stage of restarting iteration
!      !!!!!!!!

      IF (i_restart .gt. 1) THEN

       DO ina=1, Natoms
        x_pos(1,ina) = x_ave_q(ina)
        y_pos(1,ina) = y_ave_q(ina)
        z_pos(1,ina) = z_ave_q(ina)

        x_pos(2:Nwalker,ina) = x_pos(1,ina)
        y_pos(2:Nwalker,ina) = y_pos(1,ina)
        z_pos(2:Nwalker,ina) = z_pos(1,ina)
       ENDDO

      ENDIF


      ! initialize

      atomictype = 1


      id_walker(:)= 0  
      UdiffX(:,:)=0d0
      UdiffY(:,:)=0d0
      UdiffZ(:,:)=0d0

      !! Akashi (2018/09/04)
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

      !! Akashi modified the shape of weight
      !weight(:,:)=1d0
      weight(:)=1d0
      F_func(:)=0d0
      mom1(:)=0d0
      mom1Natoms(:,:)=0d0
      prob(:)=0d0
      lc=0d0

      !! Akashi added
      x_ave_q(:) = 0d0
      y_ave_q(:) = 0d0
      z_ave_q(:) = 0d0
      x_ave_p(:) = 0d0
      y_ave_p(:) = 0d0
      z_ave_p(:) = 0d0
      sq_force(:) = 0d0
      !!/Akashi added
      stepXdist(:,:) = 0d0
      stepYdist(:,:) = 0d0
      stepZdist(:,:) = 0d0

      DO iw =1, Nwalker
       id_walker(iw) = iw
      ENDDO

      !! initial normalization of p

     

      !call  calc_rad_m2(Nwalker,Natoms,x_pos(1:Nwalker,1:Natoms),y_pos(1:Nwalker,1:Natoms),z_pos(1:Nwalker,1:Natoms),rad_m2(1:Nwalker,1:Natoms,1:Natoms))

      call calc_Udiff(Nwalker,Natoms,x_pos(1:Nwalker,1:Natoms),y_pos(1:Nwalker,1:Natoms),z_pos(1:Nwalker,1:Natoms), Vdiff0(1:Nwalker))

      IF (i_restart.eq.1) THEN

       E_start = Vdiff0(1)/Natoms   ! Initial energy of system
!IF (myrank.eq.0) THEN
!       write(*,*)"Initial Energy = ",E_start
!ENDIF

      ENDIF

! TO check the delta parameter and to change if it is nessesary
      IF (i_restart.eq.2) THEN

       IF (E_start .gt. Vdiff0(1)/Natoms) THEN  

        Nstep = Nstep*N_restart

        GOTO 200
       ENDIF

      ENDIF
!  \check the delta parameter 
 

      Vdiff0(:) = (1d0 - ratio)*Vdiff0(:)

      ! prevent overflow problem
      !V_min = minval(Vdiff0(1:Nwalker))         ! prevent overflow problem 
      V_min_pr = minval(Vdiff0(1:Nwalker))         ! prevent overflow problem 
      call MPI_ALLREDUCE(V_min_pr, V_min, 1, MPI_REAL8, & 
   &                     MPI_MIN, MPI_COMM_WORLD,ierr)
      Vdiff0(1:Nwalker)=Vdiff0(1:Nwalker) - V_min
      ! prevent overflow problem

!  20 continue 

!    DO iw=1,Nwalker
!     DO ina=1, Natoms

      Vexp(1:Nwalker) = dexp(-1d0 * Vdiff0(1:Nwalker)/temp)

!     ENDDO
!    ENDDO

!Energy = minval(Vexp(1:Nwalker,1:Natoms))

!write(*,*)"min of Vexp = ",Energy
!write(*,*)"(minimum of Vexp = 0d0)/True/False = ",(Energy.eq.0d0)

!IF ((Energy.eq.0d0).and.(temp.lt.1d0)) THEN

!temp = temp * 5d0
!write(*,*)"NEW !!! temperature T = ",temp

!GOTO 20
!ENDIF

IF (myrank.eq.0) THEN
      write(*,*)"INITIALIZATION!!!!!!!!!"
      write(*,*)"V_min = ",V_min
      write(*,*)"Temperature of simulation T = ",temp

      write(*,*)"Vdiff_0(1) = "
 !   DO iw=1,Nwalker
      write(*,*)Vdiff0(1)+V_min
 !   ENDDO

      write(*,*)"V_exp(1) = "
!    DO iw=1,Nwalker
      write(*,*)Vexp(1)
!    ENDDO

ENDIF 


!Vexp(1:Nwalker,1:Natoms) = dexp(-1d0 * Vdiff0X(1:Nwalker,1:Natoms)/temp)

!IF (Energy.eq.0d0) THEN

!write(*,*)"programm will stop, because of overflow problem"

!stop
!ENDIF
     !! Akashi
      weight(:)=1d0
      !call calc_moment1(Vexp(1:Nwalker),Nwalker,weight,mom1(1:Nwalker),Vexp_ave)
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
ENDIF

!     Vdiff0X(:,:)=0d0
!     Vdiff0Y(:,:)=0d0
!     Vdiff0Z(:,:)=0d0
     !!/END of initial normalization of p
     !/END of initialization


     ! START of ITERATIONS     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      B_coeff=sqrt(2d0*temp*inv_G)

      Voff = .false.
      call system_clock(ct0, count_rate, count_max)
      DO istep =1, Nstep
       !debug
       !call system_clock(ct0, count_rate, count_max)
       !IF(istep.gt.350)stop !! debug
       weight(:)=1d0

       IF(Vswitch .and. (istep .le. 30) .and. (i_restart .gt. 1)) THEN
        Voff = .true.
        Vswitch = .false.
       ENDIF

       IF(Voff .and. istep.gt.30)THEN !! first 30 steps V is off
        Voff = .false.
        Vswitch=.true.
       ENDIF

      
       IF(Vswitch) THEN
      ! Calc log constant
        call calc_Udiff(Nwalker,Natoms,x_pos(1:Nwalker,1:Natoms),y_pos(1:Nwalker,1:Natoms),z_pos(1:Nwalker,1:Natoms),&
&                       UdiffX(1:Nwalker,1:Natoms),UdiffY(1:Nwalker,1:Natoms),UdiffZ(1:Nwalker,1:Natoms))

        call calc_Udiff(Nwalker,Natoms,x_pos(1:Nwalker,1:Natoms),y_pos(1:Nwalker,1:Natoms),z_pos(1:Nwalker,1:Natoms),&
&                       VdiffX(1:Nwalker,1:Natoms),VdiffY(1:Nwalker,1:Natoms),VdiffZ(1:Nwalker,1:Natoms),                                                  &
&                       Vdiff2X(1:Nwalker,1:Natoms),Vdiff2Y(1:Nwalker,1:Natoms),Vdiff2Z(1:Nwalker,1:Natoms) )

        VdiffX(:,:) = (1d0 - ratio)*VdiffX(:,:)
        VdiffY(:,:) = (1d0 - ratio)*VdiffY(:,:)
        VdiffZ(:,:) = (1d0 - ratio)*VdiffZ(:,:)

        Vdiff2X(:,:) = (1d0 - ratio)*Vdiff2X(:,:)
        Vdiff2Y(:,:) = (1d0 - ratio)*Vdiff2Y(:,:)
        Vdiff2Z(:,:) = (1d0 - ratio)*Vdiff2Z(:,:)

        ! debug
        IF (Debag) THEN
        write(*,*)"UdiffX"
        write(*,*)UdiffX(1,:)
        write(*,*)"UdiffY"
        write(*,*)UdiffY(1,:)
        write(*,*)"UdiffZ"
        write(*,*)UdiffZ(1,:)
   
        write(*,*)"VdiffX"
        write(*,*)VdiffX(1,:)
        write(*,*)"VdiffY"
        write(*,*)VdiffY(1,:)
        write(*,*)"VdiffZ"
        write(*,*)VdiffZ(1,:)
   
        write(*,*)"Vdiff2X"
        write(*,*)Vdiff2X(1,:)
        write(*,*)"Vdiff2Y"
        write(*,*)Vdiff2Y(1,:)
        write(*,*)"VdiffZ"
        write(*,*)Vdiff2Z(1,:)
        ENDIF
        !stop
        !/debug

        F_func(1:Nwalker) =0d0
        DO ina=1, Natoms 
         F_func(1:Nwalker) = F_func(1:Nwalker) + &
&                           Vdiff2X(1:Nwalker,ina)+1d0/temp*VdiffX(1:Nwalker,ina)*(VdiffX(1:Nwalker,ina)-UdiffX(1:Nwalker,ina)) + &
&                           Vdiff2Y(1:Nwalker,ina)+1d0/temp*VdiffY(1:Nwalker,ina)*(VdiffY(1:Nwalker,ina)-UdiffY(1:Nwalker,ina)) + &
&                           Vdiff2Z(1:Nwalker,ina)+1d0/temp*VdiffZ(1:Nwalker,ina)*(VdiffZ(1:Nwalker,ina)-UdiffZ(1:Nwalker,ina)) 
        ENDDO

        lc = lc + inv_G*dt*F_ave

        !/Calc log constant

        ! birth-death process 1
        !
        !
        !! calculate factor P

        call calc_Udiff(Nwalker,Natoms,x_pos(1:Nwalker,1:Natoms),y_pos(1:Nwalker,1:Natoms),z_pos(1:Nwalker,1:Natoms),&
&                    UdiffX(1:Nwalker,1:Natoms),UdiffY(1:Nwalker,1:Natoms),UdiffZ(1:Nwalker,1:Natoms))
        call calc_Udiff(Nwalker,Natoms,x_pos(1:Nwalker,1:Natoms),y_pos(1:Nwalker,1:Natoms),z_pos(1:Nwalker,1:Natoms),&
&                    VdiffX(1:Nwalker,1:Natoms),VdiffY(1:Nwalker,1:Natoms),VdiffZ(1:Nwalker,1:Natoms),                                                     &
&                    Vdiff2X(1:Nwalker,1:Natoms),Vdiff2Y(1:Nwalker,1:Natoms),Vdiff2Z(1:Nwalker,1:Natoms) )


        VdiffX(:,:) = (1d0 - ratio)*VdiffX(:,:)
        VdiffY(:,:) = (1d0 - ratio)*VdiffY(:,:)
        VdiffZ(:,:) = (1d0 - ratio)*VdiffZ(:,:)

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
        call calc_moment1(F_func(1:Nwalker),Nwalker,weight,mom1(1:Nwalker),F_ave_pr)
        call MPI_ALLREDUCE(F_ave_pr, F_ave, 1, MPI_REAL8, & 
   &                     MPI_SUM, MPI_COMM_WORLD,ierr)
        F_ave = F_ave/dble(nprocs)
        mom1(1:Nwalker) = mom1(1:Nwalker) + F_ave_pr - F_ave
        prob(1:Nwalker)=exp(0.5d0*dt*inv_G*mom1(1:Nwalker))

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
   &      pr_accum((npr+1):(npr+Nwalker)) + pr_accum(npr)
        ENDDO

        pr_accum(:) = pr_accum(:)/pr_accum(Nwalker*nprocs)

        !IF(myrank.eq.0)THEN
        ! DO iw = 1, Nwalker*nprocs
        !  write(*,*)"iw, pr_accum", iw, pr_accum(iw), myrank
        ! ENDDO
        !ENDIF
        !call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        !IF(myrank.eq.1)THEN
        ! DO iw = 1, Nwalker*nprocs
        !  write(*,*)"iw, pr_accum", iw, pr_accum(iw), myrank
        ! ENDDO
        !ENDIF

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
           myincr = iw
           exit
          ENDIF
         ENDDO
        ENDIF

!        write(*,*)"proc,myincr",myrank, myincr
        iw   =  1
        incr = myincr + 1

        DO while (iw .le. Nwalker)
         IF( ((dble(iw-1)+ myrank*Nwalker + rand_p )/dble(Nwalker*nprocs)).lt.pr_accum(incr) )THEN
          id_walker(iw) = incr
          iw = iw + 1
         ELSE
          incr = incr + 1
         ENDIF

         !IF(incr .gt. (myrank+1)*Nwalker )THEN
         ! write(*,*)"error: incremental"
         ! stop
         !ENDIF
        ENDDO
!        write(*,*)"proc, incr", myrank, incr

        !! allreduce can be used like allgather

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
      !/birth-death process 1

       !call system_clock(ct1)
       !write(*,*)"sec/bd1=",dble(ct1-ct0)/dble(count_rate)
      ENDIF ! Vswitch
      !!!!!!!!!!!!!!!!!!! Langevin step 
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! generate random noise
      call gen_noise(Nwalker,Natoms,dWX(1:Nwalker,1:Natoms),dWY(1:Nwalker,1:Natoms),dWZ(1:Nwalker,1:Natoms))

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!/generate random noise
      !

      ! recalculation of distances between atoms due to changes in walkers
      !call calc_rad_m2(Nwalker,Natoms,x_pos(1:Nwalker,1:Natoms),y_pos(1:Nwalker,1:Natoms),z_pos(1:Nwalker,1:Natoms),rad_m2(1:Nwalker,1:Natoms,1:Natoms))

!     print *,istep

!     DO iw=1,Nwalker
!     DO ix=1,Natoms
!     DO iy=1,Natoms
!     print *,rad_m2(1,ix,iy)," 1  ",ix,iy
!     ENDDO
!     ENDDO
!     ENDDO

      !!
      !!!!!!!!!! calculate coefficients
      call calc_Udiff(Nwalker,Natoms,x_pos(1:Nwalker,1:Natoms),y_pos(1:Nwalker,1:Natoms),z_pos(1:Nwalker,1:Natoms),&
&                     UdiffX(1:Nwalker,1:Natoms),UdiffY(1:Nwalker,1:Natoms),UdiffZ(1:Nwalker,1:Natoms))
      !! Akashi, calculate the minimum absolute force
      sq_force(1:Nwalker)=0d0
      DO ina=1, Natoms
       sq_force(1:Nwalker) = sq_force(1:Nwalker) + UdiffX(1:Nwalker, ina)*UdiffX(1:Nwalker, ina) &
      &                                          + UdiffY(1:Nwalker, ina)*UdiffY(1:Nwalker, ina) &
      &                                          + UdiffZ(1:Nwalker, ina)*UdiffZ(1:Nwalker, ina)
      ENDDO
      min_abs_force_pr = sqrt(minval(sq_force(1:Nwalker)))
      call MPI_ALLREDUCE(min_abs_force_pr, min_abs_force, 1, MPI_REAL8, & 
   &                     MPI_MIN, MPI_COMM_WORLD,ierr)
      !!/Akashi, calculate the minimum absolute force

      IF(Vswitch) THEN
       call calc_Udiff(Nwalker,Natoms,x_pos(1:Nwalker,1:Natoms),y_pos(1:Nwalker,1:Natoms),z_pos(1:Nwalker,1:Natoms),&
&                      VdiffX(1:Nwalker,1:Natoms),VdiffY(1:Nwalker,1:Natoms),VdiffZ(1:Nwalker,1:Natoms))
     
       VdiffX(:,:) = (1d0 - ratio)*VdiffX(:,:)
       VdiffY(:,:) = (1d0 - ratio)*VdiffY(:,:)
       VdiffZ(:,:) = (1d0 - ratio)*VdiffZ(:,:)
 
      ENDIF

      A_coeffX(1:Nwalker,1:Natoms) = inv_G*(-UdiffX(1:Nwalker,1:Natoms)+2d0*VdiffX(1:Nwalker,1:Natoms))
      A_coeffY(1:Nwalker,1:Natoms) = inv_G*(-UdiffY(1:Nwalker,1:Natoms)+2d0*VdiffY(1:Nwalker,1:Natoms))
      A_coeffZ(1:Nwalker,1:Natoms) = inv_G*(-UdiffZ(1:Nwalker,1:Natoms)+2d0*VdiffZ(1:Nwalker,1:Natoms))

!!!!!! Debaggin code 
  
      IF ((Debag).and.((mod(istep,steptowrite).eq.0).or.(istep.eq.1))) THEN
       write(11,'(A8)')"TIMESTEP"
       write(11,'(i0)')istep
       write(11,'(A21)')"Data for first walker"

       !write(11,'(A21)')"radius for first atom"

       !DO ina=2,Natoms
       ! write(11,'(1i4,2e18.8e4)')ina,rad_m2(1,1,ina)
       !ENDDO

       write(11,'(A55)')"Comparision between U' and V' potentials for each atom:"

       DO ina=1,Natoms
        write(11,'(1i4,1A3,2e18.8e4)')ina," X ",UdiffX(1,ina),VdiffX(1,ina)
        write(11,'(1i4,1A3,2e18.8e4)')ina," Y ",UdiffY(1,ina),VdiffY(1,ina)
        write(11,'(1i4,1A3,2e18.8e4)')ina," Z ",UdiffZ(1,ina),VdiffZ(1,ina)
       ENDDO

       write(11,'(A47)')"Check the values of V potentials for each atom:"

       !DO ina=1,Natoms
       write(11,'(1i4,1A3,1e18.8e4)')" X ",Vdiff0(1) + V_ini      
       !ENDDO  



       write(11,'(A57)')"Check the values of A*dt, dt and B*sqrt(dt)*dW for each atom:"

       DO ina=1,Natoms
        write(11,'(1i4,1A3,3e18.8e4)')ina," X ",A_coeffX(1,ina)*dt,dt,B_coeff*sqrt(dt)*dWX(1,ina)
        write(11,'(1i4,1A3,3e18.8e4)')ina," Y ",A_coeffY(1,ina)*dt,dt,B_coeff*sqrt(dt)*dWY(1,ina)
        write(11,'(1i4,1A3,3e18.8e4)')ina," Z ",A_coeffZ(1,ina)*dt,dt,B_coeff*sqrt(dt)*dWZ(1,ina)
       ENDDO


       
      

      ENDIF
!!!!!!!!!!!!!!!!!!!


! CHECK FOR TIME STEP AND if it is necessesary then to make correction !!!!!!!!!!! 
! This is a adaptive time step in the range [0:maxtimestep]

! IF ((istep.le.10).or.(mod(istep,100).eq.0)) THEN 

! 10 continue

! displacement(1:Natoms) = (abs(A_coeffX(1,1:Natoms)) + abs(A_coeffY(1,1:Natoms)) + abs(A_coeffZ(1,1:Natoms))) * dt + B_coeff*sqrt(dt) * (abs(dWX(1,1:Natoms)) + abs(dWY(1,1:Natoms)) + abs(dWZ(1,1:Natoms)))


! IF (maxval(displacement(1:Natoms)).gt.0.1d0) THEN  

! dt = dt / 2
! write(*,*)"step number is ",istep,"   and new time step is ",dt

! GOTO 10 
! ENDIF 

! IF ((maxval(displacement(1:Natoms)).lt.5.0d-2).and.(dt.lt.maxtimestep).and.(.NOT.Vswitch)) THEN

! dt = dt * 2
! write(*,*)"step number is ",istep,"   and new time step is ",dt 

! GOTO 10
! ENDIF

! ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      
      !!/calculate coefficients
      !
      !! increment
      x_pos(1:Nwalker,1:Natoms) = x_pos(1:Nwalker,1:Natoms)+A_coeffX(1:Nwalker,1:Natoms)*dt+B_coeff*sqrt(dt)*dWX(1:Nwalker,1:Natoms)
      y_pos(1:Nwalker,1:Natoms) = y_pos(1:Nwalker,1:Natoms)+A_coeffY(1:Nwalker,1:Natoms)*dt+B_coeff*sqrt(dt)*dWY(1:Nwalker,1:Natoms)
      z_pos(1:Nwalker,1:Natoms) = z_pos(1:Nwalker,1:Natoms)+A_coeffZ(1:Nwalker,1:Natoms)*dt+B_coeff*sqrt(dt)*dWZ(1:Nwalker,1:Natoms)


      ! recalculate distances, because the atomic positions are changed
      !call calc_rad_m2(Nwalker,Natoms,x_pos(1:Nwalker,1:Natoms),y_pos(1:Nwalker,1:Natoms),z_pos(1:Nwalker,1:Natoms),rad_m2(1:Nwalker,1:Natoms,1:Natoms))

      !!/increment
      !
      !
      !/Langevin step

       !call system_clock(ct2)
       !write(*,*)"sec/Lange=",dble(ct2-ct1)/dble(count_rate)

      ! birth-death process 2
      IF(Vswitch) THEN
      !

      !! calculate factor P
       call calc_Udiff(Nwalker,Natoms,x_pos(1:Nwalker,1:Natoms),y_pos(1:Nwalker,1:Natoms),z_pos(1:Nwalker,1:Natoms),&
&                      UdiffX(1:Nwalker,1:Natoms),UdiffY(1:Nwalker,1:Natoms),UdiffZ(1:Nwalker,1:Natoms))
       call calc_Udiff(Nwalker,Natoms,x_pos(1:Nwalker,1:Natoms),y_pos(1:Nwalker,1:Natoms),z_pos(1:Nwalker,1:Natoms),&
&                      VdiffX(1:Nwalker,1:Natoms),VdiffY(1:Nwalker,1:Natoms),VdiffZ(1:Nwalker,1:Natoms),                                                    &
&                      Vdiff2X(1:Nwalker,1:Natoms),Vdiff2Y(1:Nwalker,1:Natoms),Vdiff2Z(1:Nwalker,1:Natoms) )

       VdiffX(:,:) = (1d0 - ratio)*VdiffX(:,:)
       VdiffY(:,:) = (1d0 - ratio)*VdiffY(:,:)
       VdiffZ(:,:) = (1d0 - ratio)*VdiffZ(:,:)

       Vdiff2X(:,:) = (1d0 - ratio)*Vdiff2X(:,:)
       Vdiff2Y(:,:) = (1d0 - ratio)*Vdiff2Y(:,:)
       Vdiff2Z(:,:) = (1d0 - ratio)*Vdiff2Z(:,:)

       F_func(1:Nwalker) = 0d0
       DO ina=1, Natoms
        F_func(1:Nwalker) =  F_func(1:Nwalker) + &
&                             Vdiff2X(1:Nwalker,ina)+1d0/temp*VdiffX(1:Nwalker,ina)*(VdiffX(1:Nwalker,ina)-UdiffX(1:Nwalker,ina)) + & 
&                             Vdiff2Y(1:Nwalker,ina)+1d0/temp*VdiffY(1:Nwalker,ina)*(VdiffY(1:Nwalker,ina)-UdiffY(1:Nwalker,ina)) + &
&                             Vdiff2Z(1:Nwalker,ina)+1d0/temp*VdiffZ(1:Nwalker,ina)*(VdiffZ(1:Nwalker,ina)-UdiffZ(1:Nwalker,ina))
       ENDDO

       call calc_moment1(F_func(1:Nwalker),Nwalker,weight,mom1(1:Nwalker),F_ave_pr)
       call MPI_ALLREDUCE(F_ave_pr, F_ave, 1, MPI_REAL8, & 
   &                    MPI_SUM, MPI_COMM_WORLD,ierr)
       F_ave = F_ave/dble(nprocs)
       mom1(1:Nwalker) = mom1(1:Nwalker) + F_ave_pr - F_ave
       prob(1:Nwalker)=exp(0.5d0*dt*inv_G*mom1(1:Nwalker))
        !IF(myrank.eq.0)THEN
        ! DO iw = 1, Nwalker
        !  write(*,'(A10, i8,10f14.6)')"iw, x_pos", iw, x_pos(iw,:)
        ! ENDDO
        ! DO iw = 1, Nwalker
        !  write(*,*)"iw, F_func,weight", iw, F_func(iw),weight(iw), myrank
        ! ENDDO
        ! DO iw = 1, Nwalker
        !  write(*,*)"iw, prob", iw, prob(iw), myrank
        ! ENDDO
        !ENDIF

       !!/calculate factor P
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
       ! IF(myrank.eq.0)THEN
       !  DO iw = 1, Nwalker*nprocs
       !   write(*,*)"iw, pr_accum_pr", iw, pr_accum_pr(iw), myrank
       !  ENDDO
       ! ENDIF

       call MPI_ALLREDUCE(pr_accum_pr,pr_accum, Nwalker*nprocs, MPI_REAL8, & 
   &                     MPI_SUM, MPI_COMM_WORLD,ierr)
       DO iproc = 2, nprocs
        npr = (iproc-1)*Nwalker
        pr_accum((npr+1):(npr+Nwalker)) =                   &
   &      pr_accum((npr+1):(npr+Nwalker)) + pr_accum(npr)
       ENDDO
       ! IF(myrank.eq.0)THEN
       !  DO iw = 1, Nwalker*nprocs
       !   write(*,*)"iw, pr_accum", iw, pr_accum(iw), myrank
       !  ENDDO
       ! ENDIF
       ! call MPI_BARRIER(MPI_COMM_WORLD, ierr)
       ! IF(myrank.eq.1)THEN
       !  DO iw = 1, Nwalker*nprocs
       !   write(*,*)"iw, pr_accum", iw, pr_accum(iw), myrank
       !  ENDDO
       ! ENDIF

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
          myincr = iw
          exit
         ENDIF
        ENDDO
       ENDIF

!       write(*,*)"proc,myincr",myrank, myincr
       iw   =  1
       incr = myincr + 1

       DO while (iw .le. Nwalker)
        IF( ((dble(iw-1)+ myrank*Nwalker + rand_p )/dble(Nwalker*nprocs)).lt.pr_accum(incr) )THEN
         id_walker(iw) = incr
         iw = iw + 1
        ELSE
         incr = incr + 1
        ENDIF

        !IF(incr .gt. (myrank+1)*Nwalker )THEN
        ! write(*,*)"error: incremental"
        ! stop
        !ENDIF
       ENDDO
!       write(*,*)"proc, incr", myrank, incr

       !! allreduce can be used like allgather

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
       !
       !
       !call system_clock(ct3)
       !write(*,*)"sec/bd2=",dble(ct3-ct2)/dble(count_rate)
       ENDIF ! Vswitch
       !/birth-death process 2

       ! recalculation of distances between atoms due to changes in walkers
       !call calc_rad_m2(Nwalker,Natoms,x_pos(1:Nwalker,1:Natoms),y_pos(1:Nwalker,1:Natoms),z_pos(1:Nwalker,1:Natoms),rad_m2(1:Nwalker,1:Natoms,1:Natoms))


!       IF((mod(istep,steptowrite).eq.0).or.(istep.eq.1))THEN
!        write(*,*)istep
!        write(*,'(3f16.8)')maxval(x_pos(1:Nwalker,1:Natoms)),maxval(y_pos(1:Nwalker,1:Natoms)),maxval(z_pos(1:Nwalker,1:Natoms))
!        write(*,'(3f16.8)')minval(x_pos(1:Nwalker,1:Natoms)),minval(y_pos(1:Nwalker,1:Natoms)),minval(z_pos(1:Nwalker,1:Natoms))
!        write(*,*)"Vswitch,Voff",Vswitch,Voff
!       ENDIF

       IF(myrank.eq.0)THEN
!        write(21,'(2i10,1f16.8)')istep,Nwalker,lc
       ENDIF
     
                                                    ! calc distribution function
       IF((mod(istep,steptowrite).eq.0).or.(istep.eq.1))THEN

!       write(13,'(i10,f16.8)')istep,maxval(x_pos(1:Nwalker,1:Natoms))


       !! Akashi modified the shape of weight
       !weight(1:Nwalker,1:Natoms)=1d0
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

         call calc_Udiff(Nwalker,Natoms,x_pos(1:Nwalker,1:Natoms),y_pos(1:Nwalker,1:Natoms),z_pos(1:Nwalker,1:Natoms),&
&                        Vdiff0(1:Nwalker))

         Vdiff0(:) = (1d0 - ratio)*Vdiff0(:)


         DO iw=1, Nwalker
          Vdiff0(iw)=Vdiff0(iw) - V_ini
         ENDDO 
        ENDIF


        V_min_pr = Vdiff0(1)
        DO iw = 1, Nwalker
         IF (V_min_pr.gt.Vdiff0(iw)) THEN 
          V_min_pr = Vdiff0(iw)
         ENDIF
        ENDDO
        call MPI_ALLREDUCE(V_min_pr, V_min, 1, MPI_REAL8, & 
   &                    MPI_MIN, MPI_COMM_WORLD,ierr)
       
        lnweight_min=lc-V_min/temp
       !! Akashi modified the shape of weight
       !weight(1:Nwalker,1:Natoms)=exp((V_min-Vdiff0X(1:Nwalker,1:Natoms))/temp)

        DO iw = 1, Nwalker
         weight(iw)=exp((V_min-Vdiff0(iw))/temp)
        ENDDO

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
        ENDIF
!       write(*,*)"debug: xyz_ave_p(1)=",x_ave_p(1),y_ave_p(1), z_ave_p(1)
!       write(*,*)"debug: V_min=",V_min
!       write(*,*)"debug: V_min-Vdiff0X(:)=",V_min-Vdiff0X(1:Nwalker,1)
!       write(*,*)"debug: weight(:)=",weight(:)

       !!call calc_dist(x_pos(1:Nwalker,1:Natoms),y_pos(1:Nwalker,1:Natoms),z_pos(1:Nwalker,1:Natoms),weight(1:Nwalker,1:Natoms),Nwalker,Natoms,stepXdist,stepYdist,stepZdist,dist_P,Nx,Ny,Nz,lnweight_min)
       !!/calculate potential correction and exp(-bV)


!       DO ina = 1, Natoms
!       max_d_corr(:,ina) = maxloc(dist_P(:,:,:,ina))    !!!!   The Atomic positions in the special file  !!!!!!!!!!!!!!
!       max_d(:,ina)      = maxloc(dist_Q(:,:,:,ina))         !!!!  We save reaction path as the positions of the each atom with most probability  !!!!!!!!!!!!!!
!       ENDDO

!       DO ina = 1, Natoms
!         coorX_max(ina) = stepXdist(max_d_corr(1,ina),ina)
!         coorY_max(ina) = stepYdist(max_d_corr(2,ina),ina)
!         coorZ_max(ina) = stepZdist(max_d_corr(3,ina),ina)
!       ENDDO

        IF(myrank.eq.0)THEN
         call calc_Energy(Natoms,x_ave_p(1:Natoms),y_ave_p(1:Natoms),z_ave_p(1:Natoms),Energy_p,E_atoms(1:Natoms))
         call calc_Energy(Natoms,x_ave_q(1:Natoms),y_ave_q(1:Natoms),z_ave_q(1:Natoms),Energy_q,E_atoms(1:Natoms))
! call calc_Energy(Natoms,coorX_max(1:Natoms),coorY_max(1:Natoms),coorZ_max(1:Natoms),Energy)
        ENDIF
        call MPI_BCAST(Energy_q, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(Energy_p, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
! TO check the delta parameter and to change if it is nessesary
        IF ((i_restart.lt.N_restart).and.(istep.gt.400)) THEN
!IF (myrank.eq.0) THEN
!         write(*,*)"E_start = ",E_start
!         write(*,*)"Energy = ",Energy_q
!ENDIF
        ! stop !! AKASHI debug
         IF (E_start .gt. Energy_q) THEN 
!IF (myrank.eq.0) THEN 
!          write(*,*)"GO TO 200 label to change the delta parameter, because the energy is decreased! "
!ENDIF
          Nstep = Nstep*N_restart
          GOTO 200
         ENDIF

        ENDIF
!  \check the delta parameter 



! Write the reaction path with Q and P functions + Energy of system  
      !ina = 1
      !write(17,'(1i6,3e20.8e4)')istep,dist_Q(max_d(1,ina),max_d(2,ina),max_d(3,ina),ina),dist_P(max_d_corr(1,ina),max_d_corr(2,ina),max_d_corr(3,ina),ina),Energy
       IF(myrank.eq.0.and.(i_restart.lt.N_restart)) THEN 
        write(17,'(1i6,3e20.8e4,i10)')istep+restart_step,Energy_q, Energy_p,min_abs_force,Nwalker*nprocs

        write(18,'(1i6,1e20.8e4)', advance='no')istep+restart_step,Energy_q

       DO ina = 1, Natoms-1
       write(18,'(1e20.8e4)', advance='no')E_atoms(ina)
       ENDDO
       write(18,'(1e20.8e4)')E_atoms(Natoms)

       ENDIF


! To write Energy and coordinates during relaxation

       IF(i_restart.eq.N_restart) THEN

         x_ave_q(1:Natoms) = x_pos(1,1:Natoms) 
         y_ave_q(1:Natoms) = y_pos(1,1:Natoms)
         z_ave_q(1:Natoms) = z_pos(1,1:Natoms)

         call calc_Energy(Natoms,x_ave_q(1:Natoms),y_ave_q(1:Natoms),z_ave_q(1:Natoms),Energy_q,E_atoms(1:Natoms))

        iw=200+myrank
        write(iw,'(1i6,1e20.8e4)', advance='no')istep+restart_step,Energy_q

       DO ina = 1, Natoms-1
       write(iw,'(1e20.8e4)', advance='no')E_atoms(ina)
       ENDDO
       write(iw,'(1e20.8e4)')E_atoms(Natoms)

        iw=100+myrank
         call write_positions(Natoms,x_min,x_max,y_min,y_max,z_min,z_max,iw,istep+restart_step,atomictype,x_ave_q(1:Natoms),y_ave_q(1:Natoms),z_ave_q(1:Natoms))
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
       ENDIF



! write atomic coordinates for P distribution in the LAMMPS format !!!!!!!



! write atomic coordinates for P distribution in the LAMMPS format !!!!!!!

!        x_min_pr = minval(x_pos(:,:))
!        x_max_pr = maxval(x_pos(:,:))

!        y_min_pr = minval(y_pos(:,:))
!        y_max_pr = maxval(y_pos(:,:))

!        z_min_pr = minval(z_pos(:,:))
!        z_max_pr = maxval(z_pos(:,:))

!        call MPI_ALLREDUCE(x_min_pr, x_min, 1, MPI_REAL8, & 
!   &                    MPI_SUM, MPI_COMM_WORLD,ierr)
!        call MPI_ALLREDUCE(x_max_pr, x_max, 1, MPI_REAL8, & 
!   &                    MPI_SUM, MPI_COMM_WORLD,ierr)
!        call MPI_ALLREDUCE(y_min_pr, y_min, 1, MPI_REAL8, & 
!   &                    MPI_SUM, MPI_COMM_WORLD,ierr)
!        call MPI_ALLREDUCE(y_max_pr, y_max, 1, MPI_REAL8, & 
!   &                    MPI_SUM, MPI_COMM_WORLD,ierr)
!        call MPI_ALLREDUCE(z_min_pr, z_min, 1, MPI_REAL8, & 
!   &                    MPI_SUM, MPI_COMM_WORLD,ierr)
!        call MPI_ALLREDUCE(z_max_pr, z_max, 1, MPI_REAL8, & 
!   &                    MPI_SUM, MPI_COMM_WORLD,ierr)
        !call write_positions(Natoms,x_min,x_max,y_min,y_max,z_min,z_max,12,istep,atomictype,coorX_max(1:Natoms),coorY_max(1:Natoms),coorZ_max(1:Natoms))

        IF(myrank.eq.0)THEN
         call write_positions(Natoms,x_min,x_max,y_min,y_max,z_min,z_max,12,istep+restart_step,atomictype,x_ave_q(1:Natoms),y_ave_q(1:Natoms),z_ave_q(1:Natoms))
!         call write_positions(Natoms,x_min,x_max,y_min,y_max,z_min,z_max,15,istep+restart_step,atomictype,x_ave_p(1:Natoms),y_ave_p(1:Natoms),z_ave_p(1:Natoms))






! write atomic coordinates ONLY for FIST  WALKER !!! in the LAMMPS format
! !!!!!!!

!         call write_positions(Natoms,x_min,x_max,y_min,y_max,z_min,z_max,13,istep+restart_step,atomictype,x_pos(1,1:Natoms),y_pos(1,1:Natoms),z_pos(1,1:Natoms))
        ENDIF


       ENDIF !steptowrite                        !/calc distribution function


       !call system_clock(ct4)
       !write(*,*)"sec/loop=",dble(ct4-ct0)/dble(count_rate)
       !stop
       IF(istep.eq.100)THEN
        call system_clock(ct4)
        write(*,*)"sec/100loop=",dble(ct4-ct0)/dble(count_rate)
        stop
       ENDIF
      ENDDO              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   END OF ITERATION     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      restart_step = restart_step + Nstep


     ENDDO              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   END OF RESTART LOOP  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


     close(11)
     close(12)
!     close(13)
!     close(15)
     close(21)
     close(17)
     close(18)
     iw = 51 + myrank
     close(iw)

IF (myrank.eq.0) THEN
     write(*,'(A8,3e20.8e4)')"max pos.",maxval(x_pos(1:Nwalker,1:Natoms)),maxval(y_pos(1:Nwalker,1:Natoms)),maxval(z_pos(1:Nwalker,1:Natoms))
     write(*,'(A8,3e20.8e4)')"min pos.",minval(x_pos(1:Nwalker,1:Natoms)),minval(y_pos(1:Nwalker,1:Natoms)),minval(z_pos(1:Nwalker,1:Natoms))
ENDIF

   end subroutine solve
  !-----------


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
   subroutine calc_Energy(Natoms,x_pos,y_pos,z_pos,Energy,E_atoms)
   !------------------------
   !
   !This routine calculates the potential Energy for ALL atoms or for atomic system
   !
   implicit none
   
   integer, intent(in) :: Natoms 

   real(8) :: x_pos(Natoms)
   real(8) :: y_pos(Natoms)
   real(8) :: z_pos(Natoms)

   real(8), intent(out) :: Energy
   real(8), intent(out) :: E_atoms(1:Natoms)

   real(8) :: rad_m2(Natoms-1,Natoms)

   real(8) :: rad_ma(Natoms-1,Natoms)
   real(8) :: rad_mb(Natoms-1,Natoms)

   real(8),parameter  :: sigma_6=1.544804416d3    ! Angstrom^6
 
   real(8),parameter  :: sigma_12  =    2.386420683693101056d6     ! Angstrom^12

   real(8),parameter  :: epsilon4=4.0104d1 / 6d0

   integer :: ina_i,ina_j 



   Energy = 0d0
   E_atoms(1:Natoms) = 0d0

   IF (Natoms.lt.2) THEN 
                    stop
   ENDIF 

   rad_m2(:,:) = 0d0

      DO ina_i=1, Natoms-1
       DO ina_j=ina_i+1, Natoms

            rad_m2(ina_i,ina_j) = 1d0 / ((x_pos(ina_i)-x_pos(ina_j))**2d0 + (y_pos(ina_i)-y_pos(ina_j))**2d0 + (z_pos(ina_i)-z_pos(ina_j))**2d0)

       ENDDO
      ENDDO
   

    rad_ma(:,:) = rad_m2(:,:) ** 6   ! r^-12
    rad_mb(:,:) = rad_m2(:,:) ** 3   ! r^-6

      DO ina_i=1, Natoms
       DO ina_j=1, Natoms

       IF (ina_j.gt.ina_i) THEN

       E_atoms(ina_i) = E_atoms(ina_i) + epsilon4 * (sigma_12 * rad_ma(ina_i,ina_j)  - sigma_6 * rad_mb(ina_i,ina_j))

       ENDIF

       IF (ina_j.lt.ina_i) THEN      ! we use symmetry of distance matrix rad_ij = rad_ji

       E_atoms(ina_i) = E_atoms(ina_i) + epsilon4 * (sigma_12 * rad_ma(ina_j,ina_i)  - sigma_6 * rad_mb(ina_j,ina_i))

       ENDIF

       ENDDO
      ENDDO

    E_atoms(1:Natoms) = E_atoms(1:Natoms) / 2d0
    Energy = SUM(E_atoms(1:Natoms)) / Natoms

   end subroutine calc_Energy


   !------------------------
   !
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
   ! !! updated by Akashi, 2018/09/04
   !
   implicit none
   integer, intent(in) :: Nwalker
   real(8), intent(in) :: func(Nwalker)
   real(8), intent(in) :: weight(Nwalker)  !! added by Akashi
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
   ! !! updated by Akashi, 2018/06/30
   implicit none
   integer, intent(in) :: Nwalker,Natoms
   real(8), intent(in) :: func(Nwalker,Natoms)
   real(8), intent(in) :: weight(Nwalker)  !! added by Akashi
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
   subroutine write_positions(Natoms,x_min,x_max,y_min,y_max,z_min,z_max,n_file,istep,atomictype,coorX_max,coorY_max,coorZ_max)
   !------------------------
   !This subroutine write to the file the coordinates of atomic positions in the
   !format of
   !LAMMPS or OVITO  
   ! 
   implicit none
   integer, intent(in) :: Natoms
   integer, intent(in) :: atomictype,n_file,istep

   real(8), intent(in) :: x_min,x_max,y_min,y_max,z_min,z_max

   real(8), intent(in) :: coorX_max(Natoms)
   real(8), intent(in) :: coorY_max(Natoms)
   real(8), intent(in) :: coorZ_max(Natoms)

   real(8) :: coorX(Natoms)
   real(8) :: coorY(Natoms)
   real(8) :: coorZ(Natoms)


!   real(8), intent(out) :: rad_minus2(Nwalker,Natoms,Natoms)

   integer :: ina

! write atomic coordinates ONLY for FIST  WALKER !!! in the LAMMPS format
! !!!!!!!

    write(n_file,'(A14)')"ITEM: TIMESTEP"
    write(n_file,'(i0)')istep

    write(n_file,'(A21)')"ITEM: NUMBER OF ATOMS"

    write(n_file,'(i0)')Natoms

    write(n_file,'(A25)')"ITEM: BOX BOUNDS ss ss ss"

    write(n_file,'(f0.4," ",f0.4)')x_min,x_max
    write(n_file,'(f0.4," ",f0.4)')y_min,y_max
    write(n_file,'(f0.4," ",f0.4)')z_min,z_max

!!!! Recalculation of coordinates - for box coordinates  !!!



     DO ina = 1, Natoms

      coorX(ina) = (coorX_max(ina) - x_min) / (x_max - x_min)

      coorY(ina) = (coorY_max(ina) - y_min) / (y_max - y_min)

      coorZ(ina) = (coorZ_max(ina) - z_min) / (z_max - z_min)

     ENDDO

    write(n_file,'(A28)')"ITEM: ATOMS id type xs ys zs"
    DO ina = 1, Natoms

      write(n_file,'(i0," ",i0," ",f0.6," ",f0.6," ",f0.6)')ina,atomictype,coorX(ina),coorY(ina),coorZ(ina)

    ENDDO

end subroutine write_positions

program main
  use routines
  implicit none

  call initmpi()
  call stdin()
  call gen_walker()
  call solve()

  !call write_data()
end program

