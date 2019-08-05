
      module potentials
       implicit none

       interface calc_Udiff
        module procedure calc_Udiff0,calc_Udiff0_atom,calc_Udiff1,calc_Udiff12
       end interface calc_Udiff

       !interface calc_Vdiff
       ! module procedure calc_Vdiff0, calc_Vdiff1, calc_Vdiff12
       !end interface calc_Vdiff

       contains
!!!!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX!!!

        subroutine calc_Udiff0(Nwalker,Natoms,Nregions,Ncycles,x_pos,y_pos, z_pos,  pe)
         use val_mpi_lmp, only : lammps, comm_lammps, ptr
         use bistable1d_val, only : x_pos_mod, y_pos_mod, z_pos_mod
         integer, intent(in) :: Nwalker
         integer, intent(in) :: Natoms
         integer, intent(in) :: Nregions
         integer, intent(in) :: Ncycles

         real(8), intent(in) :: x_pos(Nwalker,Natoms)
         real(8), intent(in) :: y_pos(Nwalker,Natoms)
         real(8), intent(in) :: z_pos(Nwalker,Natoms)

         real(8), intent(out) :: pe(Nwalker)

         ! The analitical form of potential

         real(8) :: rad_m2(Nwalker,Natoms-1,Natoms)
         real(8) :: rad_ma(Nwalker,Natoms-1,Natoms)
         real(8) :: rad_mb(Nwalker,Natoms-1,Natoms)
   
         real(8),parameter  :: sigma_6=         1.544804416d3    ! Angstrom^6
         real(8),parameter  :: sigma_6_2=2d0  * 1.544804416d3    ! Angstrom^6
         real(8),parameter  :: sigma_12  =      2.386420683693101056d6    ! Angstrom^12
         real(8),parameter  :: sigma_12_2=2d0 * 2.386420683693101056d6    ! Angstrom^12
         real(8),parameter  :: sigma_12_7=7d0 * 2.386420683693101056d6    ! Angstrom^12
         real(8),parameter  :: epsilon24=4.0104d1   
         real(8),parameter  :: epsilon96=4.0104d1 * 4d0
         real(8),parameter  :: epsilon4= 4.0104d1 / 6d0
         integer :: iw,ina_i,ina_j


         pe(:) = 0d0


         call calc_rad_m2(Nwalker,Natoms,x_pos(1:Nwalker,1:Natoms),&
     &                                   y_pos(1:Nwalker,1:Natoms),&
     &                                   z_pos(1:Nwalker,1:Natoms),&
     &                        rad_m2(1:Nwalker,1:Natoms-1,1:Natoms))

         rad_ma(:,:,:) = rad_m2(:,:,:) ** 6   ! r^-12
         rad_mb(:,:,:) = rad_m2(:,:,:) ** 3   ! r^-6


         DO ina_i=1, Natoms-1
          DO ina_j=ina_i+1, Natoms

           pe(:) = pe(:) + epsilon4 * (sigma_12 * rad_ma(:,ina_i,ina_j)  - sigma_6 * rad_mb(:,ina_i,ina_j))

          ENDDO
         ENDDO

        end subroutine calc_Udiff0

!!!!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX!!!

        subroutine calc_Udiff0_atom(Nwalker,Natoms,Nregions,Ncycles,&
     &              x_pos, y_pos, z_pos,  pe_atom)
         use val_mpi_lmp, only : lammps, comm_lammps, ptr
         use bistable1d_val, only : x_pos_mod, y_pos_mod, z_pos_mod
         integer, intent(in) :: Nwalker
         integer, intent(in) :: Natoms
         integer, intent(in) :: Nregions
         integer, intent(in) :: Ncycles

         real(8), intent(in) :: x_pos(Nwalker,Natoms)
         real(8), intent(in) :: y_pos(Nwalker,Natoms)
         real(8), intent(in) :: z_pos(Nwalker,Natoms)

         real(8), intent(out) :: pe_atom(Nwalker,Natoms)

         ! The analitical form of potential

         real(8) :: rad_m2(Nwalker,Natoms-1,Natoms)
         real(8) :: rad_ma(Nwalker,Natoms-1,Natoms)
         real(8) :: rad_mb(Nwalker,Natoms-1,Natoms)
         real(8) :: pij(Nwalker)

         real(8),parameter  :: sigma_6=         1.544804416d3    ! Angstrom^6
         real(8),parameter  :: sigma_6_2=2d0  * 1.544804416d3    ! Angstrom^6
         real(8),parameter  :: sigma_12  =      2.386420683693101056d6  ! Angstrom^12
         real(8),parameter  :: sigma_12_2=2d0 * 2.386420683693101056d6  ! Angstrom^12
         real(8),parameter  :: sigma_12_7=7d0 * 2.386420683693101056d6  ! Angstrom^12

         real(8),parameter  :: epsilon24=4.0104d1
         real(8),parameter  :: epsilon96=4.0104d1 * 4d0
         real(8),parameter  :: epsilon4= 4.0104d1 / 6d0
         integer :: iw,ina_i,ina_j

         pij(:) = 0.0d0
         pe_atom(:,:) = 0.0d0

         call calc_rad_m2(Nwalker,Natoms,x_pos(1:Nwalker,1:Natoms),&
     &                                   y_pos(1:Nwalker,1:Natoms),&
     &                                   z_pos(1:Nwalker,1:Natoms),&
     &                        rad_m2(1:Nwalker,1:Natoms-1,1:Natoms))

         rad_ma(:,:,:) = rad_m2(:,:,:) ** 6   ! r^-12
         rad_mb(:,:,:) = rad_m2(:,:,:) ** 3   ! r^-6


         DO ina_i=1, Natoms-1
          DO ina_j=ina_i+1, Natoms

           ! pe(:) = pe(:) + epsilon4 * (sigma_12 * rad_ma(:,ina_i,ina_j) - sigma_6 * rad_mb(:,ina_i,ina_j))
           
           pij(:) = epsilon4 * (sigma_12 * rad_ma(:,ina_i,ina_j) - sigma_6 * rad_mb(:,ina_i,ina_j))

           pe_atom(:,ina_i) = pe_atom(:,ina_i) + pij(:)   !  epsilon4 * (sigma_12 * rad_ma(:,ina_i,ina_j) - sigma_6 * rad_mb(:,ina_i,ina_j))

           pe_atom(:,ina_j) = pe_atom(:,ina_j) + pij(:)   !  epsilon4 * (sigma_12 * rad_ma(:,ina_i,ina_j) - sigma_6 * rad_mb(:,ina_i,ina_j))

          ENDDO
         ENDDO


!          pe_atom(icy,1:Natoms) = pe_tmp(1:Natoms)

        end subroutine calc_Udiff0_atom
!!!!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX!!!

        subroutine calc_Udiff1(Nwalker,Natoms,Nregions, Ncycles,&
     &                         x_pos,y_pos,z_pos,dx_pe,dy_pe, dz_pe)
         use val_mpi_lmp, only : lammps, comm_lammps, ptr
         use bistable1d_val, only : x_pos_mod, y_pos_mod, z_pos_mod
         integer, intent(in) :: Nwalker
         integer, intent(in) :: Natoms
         integer, intent(in) :: Nregions
         integer, intent(in) :: Ncycles
   
         real(8), intent(in) :: x_pos(Nwalker,Natoms)
         real(8), intent(in) :: y_pos(Nwalker,Natoms)
         real(8), intent(in) :: z_pos(Nwalker,Natoms)


         real(8), intent(out) :: dx_pe(Nwalker,Natoms)
         real(8), intent(out) :: dy_pe(Nwalker,Natoms)
         real(8), intent(out) :: dz_pe(Nwalker,Natoms)

         real(8) :: rad_m2(Nwalker,Natoms-1,Natoms)
         real(8) :: rad_ma(Nwalker,Natoms-1,Natoms)
         real(8) :: rad_mb(Nwalker,Natoms-1,Natoms)


         real(8) :: VdiffX(Nwalker,Natoms), VdiffY(Nwalker,Natoms),VdiffZ(Nwalker,Natoms)

!         real(8) :: x_pos_p(Nwalker,Natoms)
!         real(8) :: x_pos_m(Nwalker,Natoms)
!         real(8) :: y_pos_p(Nwalker,Natoms)
!         real(8) :: y_pos_m(Nwalker,Natoms)
!         real(8) :: z_pos_p(Nwalker,Natoms)
!         real(8) :: z_pos_m(Nwalker,Natoms)

!         real(8) :: pe_x_p(Nwalker)
!         real(8) :: pe_x_m(Nwalker)
!         real(8) :: pe_y_p(Nwalker)
!         real(8) :: pe_y_m(Nwalker)
!         real(8) :: pe_z_p(Nwalker)
!         real(8) :: pe_z_m(Nwalker)


         !real(8),parameter  :: ratio=0.02d0 ! T/T0
         !real(8),parameter  :: ratio=0.1d0 ! T/T0
         !real(8),parameter  :: ratio=0.25d0 ! T/T0
         !real(8),parameter  :: ratio=0.4d0 ! T/T0
   
         real(8),parameter  :: sigma_6=         1.544804416d3    ! Angstrom^6
         real(8),parameter  :: sigma_6_2=2d0  * 1.544804416d3    ! Angstrom^6
         real(8),parameter  :: sigma_12  =      2.386420683693101056d6    ! Angstrom^12
         real(8),parameter  :: sigma_12_2=2d0 * 2.386420683693101056d6    ! Angstrom^12
         real(8),parameter  :: sigma_12_7=7d0 * 2.386420683693101056d6    ! Angstrom^12
         real(8),parameter  :: epsilon24=4.0104d1   
         real(8),parameter  :: epsilon96=4.0104d1 * 4d0
         real(8),parameter  :: epsilon4= 4.0104d1 / 6d0
         real(8), parameter :: del = 0.00001d0
         integer :: iw,ina_i,ina_j, ina

         !integer :: ct0, ct1, count_rate, count_max


         call calc_rad_m2(Nwalker,Natoms,x_pos(1:Nwalker,1:Natoms),&
     &                                   y_pos(1:Nwalker,1:Natoms),&
     &                                   z_pos(1:Nwalker,1:Natoms),&
     &                        rad_m2(1:Nwalker,1:Natoms-1,1:Natoms))


    rad_ma(:,:,:) = rad_m2(:,:,:) ** 4   ! r^-8
    rad_mb(:,:,:) = rad_m2(:,:,:) ** 7   ! r^-14


    VdiffX(:,:) = 0.0d0 
    VdiffY(:,:) = 0.0d0
    VdiffZ(:,:) = 0.0d0

!    DO iw=1, Nwalker

      ! V'_x = 24epsilon * [ sima^6/r^8 - 2*sigma^12/r^14  ] * (x_i -
      ! x_j)  

      DO ina_i=1, Natoms
       DO ina_j=1, Natoms

       IF (ina_j.gt.ina_i) THEN
               VdiffX(:,ina_i) = VdiffX(:,ina_i) + epsilon24 * (sigma_6 * rad_ma(:,ina_i,ina_j)  - sigma_12_2 * rad_mb(:,ina_i,ina_j)) * (x_pos(:,ina_i) - x_pos(:,ina_j))
               VdiffY(:,ina_i) = VdiffY(:,ina_i) + epsilon24 * (sigma_6 * rad_ma(:,ina_i,ina_j)  - sigma_12_2 * rad_mb(:,ina_i,ina_j)) * (y_pos(:,ina_i) - y_pos(:,ina_j))

               VdiffZ(:,ina_i) = VdiffZ(:,ina_i) + epsilon24 * (sigma_6 * rad_ma(:,ina_i,ina_j)  - sigma_12_2 * rad_mb(:,ina_i,ina_j)) * (z_pos(:,ina_i) - z_pos(:,ina_j))
       ENDIF

       IF (ina_j.lt.ina_i) THEN      ! we use symmetry of distance matrix rad_ij = rad_ji
               VdiffX(:,ina_i) = VdiffX(:,ina_i) + epsilon24 * (sigma_6 * rad_ma(:,ina_j,ina_i)  - sigma_12_2 * rad_mb(:,ina_j,ina_i)) * (x_pos(:,ina_i) - x_pos(:,ina_j))
               VdiffY(:,ina_i) = VdiffY(:,ina_i) + epsilon24 * (sigma_6 * rad_ma(:,ina_j,ina_i)  - sigma_12_2 * rad_mb(:,ina_j,ina_i)) * (y_pos(:,ina_i) - y_pos(:,ina_j))

               VdiffZ(:,ina_i) = VdiffZ(:,ina_i) + epsilon24 * (sigma_6 * rad_ma(:,ina_j,ina_i)  - sigma_12_2 * rad_mb(:,ina_j,ina_i)) * (z_pos(:,ina_i) - z_pos(:,ina_j))
       ENDIF

     ENDDO
    ENDDO
!  ENDDO

    dx_pe(:,:) = VdiffX(:,:) 
    dy_pe(:,:) = VdiffY(:,:)
    dz_pe(:,:) = VdiffZ(:,:)



         !call system_clock(ct0, count_rate, count_max)

!         x_pos_p(:,:) = x_pos(:,:) 
!         x_pos_m(:,:) = x_pos(:,:) 
!         y_pos_p(:,:) = y_pos(:,:) 
!         y_pos_m(:,:) = y_pos(:,:) 
!         z_pos_p(:,:) = z_pos(:,:) 
!         z_pos_m(:,:) = z_pos(:,:) 

!         DO ina = 1, Natoms

!          x_pos_p(:,ina) = x_pos_p(:,ina) + del
!          x_pos_m(:,ina) = x_pos_m(:,ina) - del
!          y_pos_p(:,ina) = y_pos_p(:,ina) + del
!          y_pos_m(:,ina) = y_pos_m(:,ina) - del
!          z_pos_p(:,ina) = z_pos_p(:,ina) + del
!          z_pos_m(:,ina) = z_pos_m(:,ina) - del

!          call calc_Udiff0(Nwalker,Natoms,Nregions, Ncycles,x_pos_p, y_pos  , z_pos  ,pe_x_p)
!          call calc_Udiff0(Nwalker,Natoms,Nregions, Ncycles,x_pos_m, y_pos  , z_pos  ,pe_x_m)
!          call calc_Udiff0(Nwalker,Natoms,Nregions, Ncycles,x_pos  , y_pos_p, z_pos  ,pe_y_p)
!          call calc_Udiff0(Nwalker,Natoms,Nregions, Ncycles,x_pos  , y_pos_m, z_pos  ,pe_y_m)
!          call calc_Udiff0(Nwalker,Natoms,Nregions, Ncycles,x_pos  , y_pos  , z_pos_p,pe_z_p)
!          call calc_Udiff0(Nwalker,Natoms,Nregions, Ncycles,x_pos  , y_pos  , z_pos_m,pe_z_m)

!          dx_pe(:,ina) = (pe_x_p(:) - pe_x_m(:))/(2d0*del)
!          dy_pe(:,ina) = (pe_y_p(:) - pe_y_m(:))/(2d0*del)
!          dz_pe(:,ina) = (pe_z_p(:) - pe_z_m(:))/(2d0*del)

!          x_pos_p(:,ina) = x_pos_p(:,ina) - del
!          x_pos_m(:,ina) = x_pos_m(:,ina) + del
!          y_pos_p(:,ina) = y_pos_p(:,ina) - del
!          y_pos_m(:,ina) = y_pos_m(:,ina) + del
!          z_pos_p(:,ina) = z_pos_p(:,ina) - del
!          z_pos_m(:,ina) = z_pos_m(:,ina) + del
!         ENDDO

        !call system_clock(ct1)
        !write(*,*)"sec/Udiff1=",dble(ct1-ct0)/dble(count_rate)

        end subroutine calc_Udiff1

!!!!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX!!!

        subroutine calc_Udiff12(Nwalker,Natoms,Nregions,Ncycles,&
     &    x_pos,y_pos,z_pos,dx_pe,dy_pe, dz_pe,ddx_pe, ddy_pe, ddz_pe)
         use val_mpi_lmp, only : lammps, comm_lammps, ptr
         use bistable1d_val, only : x_pos_mod, y_pos_mod, z_pos_mod

         integer, intent(in) :: Nwalker
         integer, intent(in) :: Natoms
         integer, intent(in) :: Nregions
         integer, intent(in) :: Ncycles
   
         real(8), intent(in) :: x_pos(Nwalker,Natoms)
         real(8), intent(in) :: y_pos(Nwalker,Natoms)
         real(8), intent(in) :: z_pos(Nwalker,Natoms)


         real(8), intent(out) :: dx_pe(Nwalker,Natoms)
         real(8), intent(out) :: dy_pe(Nwalker,Natoms)
         real(8), intent(out) :: dz_pe(Nwalker,Natoms)

         real(8), intent(out) :: ddx_pe(Nwalker,Natoms)
         real(8), intent(out) :: ddy_pe(Nwalker,Natoms)
         real(8), intent(out) :: ddz_pe(Nwalker,Natoms)

         real(8) :: rad_m2(Nwalker,Natoms-1,Natoms)
         real(8) :: rad_ma(Nwalker,Natoms-1,Natoms)
         real(8) :: rad_mb(Nwalker,Natoms-1,Natoms)

         real(8) :: VdiffX(Nwalker,Natoms)
         real(8) :: VdiffY(Nwalker,Natoms)
         real(8) :: VdiffZ(Nwalker,Natoms)

  
         real(8),parameter  :: sigma_6=         1.544804416d3    ! Angstrom^6
         real(8),parameter  :: sigma_6_2=2d0  * 1.544804416d3    ! Angstrom^6
         real(8),parameter  :: sigma_12  =      2.386420683693101056d6    ! Angstrom^12
         real(8),parameter  :: sigma_12_2=2d0 * 2.386420683693101056d6    ! Angstrom^12
         real(8),parameter  :: sigma_12_7=7d0 * 2.386420683693101056d6    ! Angstrom^12
         real(8),parameter  :: epsilon24=4.0104d1   
         real(8),parameter  :: epsilon96=4.0104d1 * 4d0
         real(8),parameter  :: epsilon4= 4.0104d1 / 6d0
         real(8), parameter :: del = 0.00001d0
         integer :: iw,ina_i,ina_j, ina

         !integer :: ct0, ct1, count_rate, count_max

         !call system_clock(ct0, count_rate, count_max)

          ddx_pe(:,:) = 0.0d0
          ddy_pe(:,:) = 0.0d0
          ddz_pe(:,:) = 0.0d0

          dx_pe(:,:) = 0.0d0
          dy_pe(:,:) = 0.0d0
          dz_pe(:,:) = 0.0d0

          VdiffX(:,:) = 0.0d0
          VdiffY(:,:) = 0.0d0
          VdiffZ(:,:) = 0.0d0


          rad_m2(:,:,:) = 0.0d0 

          call calc_rad_m2(Nwalker,Natoms,x_pos(1:Nwalker,1:Natoms),&
                &                                   y_pos(1:Nwalker,1:Natoms),&
                &                                   z_pos(1:Nwalker,1:Natoms),&
                &                        rad_m2(1:Nwalker,1:Natoms-1,1:Natoms))


          rad_ma(:,:,:) = rad_m2(:,:,:) ** 5   ! r^-10
          rad_mb(:,:,:) = rad_m2(:,:,:) ** 8   ! r^-16

          DO iw=1, Nwalker

             ! V''_xx = 24epsilon * [ sima^6/r^8 - 2*sigma^12/r^14  ] + 96epsilon * [
             ! (7*sigma^12/r^16 - 2*sigma^6/r^10) ] * (x_i - x_j)^2 

           DO ina_i=1, Natoms
            DO ina_j=1, Natoms

                IF (ina_j.gt.ina_i) THEN

                       VdiffX(iw,ina_i) = VdiffX(iw,ina_i) + epsilon24 * (sigma_6 * rad_ma(iw,ina_i,ina_j)  - sigma_12_2 * rad_mb(iw,ina_i,ina_j)) / rad_m2(iw,ina_i,ina_j) +&
                        &epsilon96 * (sigma_12_7 * rad_mb(iw,ina_i,ina_j) - sigma_6_2 * rad_ma(iw,ina_i,ina_j) ) * ( (x_pos(iw,ina_i) - x_pos(iw,ina_j))**2 )

                       VdiffY(iw,ina_i) = VdiffY(iw,ina_i) + epsilon24 * (sigma_6 * rad_ma(iw,ina_i,ina_j)  - sigma_12_2 * rad_mb(iw,ina_i,ina_j)) / rad_m2(iw,ina_i,ina_j) +&
                        &epsilon96 * (sigma_12_7 * rad_mb(iw,ina_i,ina_j) - sigma_6_2 * rad_ma(iw,ina_i,ina_j) ) * ( (y_pos(iw,ina_i) - y_pos(iw,ina_j))**2 )

                       VdiffZ(iw,ina_i) = VdiffZ(iw,ina_i) + epsilon24 * (sigma_6 * rad_ma(iw,ina_i,ina_j)  - sigma_12_2 * rad_mb(iw,ina_i,ina_j)) / rad_m2(iw,ina_i,ina_j) +&
                        &epsilon96 * (sigma_12_7 * rad_mb(iw,ina_i,ina_j) - sigma_6_2 * rad_ma(iw,ina_i,ina_j) ) * ( (z_pos(iw,ina_i) - z_pos(iw,ina_j))**2 )

               ENDIF

               IF (ina_j.lt.ina_i) THEN      ! we use symmetry of distance matrix rad_ij = rad_ji

                       VdiffX(iw,ina_i) = VdiffX(iw,ina_i) + epsilon24 * (sigma_6 * rad_ma(iw,ina_j,ina_i)  - sigma_12_2 * rad_mb(iw,ina_j,ina_i)) / rad_m2(iw,ina_j,ina_i) +&
                        &epsilon96 * (sigma_12_7 * rad_mb(iw,ina_j,ina_i) - sigma_6_2 * rad_ma(iw,ina_j,ina_i) ) * ( (x_pos(iw,ina_i) - x_pos(iw,ina_j))**2 )

                       VdiffY(iw,ina_i) = VdiffY(iw,ina_i) + epsilon24 * (sigma_6 * rad_ma(iw,ina_j,ina_i)  - sigma_12_2 * rad_mb(iw,ina_j,ina_i)) / rad_m2(iw,ina_j,ina_i) +&
                        &epsilon96 * (sigma_12_7 * rad_mb(iw,ina_j,ina_i) - sigma_6_2 * rad_ma(iw,ina_j,ina_i) ) * ( (y_pos(iw,ina_i) - y_pos(iw,ina_j))**2 )

                       VdiffZ(iw,ina_i) = VdiffZ(iw,ina_i) + epsilon24 * (sigma_6 * rad_ma(iw,ina_j,ina_i)  - sigma_12_2 * rad_mb(iw,ina_j,ina_i)) / rad_m2(iw,ina_j,ina_i) +&
                        &epsilon96 * (sigma_12_7 * rad_mb(iw,ina_j,ina_i) - sigma_6_2 * rad_ma(iw,ina_j,ina_i) ) * ( (z_pos(iw,ina_i) - z_pos(iw,ina_j))**2 )

               ENDIF

            ENDDO
          ENDDO
         ENDDO

   ddx_pe(:,:) = VdiffX(:,:)
   ddy_pe(:,:) = VdiffY(:,:)
   ddz_pe(:,:) = VdiffZ(:,:)


   ! The calculation of the first derivative of potential
   VdiffX(:,:) = 0.0d0
   VdiffY(:,:) = 0.0d0
   VdiffZ(:,:) = 0.0d0

   call calc_Udiff1(Nwalker,Natoms,Nregions, Ncycles,&
     &                         x_pos,y_pos,z_pos, VdiffX, VdiffY, VdiffZ)


   dx_pe(:,:) = VdiffX(:,:)
   dy_pe(:,:) = VdiffY(:,:)
   dz_pe(:,:) = VdiffZ(:,:)

        !call system_clock(ct1)
        !write(*,*)"sec/Udiff12=",dble(ct1-ct0)/dble(count_rate)

   end subroutine calc_Udiff12
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

        real(8), intent(out) :: rad_minus2(Nwalker,Natoms-1,Natoms)
   
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
       end module potentials
