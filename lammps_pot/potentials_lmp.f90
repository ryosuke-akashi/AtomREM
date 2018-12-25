      module potentials
       implicit none
       
       interface calc_Udiff
        module procedure calc_Udiff0,calc_Udiff0_atom, calc_Udiff1, calc_Udiff12
       end interface calc_Udiff

       !interface calc_Vdiff
       ! module procedure calc_Vdiff0, calc_Vdiff1, calc_Vdiff12
       !end interface calc_Vdiff

       contains

        subroutine calc_Udiff0(Nwalker,Natoms,Nregions,Ncycles,x_pos, y_pos, z_pos,  pe)
         use val_mpi_lmp, only : lammps, comm_lammps, ptr
         integer, intent(in) :: Nwalker
         integer, intent(in) :: Natoms
         integer, intent(in) :: Nregions
         integer, intent(in) :: Ncycles
   
         real(8), intent(in) :: x_pos(Nwalker,Natoms)
         real(8), intent(in) :: y_pos(Nwalker,Natoms)
         real(8), intent(in) :: z_pos(Nwalker,Natoms)

         real(8), intent(out) :: pe(Nwalker)

         real(8) :: pos(3*Nregions*Natoms)
         real(8) :: pe_tmp(Nregions)
         !real(8),parameter  :: ratio=0.02d0 ! T/T0
         !real(8),parameter  :: ratio=0.1d0 ! T/T0
         !real(8),parameter  :: ratio=0.25d0 ! T/T0
         !real(8),parameter  :: ratio=0.4d0 ! T/T0
   
         integer :: iw,ina,idx, idx2, icy
         character(16) :: id_pe(Nregions)
         character(16) :: ich
         integer :: len_id_pe(Nregions)

         real(8), parameter :: dcut2 = 40.0d0 
         !! must be consistent with dcut2 in main.f90

         DO iw = 1, Nregions
          write(ich,'(i0)')iw
          id_pe(iw) = "pe" //  trim(ich)
          len_id_pe(iw) = len(trim(id_pe(iw)))
         ENDDO

         DO icy = 1, Ncycles
          

          DO iw = 1, Nregions
           DO ina = 1, Natoms
            idx  = (iw -1)*Natoms*3 + 3*(ina-1)
            idx2 = (icy-1)*Nregions + iw
            !write(*,*) "idx, idx2", idx, idx2
            IF(idx2.gt.Nwalker) THEN  !! pudding
             pos(idx+1) = x_pos(Nwalker,ina) + dcut2*dble(iw-1) !! test
             pos(idx+2) = y_pos(Nwalker,ina)
             pos(idx+3) = z_pos(Nwalker,ina)
            ELSE
             pos(idx+1) = x_pos(idx2,ina) + dcut2*dble(iw-1) !! test
             pos(idx+2) = y_pos(idx2,ina)
             pos(idx+3) = z_pos(idx2,ina)
            ENDIF
           ENDDO
          ENDDO
          call lammps_scatter_atoms(ptr, 'x', 1, pos(:))
          call lammps_command(ptr, 'run 0', 5)
          DO iw = 1, Nregions
           call lammps_extract_compute(ptr,trim(id_pe(iw)) , len_id_pe(iw), pe_tmp(iw))
          ENDDO

          DO iw = 1, Nregions
           idx2 = (icy-1)*Nregions + iw
           IF(idx2.gt.Nwalker) return  !! pudding
           pe(idx2) = pe_tmp(iw)
          ENDDO

         ENDDO !icy


        end subroutine calc_Udiff0


        subroutine calc_Udiff0_atom(Nwalker,Natoms,Nregions,Ncycles,    &
     &              x_pos, y_pos, z_pos,  pe_atom)
         use val_mpi_lmp, only : lammps, comm_lammps, ptr
         integer, intent(in) :: Nwalker
         integer, intent(in) :: Natoms
         integer, intent(in) :: Nregions
         integer, intent(in) :: Ncycles
   
         real(8), intent(in) :: x_pos(Nwalker,Natoms)
         real(8), intent(in) :: y_pos(Nwalker,Natoms)
         real(8), intent(in) :: z_pos(Nwalker,Natoms)

         real(8), intent(out) :: pe_atom(Nwalker,Natoms)

         real(8) :: pos(3*Nregions*Natoms)
         real(8) :: pe_tmp(Natoms*Nregions)
         !real(8),parameter  :: ratio=0.02d0 ! T/T0
         !real(8),parameter  :: ratio=0.1d0 ! T/T0
         !real(8),parameter  :: ratio=0.25d0 ! T/T0
         !real(8),parameter  :: ratio=0.4d0 ! T/T0
   
         integer :: iw,ina,idx, idx2, icy
         character(16) :: id_pe
         character(16) :: ich
         integer :: len_id_pe

         real(8), parameter :: dcut2 = 40.0d0 
         !! must be consistent with dcut2 in main.f90

         id_pe = "peatom"
         len_id_pe = len(trim(id_pe))

         pe_tmp(:) = 0d0

         write(*,*)"id_pe=", id_pe(:)
         DO icy = 1, Ncycles
          

          DO iw = 1, Nregions
           DO ina = 1, Natoms
            idx  = (iw -1)*Natoms*3 + 3*(ina-1)
            idx2 = (icy-1)*Nregions + iw
            !write(*,*) "idx, idx2", idx, idx2
            IF(idx2.gt.Nwalker) THEN  !! pudding
             pos(idx+1) = x_pos(Nwalker,ina) + dcut2*dble(iw-1) !! test
             pos(idx+2) = y_pos(Nwalker,ina)
             pos(idx+3) = z_pos(Nwalker,ina)
            ELSE
             pos(idx+1) = x_pos(idx2,ina) + dcut2*dble(iw-1) !! test
             pos(idx+2) = y_pos(idx2,ina)
             pos(idx+3) = z_pos(idx2,ina)
            ENDIF
           ENDDO
          ENDDO
          call lammps_scatter_atoms(ptr, 'x', 1, pos(:))
          call lammps_command(ptr, 'run 0', 5)
          call lammps_extract_compute_atom_vec(ptr,trim(id_pe) , len_id_pe, pe_tmp(:), Natoms*Nregions)


          !write(*,*)"pe_tmp",pe_tmp(:)
          !stop

          DO iw = 1, Nregions
           idx2 = (icy-1)*Nregions + iw
           IF(idx2.gt.Nwalker) return  !! pudding
           DO ina = 1, Natoms
            pe_atom(idx2,ina) = pe_tmp((iw-1)*Natoms+ina)
           ENDDO
          ENDDO

         ENDDO !icy
         !DO iw = 1, Nwalker
         ! DO ina = 1, Natoms
         !  write(*,*)"pe_atom",pe_atom(iw,ina)
         ! ENDDO
         !ENDDO
         !stop

        end subroutine calc_Udiff0_atom

        subroutine calc_Udiff1(Nwalker,Natoms,Nregions, Ncycles,         &
     &                         x_pos,y_pos,z_pos,dx_pe,dy_pe, dz_pe)
         use val_mpi_lmp, only : lammps, comm_lammps, ptr
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

         real(8) :: pos(3*Nregions*Natoms)
         real(8) :: pe_x_p(Nregions)
         real(8) :: pe_x_m(Nregions)
         real(8) :: pe_y_p(Nregions)
         real(8) :: pe_y_m(Nregions)
         real(8) :: pe_z_p(Nregions)
         real(8) :: pe_z_m(Nregions)


         !real(8),parameter  :: ratio=0.02d0 ! T/T0
         !real(8),parameter  :: ratio=0.1d0 ! T/T0
         !real(8),parameter  :: ratio=0.25d0 ! T/T0
         !real(8),parameter  :: ratio=0.4d0 ! T/T0
   
         real(8), parameter :: del = 0.00001d0
         integer :: iw,ina_i,ina_j, ina,  idx, idx2, icy
         character(16) :: id_pe(Nregions)
         character(16) :: ich
         integer :: len_id_pe(Nregions)

         integer :: ct0, ct1, ct2, ct3, count_rate, count_max

         real(8), parameter :: dcut2 = 40.0d0 
         !! must be consistent with dcut2 in main.f90


         !call system_clock(ct0, count_rate, count_max)

         DO iw = 1, Nregions
          write(ich,'(i0)')iw
          id_pe(iw) = "pe" //  trim(ich)
          len_id_pe(iw) = len(trim(id_pe(iw)))
         ENDDO

         DO icy = 1, Ncycles
          

          DO iw = 1, Nregions
           DO ina = 1, Natoms
            idx  = (iw -1)*Natoms*3 + 3*(ina-1)
            idx2 = (icy-1)*Nregions + iw
            !write(*,*) "idx, idx2", idx, idx2
            IF(idx2.gt.Nwalker) THEN  !! pudding
             pos(idx+1) = x_pos(Nwalker,ina) + dcut2*dble(iw-1) !! test
             pos(idx+2) = y_pos(Nwalker,ina)
             pos(idx+3) = z_pos(Nwalker,ina)
            ELSE
             pos(idx+1) = x_pos(idx2,ina) + dcut2*dble(iw-1) !! test
             pos(idx+2) = y_pos(idx2,ina)
             pos(idx+3) = z_pos(idx2,ina)
            ENDIF
           ENDDO
          ENDDO



          DO ina = 1, Natoms

           DO iw = 1, Nregions
            idx = (iw-1)*Natoms*3 + 3*(ina-1)
            pos(idx+1) = pos(idx+1) + del
           ENDDO

           call lammps_scatter_atoms(ptr, 'x', 1, pos(:))
           call lammps_command(ptr, 'run 0', 5)
           DO iw = 1, Nregions
            call lammps_extract_compute(ptr,trim(id_pe(iw)) , len_id_pe(iw), pe_x_p(iw))
           ENDDO

           DO iw = 1, Nregions
            idx = (iw-1)*Natoms*3 + 3*(ina-1)
            pos(idx+1) = pos(idx+1) - 2d0*del
           ENDDO

           call lammps_scatter_atoms(ptr, 'x', 1, pos(:))
           call lammps_command(ptr, 'run 0', 5)
           DO iw = 1, Nregions
            call lammps_extract_compute(ptr,trim(id_pe(iw)) , len_id_pe(iw), pe_x_m(iw))
           ENDDO

           DO iw = 1, Nregions
            idx = (iw-1)*Natoms*3 + 3*(ina-1)
            pos(idx+1) = pos(idx+1) + del
           ENDDO


           DO iw = 1, Nregions
            idx = (iw-1)*Natoms*3 + 3*(ina-1)
            pos(idx+2) = pos(idx+2) + del
           ENDDO
           call lammps_scatter_atoms(ptr, 'x', 1, pos(:))
           call lammps_command(ptr, 'run 0', 5)
           DO iw = 1, Nregions
            call lammps_extract_compute(ptr,trim(id_pe(iw)) , len_id_pe(iw), pe_y_p(iw))
           ENDDO
           DO iw = 1, Nregions
            idx = (iw-1)*Natoms*3 + 3*(ina-1)
            pos(idx+2) = pos(idx+2) - 2d0*del
           ENDDO
           call lammps_scatter_atoms(ptr, 'x', 1, pos(:))
           call lammps_command(ptr, 'run 0', 5)
           DO iw = 1, Nregions
            call lammps_extract_compute(ptr,trim(id_pe(iw)) , len_id_pe(iw), pe_y_m(iw))
           ENDDO
           DO iw = 1, Nregions
            idx = (iw-1)*Natoms*3 + 3*(ina-1)
            pos(idx+2) = pos(idx+2) + del
           ENDDO

           DO iw = 1, Nregions
            idx = (iw-1)*Natoms*3 + 3*(ina-1)
            pos(idx+3) = pos(idx+3) + del
           ENDDO
           call lammps_scatter_atoms(ptr, 'x', 1, pos(:))
           call lammps_command(ptr, 'run 0', 5)
           DO iw = 1, Nregions
            call lammps_extract_compute(ptr,trim(id_pe(iw)) , len_id_pe(iw), pe_z_p(iw))
           ENDDO
           DO iw = 1, Nregions
            idx = (iw-1)*Natoms*3 + 3*(ina-1)
            pos(idx+3) = pos(idx+3) - 2d0*del
           ENDDO
           call lammps_scatter_atoms(ptr, 'x', 1, pos(:))
           call lammps_command(ptr, 'run 0', 5)
           DO iw = 1, Nregions
            call lammps_extract_compute(ptr,trim(id_pe(iw)) , len_id_pe(iw), pe_z_m(iw))
           ENDDO
           DO iw = 1, Nregions
            idx = (iw-1)*Natoms*3 + 3*(ina-1)
            pos(idx+3) = pos(idx+3) + del
           ENDDO

           DO iw = 1, Nregions
            idx2 = (icy-1)*Nregions + iw
            IF(idx2.gt.Nwalker) return  !! pudding
            dx_pe(idx2,ina) = (pe_x_p(iw) - pe_x_m(iw))/(2d0*del)
            dy_pe(idx2,ina) = (pe_y_p(iw) - pe_y_m(iw))/(2d0*del)
            dz_pe(idx2,ina) = (pe_z_p(iw) - pe_z_m(iw))/(2d0*del)
           ENDDO
          ENDDO ! ina
         ENDDO ! icy

         !write(*,*)"pe=", pe_x_p(:)

         !call system_clock(ct1)
         !write(1002,*)"sec/allwalkers=",dble(ct1-ct0)/dble(count_rate)
         !stop


        end subroutine calc_Udiff1

        subroutine calc_Udiff12(Nwalker,Natoms,Nregions,Ncycles,&
     &    x_pos,y_pos,z_pos,dx_pe,dy_pe, dz_pe,ddx_pe, ddy_pe, ddz_pe)
         use val_mpi_lmp, only : lammps, comm_lammps, ptr
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


         real(8) :: pos(3*Nregions*Natoms)
         real(8) :: pe(Nregions)
         real(8) :: pe_x_p(Nregions), pe_y_p(Nregions), pe_z_p(Nregions)
         real(8) :: pe_x_m(Nregions), pe_y_m(Nregions), pe_z_m(Nregions)

         !real(8),parameter  :: ratio=0.02d0 ! T/T0
         !real(8),parameter  :: ratio=0.1d0 ! T/T0
         !real(8),parameter  :: ratio=0.25d0 ! T/T0
         !real(8),parameter  :: ratio=0.4d0 ! T/T0
   
         real(8), parameter :: del = 0.00001d0
         integer :: iw,ina_i,ina_j, ina, idx, idx2, icy
         character(16) :: id_pe(Nregions)
         character(16) :: ich
         integer :: len_id_pe(Nregions)

         real(8), parameter :: dcut2 = 40.0d0 
         !! must be consistent with dcut2 in main.f90

         !integer :: ct0, ct1, count_rate, count_max

         !call system_clock(ct0, count_rate, count_max)
         DO iw = 1, Nregions
          write(ich,'(i0)')iw
          id_pe(iw) = "pe" //  trim(ich)
          len_id_pe(iw) = len(trim(id_pe(iw)))
         ENDDO

         DO icy = 1, Ncycles
          

          DO iw = 1, Nregions
           DO ina = 1, Natoms
            idx  = (iw -1)*Natoms*3 + 3*(ina-1)
            idx2 = (icy-1)*Nregions + iw
            !write(*,*) "idx, idx2", idx, idx2
            IF(idx2.gt.Nwalker) THEN  !! pudding
             pos(idx+1) = x_pos(Nwalker,ina) + dcut2*dble(iw-1) !! test
             pos(idx+2) = y_pos(Nwalker,ina)
             pos(idx+3) = z_pos(Nwalker,ina)
            ELSE
             pos(idx+1) = x_pos(idx2,ina) + dcut2*dble(iw-1) !! test
             pos(idx+2) = y_pos(idx2,ina)
             pos(idx+3) = z_pos(idx2,ina)
            ENDIF
           ENDDO
          ENDDO
          call lammps_scatter_atoms(ptr, 'x', 1, pos(:))
          call lammps_command(ptr, 'run 0', 5)
          DO iw = 1, Nregions
           call lammps_extract_compute(ptr,trim(id_pe(iw)) , len_id_pe(iw), pe(iw))
          ENDDO



          DO ina = 1, Natoms

           DO iw = 1, Nregions
            idx = (iw-1)*Natoms*3 + 3*(ina-1)
            pos(idx+1) = pos(idx+1) + del
           ENDDO

           call lammps_scatter_atoms(ptr, 'x', 1, pos(:))
           call lammps_command(ptr, 'run 0', 5)
           DO iw = 1, Nregions
            call lammps_extract_compute(ptr,trim(id_pe(iw)) , len_id_pe(iw), pe_x_p(iw))
           ENDDO

           DO iw = 1, Nregions
            idx = (iw-1)*Natoms*3 + 3*(ina-1)
            pos(idx+1) = pos(idx+1) - 2d0*del
           ENDDO

           call lammps_scatter_atoms(ptr, 'x', 1, pos(:))
           call lammps_command(ptr, 'run 0', 5)
           DO iw = 1, Nregions
            call lammps_extract_compute(ptr,trim(id_pe(iw)) , len_id_pe(iw), pe_x_m(iw))
           ENDDO

           DO iw = 1, Nregions
            idx = (iw-1)*Natoms*3 + 3*(ina-1)
            pos(idx+1) = pos(idx+1) + del
           ENDDO


           DO iw = 1, Nregions
            idx = (iw-1)*Natoms*3 + 3*(ina-1)
            pos(idx+2) = pos(idx+2) + del
           ENDDO
           call lammps_scatter_atoms(ptr, 'x', 1, pos(:))
           call lammps_command(ptr, 'run 0', 5)
           DO iw = 1, Nregions
            call lammps_extract_compute(ptr,trim(id_pe(iw)) , len_id_pe(iw), pe_y_p(iw))
           ENDDO
           DO iw = 1, Nregions
            idx = (iw-1)*Natoms*3 + 3*(ina-1)
            pos(idx+2) = pos(idx+2) - 2d0*del
           ENDDO
           call lammps_scatter_atoms(ptr, 'x', 1, pos(:))
           call lammps_command(ptr, 'run 0', 5)
           DO iw = 1, Nregions
            call lammps_extract_compute(ptr,trim(id_pe(iw)) , len_id_pe(iw), pe_y_m(iw))
           ENDDO
           DO iw = 1, Nregions
            idx = (iw-1)*Natoms*3 + 3*(ina-1)
            pos(idx+2) = pos(idx+2) + del
           ENDDO

           DO iw = 1, Nregions
            idx = (iw-1)*Natoms*3 + 3*(ina-1)
            pos(idx+3) = pos(idx+3) + del
           ENDDO
           call lammps_scatter_atoms(ptr, 'x', 1, pos(:))
           call lammps_command(ptr, 'run 0', 5)
           DO iw = 1, Nregions
            call lammps_extract_compute(ptr,trim(id_pe(iw)) , len_id_pe(iw), pe_z_p(iw))
           ENDDO
           DO iw = 1, Nregions
            idx = (iw-1)*Natoms*3 + 3*(ina-1)
            pos(idx+3) = pos(idx+3) - 2d0*del
           ENDDO
           call lammps_scatter_atoms(ptr, 'x', 1, pos(:))
           call lammps_command(ptr, 'run 0', 5)
           DO iw = 1, Nregions
            call lammps_extract_compute(ptr,trim(id_pe(iw)) , len_id_pe(iw), pe_z_m(iw))
           ENDDO
           DO iw = 1, Nregions
            idx = (iw-1)*Natoms*3 + 3*(ina-1)
            pos(idx+3) = pos(idx+3) + del
           ENDDO

           DO iw = 1, Nregions
            idx2 = (icy-1)*Nregions + iw
            IF(idx2.gt.Nwalker) return  !! pudding
            dx_pe(idx2,ina) = (pe_x_p(iw) - pe_x_m(iw))/(2d0*del)
            dy_pe(idx2,ina) = (pe_y_p(iw) - pe_y_m(iw))/(2d0*del)
            dz_pe(idx2,ina) = (pe_z_p(iw) - pe_z_m(iw))/(2d0*del)
            ddx_pe(idx2,ina) = (pe_x_p(iw) -2d0*pe(iw) + pe_x_m(iw))/(del*del)
            ddy_pe(idx2,ina) = (pe_y_p(iw) -2d0*pe(iw) + pe_y_m(iw))/(del*del)
            ddz_pe(idx2,ina) = (pe_z_p(iw) -2d0*pe(iw) + pe_z_m(iw))/(del*del)
           ENDDO
          ENDDO ! ina
         ENDDO ! icy

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
