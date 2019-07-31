
      module potentials
       implicit none

       interface calc_Udiff
        module procedure calc_Udiff0,calc_Udiff0_atom, calc_Udiff1,calc_Udiff12
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

         real(8) :: pos(3*Nregions*Natoms)
         real(8) :: pe_tmp(Natoms) 

         integer :: iw,ina,idx, idx2, icy
         character(16) :: id_pe
         integer :: len_id_pe


          id_pe = "peatom"
          len_id_pe = len(trim(id_pe))

          pe_tmp(:) = 0d0

         DO icy = 1, Nwalker

           DO ina = 1, Natoms

             idx  = 3*(ina-1)

             pos(idx+1) = x_pos(icy,ina) 
             pos(idx+2) = y_pos(icy,ina)
             pos(idx+3) = z_pos(icy,ina)

           ENDDO

          call lammps_scatter_atoms(ptr, 'x', 1, pos(:))
          call lammps_command(ptr, 'run 0', 5)
          call lammps_extract_compute_atom_vec(ptr, trim(id_pe), len_id_pe, pe_tmp(:), Natoms)


          pe(icy) = sum(pe_tmp(1:Natoms))


         ENDDO !icy


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

         real(8) :: pos(3*Natoms)
         real(8) :: pe_tmp(Natoms)
         
         integer :: iw,ina,idx, idx2, icy
         character(16) :: id_pe
         
         integer :: len_id_pe

         id_pe = "peatom"
         len_id_pe = len(trim(id_pe))

         pe_tmp(:) = 0d0


         DO icy = 1, Nwalker

           DO ina = 1, Natoms

             idx  = 3*(ina-1)

             pos(idx+1) = x_pos(icy,ina)
             pos(idx+2) = y_pos(icy,ina)
             pos(idx+3) = z_pos(icy,ina)

           ENDDO

          call lammps_scatter_atoms(ptr, 'x', 1, pos(:))
          call lammps_command(ptr, 'run 0', 5)
          call lammps_extract_compute_atom_vec(ptr, trim(id_pe), len_id_pe, pe_tmp(:), Natoms)

          pe_atom(icy,1:Natoms) = pe_tmp(1:Natoms)


         ENDDO !icy

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

         real(8) :: pos(3*Natoms)

         real(8) :: f_x(Natoms)
         real(8) :: f_y(Natoms)
         real(8) :: f_z(Natoms)

         integer :: iw,ina_i,ina_j, ina,  idx, idx2, icy
         character(16) :: id_X, id_Y, id_Z
         
         integer :: len_id_X, len_id_Y, len_id_Z




         id_X = "frcX"
         id_Y = "frcY"
         id_Z = "frcZ"
         len_id_X = len(trim(id_X))
         len_id_Y = len(trim(id_Y))
         len_id_Z = len(trim(id_Z))


         dx_pe(1:Nwalker,1:Natoms) = 0.0d0
         dy_pe(1:Nwalker,1:Natoms) = 0.0d0
         dz_pe(1:Nwalker,1:Natoms) = 0.0d0
         
         f_x(1:Natoms) = 0.0d0
         f_y(1:Natoms) = 0.0d0
         f_z(1:Natoms) = 0.0d0


         DO icy = 1, Nwalker

           DO ina = 1, Natoms

             idx  = 3*(ina-1)

             pos(idx+1) = x_pos(icy,ina)
             pos(idx+2) = y_pos(icy,ina)
             pos(idx+3) = z_pos(icy,ina)

           ENDDO

          call lammps_scatter_atoms(ptr, 'x', 1, pos(:))
          call lammps_command(ptr, 'run 0', 5)
          call lammps_extract_compute_atom_vec(ptr, trim(id_X), len_id_X, f_x(:), Natoms)
          call lammps_extract_compute_atom_vec(ptr, trim(id_Y), len_id_Y, f_y(:), Natoms)
          call lammps_extract_compute_atom_vec(ptr, trim(id_Z), len_id_Z, f_z(:), Natoms)

          dx_pe(icy,1:Natoms) = -1.0d0 * f_x(1:Natoms)
          dy_pe(icy,1:Natoms) = -1.0d0 * f_y(1:Natoms)
          dz_pe(icy,1:Natoms) = -1.0d0 * f_z(1:Natoms)

         ENDDO !icy


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

         real(8) :: pos(3*Natoms)

         real(8) :: f_x(Natoms), hs_x(Natoms)
         real(8) :: f_y(Natoms), hs_y(Natoms)
         real(8) :: f_z(Natoms), hs_z(Natoms)

         integer :: iw,ina_i,ina_j, ina,  idx, idx2, icy
         character(16) :: id_X, id_Y, id_Z, id_hes
        
         integer :: len_id_X, len_id_Y, len_id_Z, len_id_hes


         id_X = "frcX"
         id_Y = "frcY"
         id_Z = "frcZ"
         len_id_X = len(trim(id_X))
         len_id_Y = len(trim(id_Y))
         len_id_Z = len(trim(id_Z))


         id_hes = "hsmtrxD"



         len_id_hes = len(trim(id_hes))

         dx_pe(1:Nwalker,1:Natoms) = 0.0d0
         dy_pe(1:Nwalker,1:Natoms) = 0.0d0
         dz_pe(1:Nwalker,1:Natoms) = 0.0d0

         f_x(1:Natoms) = 0.0d0
         f_y(1:Natoms) = 0.0d0
         f_z(1:Natoms) = 0.0d0


         DO icy = 1, Nwalker

           DO ina = 1, Natoms

             idx  = 3*(ina-1)

             pos(idx+1) = x_pos(icy,ina)
             pos(idx+2) = y_pos(icy,ina)
             pos(idx+3) = z_pos(icy,ina)

           ENDDO

          call lammps_scatter_atoms(ptr, 'x', 1, pos(:))
          call lammps_command(ptr, 'run 0', 5)

          call lammps_extract_compute_atom_vec(ptr, trim(id_X), len_id_X, f_x(:), Natoms)
          call lammps_extract_compute_atom_vec(ptr, trim(id_Y), len_id_Y, f_y(:), Natoms)
          call lammps_extract_compute_atom_vec(ptr, trim(id_Z), len_id_Z, f_z(:), Natoms)

          dx_pe(icy,1:Natoms) = -1.0d0 * f_x(1:Natoms)
          dy_pe(icy,1:Natoms) = -1.0d0 * f_y(1:Natoms)
          dz_pe(icy,1:Natoms) = -1.0d0 * f_z(1:Natoms)

!          call lammps_extract_variable_atom_vec(ptr, trim(id_hes), len_id_hes, hs_x(:), Natoms)

          call lammps_extract_compute_hessian(ptr, trim(id_hes), len_id_hes, hs_x(:), hs_y(:), hs_z(:), Natoms )

          ddx_pe(icy,1:Natoms) = hs_x(1:Natoms)
          ddy_pe(icy,1:Natoms) = hs_y(1:Natoms)
          ddz_pe(icy,1:Natoms) = hs_z(1:Natoms)

         ENDDO !icy

        end subroutine calc_Udiff12

!!!!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX!!!
!------------------------  
       end module potentials 



