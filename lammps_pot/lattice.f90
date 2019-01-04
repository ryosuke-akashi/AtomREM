      subroutine lattice(Nwalker, Natoms, x_pos, y_pos, z_pos,          &
     &                   x_pos_mod, y_pos_mod, z_pos_mod)
       use bistable1d_val, only : a_vec, b_vec, a_orig, bounds

       integer, intent(in) :: Nwalker
       integer, intent(in) :: Natoms
       real(8), intent(in) :: x_pos(Nwalker, Natoms)
       real(8), intent(in) :: y_pos(Nwalker, Natoms)
       real(8), intent(in) :: z_pos(Nwalker, Natoms)
       real(8), intent(out) :: x_pos_mod(Nwalker, Natoms)
       real(8), intent(out) :: y_pos_mod(Nwalker, Natoms)
       real(8), intent(out) :: z_pos_mod(Nwalker, Natoms)

       ! internal variables
       !real(8) x_pos_tmp(Nwalker, Natoms)
       !real(8) y_pos_tmp(Nwalker, Natoms)
       !real(8) z_pos_tmp(Nwalker, Natoms)
       real(8) n1(Nwalker, Natoms)
       real(8) n2(Nwalker, Natoms)
       real(8) n3(Nwalker, Natoms)
       integer ina
       !/internal variables

       x_pos_mod(:,:) = x_pos(:,:) - a_orig(1)
       y_pos_mod(:,:) = y_pos(:,:) - a_orig(2)
       z_pos_mod(:,:) = z_pos(:,:) - a_orig(3)

       ! vec = n1*a1 +n2*a2 + n3*a3

       DO ina = 1, Natoms
        n1(1:Nwalker,ina) = x_pos_mod(1:Nwalker, ina)*b_vec(1,1)        &
     &                    + y_pos_mod(1:Nwalker, ina)*b_vec(1,2)        &
     &                    + z_pos_mod(1:Nwalker, ina)*b_vec(1,3)
        n2(1:Nwalker,ina) = x_pos_mod(1:Nwalker, ina)*b_vec(2,1)        &
     &                    + y_pos_mod(1:Nwalker, ina)*b_vec(2,2)        &
     &                    + z_pos_mod(1:Nwalker, ina)*b_vec(2,3)
        n3(1:Nwalker,ina) = x_pos_mod(1:Nwalker, ina)*b_vec(3,1)        &
     &                    + y_pos_mod(1:Nwalker, ina)*b_vec(3,2)        &
     &                    + z_pos_mod(1:Nwalker, ina)*b_vec(3,3)
       ENDDO

       ! pull back the positions within the given unitcell

       !n1(1:Nwalker, 1:Natoms) = 0d0
       !n2(1:Nwalker, 1:Natoms) = 0d0
       !n3(1:Nwalker, 1:Natoms) = 0d0

       IF(bounds(1).eq."p")THEN
        n1(1:Nwalker, 1:Natoms) = dble(floor(n1(1:Nwalker, 1:Natoms)))
       ENDIF

       IF(bounds(2).eq."p")THEN
        n2(1:Nwalker, 1:Natoms) = dble(floor(n2(1:Nwalker, 1:Natoms)))
       ENDIF

       IF(bounds(2).eq."p")THEN
        n3(1:Nwalker, 1:Natoms) = dble(floor(n3(1:Nwalker, 1:Natoms)))
       ENDIF


       x_pos_mod(1:Nwalker,1:Natoms) = x_pos(1:Nwalker,1:Natoms)        &
     &                           - n1(1:Nwalker,1:Natoms)*a_vec(1,1)    &
     &                           - n2(1:Nwalker,1:Natoms)*a_vec(2,1)    &
     &                           - n3(1:Nwalker,1:Natoms)*a_vec(3,1) 

       y_pos_mod(1:Nwalker,1:Natoms) = y_pos(1:Nwalker,1:Natoms)        &
     &                           - n1(1:Nwalker,1:Natoms)*a_vec(1,2)    &
     &                           - n2(1:Nwalker,1:Natoms)*a_vec(2,2)    &
     &                           - n3(1:Nwalker,1:Natoms)*a_vec(3,2) 

       z_pos_mod(1:Nwalker,1:Natoms) = z_pos(1:Nwalker,1:Natoms)        &
     &                           - n1(1:Nwalker,1:Natoms)*a_vec(1,3)    &
     &                           - n2(1:Nwalker,1:Natoms)*a_vec(2,3)    &
     &                           - n3(1:Nwalker,1:Natoms)*a_vec(3,3) 


      end subroutine
   
