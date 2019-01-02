       program make_data
       ! a small program to make FCC or HCP structures in Hex supercell.
       implicit none

       real(8), parameter :: bond=1.08736d0  ! interatomic distance
       integer, parameter :: repeat_x=3 ! must be larger than 1
       integer, parameter :: repeat_y=3 ! must be larger than 1
       integer, parameter :: repeat_z=6 ! must be multiple of 6
       real(8), parameter :: a_orig(3) = (/-0.5, -0.5, -0.5/)
       character(3), parameter :: latt="fcc" ! "fcc" or "hcp"

       real(8) pos(3, repeat_x*repeat_y*repeat_z)
       real(8) xlo, xhi
       real(8) ylo, yhi
       real(8) zlo, zhi
       real(8) xy, xz, yz
       integer Natoms
       integer ind, iat
       integer ix, iy, iz

       !! generate the lattice parameters
       xlo = a_orig(1)
       ylo = a_orig(2)
       zlo = a_orig(3)
       xhi = xlo + bond*repeat_x
       yhi = ylo + bond*repeat_y*dsqrt(3d0)*0.5d0
       zhi = zlo + bond*dsqrt(2d0/3d0)*repeat_z

       xy = bond*0.5d0*mod((2*repeat_x-repeat_y),2)
       xz = 0d0
       yz = 0d0
       !!/generate the lattice parameters

       !! generate atomic positions
       Natoms = repeat_x*repeat_y*repeat_z

       ind = 0
       DO iz = 1, repeat_z
       DO iy = 1, repeat_y
       DO ix = 1, repeat_x
        ind = ind+1
        pos(1,ind) = bond*(ix-1) + 0.5d0*bond*(mod((iy-1),2))
        pos(2,ind) = bond*(iy-1)*dsqrt(3d0)*0.5d0             &
     &                    + dsqrt(1d0/3d0)*bond*((mod((iz-1),3)))
        pos(3,ind) = bond*(iz-1)*dsqrt(2d0/3d0)
       ENDDO
       ENDDO
       ENDDO
       pos(1,:) = pos(1,:)- a_orig(1)
       pos(2,:) = pos(2,:)- a_orig(2)
       pos(3,:) = pos(3,:)- a_orig(3)
       write(*,*)"ind=", ind, "Natoms=", Natoms


       write(*,'(2f12.6,a9)')xlo, xhi, "xlo xhi"
       write(*,'(2f12.6,a9)')ylo, yhi, "ylo yhi"
       write(*,'(2f12.6,a9)')zlo, zhi, "zlo zhi"
       write(*,'(3f12.6,a10)')xy, xz, yz, "xy xz yz"
       open(unit=10, file="atoms.dat", status="unknown")
       DO iat = 1, Natoms
        write(10,'(i10, 3f15.6)') 1, pos(:,iat)
       ENDDO
       close(10)
       
       !!/generate atomic positions

       end program
