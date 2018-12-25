      program convert
      implicit none

      character(256), parameter :: fin="lammpstrj.dat"
      character(256), parameter :: fout="atoms.dat"
      integer, parameter :: fp = 20

      integer Natoms
      integer ina
      integer, allocatable :: itype(:)
      real(8), allocatable :: pos(:,:)

      real(8) xlo, xhi, ylo, yhi, zlo, zhi
      real(8) dummy

      open(unit=fp, file=fin, status="unknown")

      read(fp, *)Natoms
      allocate(itype(Natoms), pos(3,Natoms))
      read(fp, *)xlo, xhi
      read(fp, *)ylo, yhi
      read(fp, *)zlo, zhi
      DO ina = 1, Natoms
       read(fp, *)dummy, itype(ina), pos(1:3,ina)
      ENDDO

      close(fp)

      pos(1,:) = xlo + pos(1,:)*(xhi-xlo)
      pos(2,:) = ylo + pos(2,:)*(yhi-ylo)
      pos(3,:) = zlo + pos(3,:)*(zhi-zlo)

      open(unit=fp, file=fout, status="unknown")
      DO ina = 1, Natoms
       write(fp, '(i8, 3f16.8)')itype(ina), pos(1:3,ina)
      ENDDO
      close(fp)
      end program
