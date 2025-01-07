! This program perform bestfit two pdb strcturues
! based on all atoms
! You can specify the atom range to be fited.

program axis_rotation
  use bestfit
  use quaternion
  use pdb

  implicit none

  integer i,j,N ! num of atoms

!  type (PDB_atom), allocatable, dimension(:) :: PA
  double precision, allocatable, dimension(:,:) :: coord
!  double precision :: coord(3)
  double precision, allocatable, dimension(:) :: masses

  double precision :: temp(3)
  
  character(30) line

  double precision t, x(3), l, pi, center(3)
  double precision :: q(0:3)

  character(50)  pdbfilename

  character(30) :: argv, programname, d
  integer :: argc, iargc
  
  argc=iargc()
  call getarg(0,programname)
  if (argc < 5) then
     call usage(programname)
  else
     call getarg(1,pdbfilename)
     call getarg(2,d)
     read(d,*),x(1)
     call getarg(3,d)
     read(d,*),x(2)
     call getarg(4,d)
     read(d,*),x(3)
     call getarg(5,d)
     read(d,*),t
     call getarg(6,d)
     read(d,*),N
!     call getarg(6,d)
!     read(d,*),center(1)
!     call getarg(7,d)
!     read(d,*),center(2)
!     call getarg(8,d)
!     read(d,*),center(3)
  end if

  allocate(coord(3,1))

  open(17,file=pdbfilename,status='old')
  do i=1,N,1
     read(17,'(A30,3F8.3)') line,coord(1,1),coord(2,1),coord(3,1)
     center(1)=center(1)+coord(1,1)
     center(2)=center(2)+coord(2,1)
     center(3)=center(3)+coord(3,1)
  end do
  center(1)=center(1)/N
  center(2)=center(2)/N
  center(3)=center(3)/N
  !  call read_PDB (17, PA, N)
  close(17)

  pi=acos(-1.0d0)
  
  t=t/180.0d0*pi
  q(0) = cos(t*0.5d0)

  do i=1,3,1
     l = l + x(i)**2
  end do
  l = sqrt(l)
  
  do i=1,3,1
     q(i) = x(i)/l
  end do

  do i=1,3,1
     q(i) = sin(t*0.5d0)*q(i)
  end do

  open(17,file=pdbfilename,status='old')
  do i=1,N,1
     read(17,'(A30,3F8.3)') line,coord(1,1),coord(2,1),coord(3,1)

!     call com_shift(N, coord, masses, center)
!     do i = 1, N
     coord(:, 1) = coord(:, 1) - center(:)
!     end do

!     do i = 1, N
     call rotate(coord(:, 1), q, temp)
     coord(:, 1) = temp(:)
!    end do

!     call com_unshift(N, coord, masses, center)
!     do i = 1, N
     coord(:, 1) = coord(:, 1) + center(:)
!     end do

     write(*,'(A30,3F8.2,1X)') line,coord(1,1),coord(2,1),coord(3,1)
  end do

  close(17)
  
!  allocate(masses(N))
!  masses=1.0d0
  

!  deallocate(masses)
  
!  do i=1,N,1
!     write(*,'("ATOM  ",I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.2,1X)'),&
!          PA(i)%serial,PA(i)%atm_name,&
!          PA(i)%alte_l_i,PA(i)%res_name,&
!          PA(i)%chain_i,mod(PA(i)%res_seq_num,9999),&
!          PA(i)%code_for_insertion_res,&
!          coord(1,i),coord(2,i),coord(3,i)
!  enddo
  
!  deallocate(PA)
!  deallocate(coord)

contains
  subroutine usage(programname)
    implicit none
    character(10) :: programname
    print *,"Usage: ",programname,"PDB x y z t n"
    stop
  end subroutine usage

end program axis_rotation
