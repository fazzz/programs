! This program perform bestfit two pdb strcturues
! based on all atoms
! You can specify the atom range to be fited.

program mix_pdb
  implicit none

  integer i,j,N ! num of atoms

  double precision :: coord1(3), coord2(3), w
  
  character(30) line
  character(60) pdbfilename1, pdbfilename2

  character(30) :: argv, programname, d
  integer :: argc, iargc
  
  argc=iargc()
  call getarg(0,programname)
  if (argc < 4) then
     call usage(programname)
  else
     call getarg(1,pdbfilename1)
     call getarg(2,pdbfilename2)
     call getarg(3,d)
     read(d,*),w
     call getarg(4,d)
     read(d,*),N
  end if

  open(16,file=pdbfilename1,status='old')
  open(17,file=pdbfilename2,status='old')
  do i=1,N,1
     read(16,'(A30,3F8.3)') line,coord1(1),coord1(2),coord1(3)
     read(17,'(A30,3F8.3)') line,coord2(1),coord2(2),coord2(3)

     coord1(:) = w * coord1(:) + (1.0-w) * coord2(:) 

     write(*,'(A30,3F8.2,1X)') line,coord1(1),coord1(2),coord1(3)
     
  end do
  close(16)
  close(17)
  
contains
  subroutine usage(programname)
    implicit none
    character(10) :: programname
    print *,"Usage: ",programname,"PDB1 PDB2 w(weight for PDB1) N(num of atoms)"
    stop
  end subroutine usage

end program mix_pdb
