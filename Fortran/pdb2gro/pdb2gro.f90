program pdb2gro
  use pdb
  use gro

  implicit none

  integer N ! num of atoms

  integer i

  real boxsize_x, boxsize_y, boxsize_z ! length of box ( nm )

  type (PDB_atom), allocatable, dimension(:) :: PA
  type (gro_atom), allocatable, dimension(:) :: GA

  character(50) pdb_filename
  character(50) gro_filename

  character(10) :: argv, programname
  integer :: argc, iargc

  argc=iargc()
  call getarg(0,programname)
  if (argc < 5) then
     call usage(programname)
  else
     call getarg(1,pdb_filename)
     call getarg(2,gro_filename)
     call getarg(3,argv)
     read(argv,*) boxsize_x
     call getarg(4,argv)
     read(argv,*) boxsize_y
     call getarg(5,argv)
     read(argv,*) boxsize_z
  end if

  open(17,file=pdb_filename,status='old')
  call read_PDB (17, PA, N)
  close(17)

!  print *,"yes"

  allocate(GA(N))

  call cp_pdb2gro(PA, GA, N)

!  print *,"yes"

  open(19,file=gro_filename,position='append')
  call write_GRO (19, GA, N, boxsize_x, boxsize_y, boxsize_z)
  close(19)

!  print *,"yes"

  deallocate(PA)
  deallocate(GA)

contains
  subroutine cp_pdb2gro(PA, GA, N)
    type (PDB_atom),allocatable,dimension(:),intent(in) :: PA
    type (GRO_atom),allocatable,dimension(:),intent(inout) :: GA
    integer,intent(in) :: N

    integer i, j, l

    character(4) TEMP

    do i=1,N
       GA(i)%res_serial=PA(i)%res_seq_num
       GA(i)%atm_serial=i !PA(i)%serial

       GA(i)%res_name=PA(i)%res_name
       GA(i)%atm_name=PA(i)%atm_name

       GA(i)%x=PA(i)%x / 10.0e0
       GA(i)%y=PA(i)%y / 10.0e0
       GA(i)%z=PA(i)%z / 10.0e0
    end do
  end subroutine cp_pdb2gro

  subroutine usage(programname)
    implicit none
    character(10) :: programname

    print *,"Usage: ",programname," PDB_file gro_file boxsize_x boxsize_y boxsize_z"
    stop
  end subroutine usage

end program pdb2gro
