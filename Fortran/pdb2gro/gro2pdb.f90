program gro2pdb
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
  if (argc < 2) then
     call usage(programname)
  else
     call getarg(1,gro_filename)
     call getarg(2,pdb_filename)
  end if

  open(17,file=gro_filename,status='old')
  call read_GRO (17, GA, N, boxsize_x, boxsize_y, boxsize_z)
  close(17)

  print *,"yes"

  allocate(PA(N))

  call cp_gro2pdb(GA, PA, N)

  print *,"yes"

  open(19,file=pdb_filename,position='append')
  call write_PDB (19, PA, N)
  close(19)

  print *,"yes"

  deallocate(PA)
  deallocate(GA)

contains
  subroutine cp_gro2pdb(GA, PA, N)
    type (GRO_atom),allocatable,dimension(:),intent(in) :: GA
    type (PDB_atom),allocatable,dimension(:),intent(inout) :: PA
    integer,intent(in) :: N

    integer i

    do i=1,N
!       PA(i)%nameRecord="ATOM  "
       PA(i)%serial=GA(i)%atm_serial
       PA(i)%atm_name=adjustl(GA(i)%atm_name)
!       PA(i)%alte_l_i=" "
       PA(i)%res_name=GA(i)%res_name
!       PA(i)%chain_i=" "

       PA(i)%res_seq_num=GA(i)%res_serial
!       PA(i)%code_for_insertion_res=" "
       PA(i)%x=GA(i)%x * 10.0e0
       PA(i)%y=GA(i)%y * 10.0e0
       PA(i)%z=GA(i)%z * 10.0e0
!       PA(i)%occupancy=1.0e0
!       PA(i)%temperature_factor=10.0e0
!       PA(i)%seq_i=" "
!       PA(i)%element=" "
    end do
  end subroutine cp_gro2pdb

  subroutine usage(programname)
    implicit none
    character(10) :: programname

    print *,"Usage: ",programname," gro_file PDB_file"
    stop
  end subroutine usage

end program gro2pdb
