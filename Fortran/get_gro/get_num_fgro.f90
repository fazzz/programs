program get_num_fgro
  use gro

  implicit none

  integer N ! num of atoms

  real boxsize_x, boxsize_y, boxsize_z ! length of box ( nm )

  type (gro_atom), allocatable, dimension(:) :: GA

  character(50) gro_filename

  character(10) :: argv, programname
  integer :: argc, iargc

  argc=iargc()
  call getarg(0,programname)
  if (argc < 1) then
     call usage(programname)
  else
     call getarg(1,gro_filename)
  end if

  open(17,file=gro_filename,status='old')
  call read_GRO (17, GA, N, boxsize_x, boxsize_y, boxsize_z)
  close(17)

  write (*,*), N

  deallocate(GA)

contains
  subroutine usage(programname)
    implicit none
    character(10) :: programname

    print *,"Usage: ",programname," gro_file"
  end subroutine usage

end program get_num_fgro
