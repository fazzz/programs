program read_list_Boltzmann_Medal_Winners
  implicit none
  integer Num
  double precision charge, LJ1, LJ2
  character Atom, filename*50

  character(50) :: argv
  integer :: argc

  argc=iargc()
  if (argc < 1) then
     call getarg(0,argv)
     call usage(argv)
  else
     call getarg(1,filename)
  end if
  
  open(17,file=filename,status='old')
  do 
     read(17,'(I4,2X,A20,1X,A50)',end=999),year,name,belonging
     write(*,'(I4,2X,A20,1X,A50)'),year,name,belonging
  end do
  999 close(17)

end program read_list_Boltzmann_Medal_Winners
