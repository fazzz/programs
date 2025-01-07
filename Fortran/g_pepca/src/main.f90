! This program performs Potential Energy Principal Component Analysis (PEPCA)
! for Energy Trajectory.

program g_pepca
  use pepca

  implicit none

  integer :: i, j, k, l, m
  integer :: nframe, nene ! # of frame, # of ene
  
  double precision, allocatable, dimension(:)     :: ev       ! eigen value
  double precision, allocatable, dimension(:,:)   :: varm     ! variance-covariance matrix
  double precision, allocatable, dimension(:,:) :: ene(:,:), B(:,:)
  double precision, allocatable, dimension(:,:) :: pc
  double precision v(5)

  double precision sumev,cont

  character ene_filename*50
  character out_filename*50
  character eig_filename*50

  character(10) :: argv, programname
  integer :: argc, iargc

  NAMELIST/mdinfo/nene
  open(16,file='parameters_pepca')
  read(16,mdinfo)
  close(16)
  
  argc=iargc()
  call getarg(0,programname)
  if (argc < 3) then
     call usage(programname)
  else
     call getarg(1,ene_filename)
     call getarg(2,out_filename)
     call getarg(3,eig_filename)
  end if

  allocate(ene(nene,1))

  i = 1 ! index for frame
  j = 0 ! index for dimension
  open(21,file=ene_filename,status='old')
  do 
     ! input file format is as below.
     read(21,'(F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4)',end=999),v(1),v(2),v(3),v(4),v(5)
     do k = 1,5,1
        j = j + 1
        if ( j > nene ) then
           allocate(B(nene,i))
           B = ene
           i = i + 1
           j = 1
           deallocate(ene)
           allocate(ene(nene,i))
           ene(1:nene,1:i-1) = B(1:nene,1:i-1)
           deallocate(B)
           ene(j,i) = v(k)
!           write(*,'(F8.3)') ene(j,i)
        else
           ene(j,i) = v(k)
!           write(*,'(F8.3)') ene(j,i)
        end if
     end do
  end do
999 close(21)

  nframe = i - 1

!  do i = 1, ndih
!     write(*,'(F8.3,1x)',advance='NO') dih(i,1)
!  end do
     
!  write(*,*)

!  do i = 1, ndih
!     write(*,'(F8.3,1x)',advance='NO') dih(i,nframe)
!  end do

!  write(*,*)

!  stop

!  stop

  ! calculate variance covariance matrix
  call enevc(nene, ene, varm,  nframe)

!  write(*,*) "YES"
  
  allocate(ev(nene))

!  write(*,*) "YES"
  
  ! calculate eigen value and eigen vector
  call pepcaeigvec(nene, varm, ev)

!  write(*,*) "YES"
  
  sumev = sum(ev)
  
  cont = 0.0D0
  do i = 0, nene-1
     cont = cont + ev(nene-i)/sumev
     write (*, '(I4,1X,E12.7,1X,F10.6)') , i, ev(nene-i), cont*100
  end do

  open(19,file=eig_filename,status='new')
  do i = 1, nene
     write (19, '(I4,1X)',advance='no'), i
!     do j =1, npca
!        write (19, '(F10.6,1X)',advance='no') , varm(j,i)
!     end do
     write (19, '(F10.6,1X,F10.6,1X)',advance='no'), varm(1,i), varm(2,i)
     write(19,'(/)',advance='no')
  end do
  close(19)
  
  allocate(pc(2,nframe))

  pc = 0.0d0
  
!  write(*,*) "YES"
  
  open(20,file=out_filename,status='new')
  do i = 1, nframe
     do j = 1, nene
        l = (j-1)*2+k
        pc(1,i) = pc(1,i) + varm(1,j) * ene(j,i)
        pc(2,i) = pc(2,i) + varm(2,j) * ene(j,i)
     end do
     write (20, '(F12.4,1X,F12.4,1X)', advance='no'), pc(1,i), pc(2,i)
     write(20,'(/)',advance='no')
  end do
  close(20)
  
  deallocate(ev, varm)
  deallocate(pc)
  deallocate(ene)
  
contains

  ! Show help mesage
  subroutine usage(programname)
    implicit none
    character(10) :: programname
    
    print *,"Usage: ",programname,"ene_filename out_filename eig_filename"
    stop
  end subroutine usage
  
end program g_pepca
