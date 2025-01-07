! This program conducts Sparse Structure Learning-Principle Component Analysis (SSL-PCA)
! for xtc file (GROMACS trajectory)

program g_sslpca ! iA <- x
  use xdr, only: xtcfile
  use gettop
  use bestfit
  use standardlize
  use graphical_Lasso
  use pca

  implicit none

  integer i,j,k,l
  integer :: nframe, natom, ndim ! # of frame, # of atom, # of dimensions
  integer :: pc1, pc2

  type(xtcfile) :: xtc ! trajectory (GROMACS format)

  ! parameters for molecule
  character(5) :: atomname, resname
  integer,allocatable,dimension(:) :: nres, resid, cgnr
  double precision,allocatable,dimension(:) :: charg, mass
  character(5) :: atype, res, atom

  double precision,allocatable,dimension(:,:) :: com ! center of mass trajectory

  double precision, allocatable, dimension(:)     :: ev       ! eigen value
  double precision, allocatable, dimension(:,:)   :: varm     ! variance-covariance matrix
  double precision, allocatable, dimension(:,:,:) :: coord, bfcoord
  double precision, allocatable, dimension(:,:)   :: pc

  double precision sumev,cont

  double precision,allocatable,dimension(:,:) :: x     ! (n x m), n times, m dimension
  double precision,allocatable,dimension(:)   :: avx   ! (m), average of data
  double precision,allocatable,dimension(:,:) :: A, iA ! (m x m), presicion matrix
  double precision,allocatable,dimension(:,:) :: S     ! (m x m), variance-covariance matrix
  double precision,allocatable,dimension(:,:) :: B     ! dummy matrix
  double precision rou                                 ! degree of sparsity
  double precision v(5)                                ! dummy values

  character top_filename*50
  character xtc_filename*50
  character out_filename*50
  character eig_filename*50
  character smt_filename*50

  character(10) :: argv, programname
  integer :: argc, iargc

  NAMELIST/sslpcacontrols/natom, nframe, rou, pc1, pc2
  open(16,file='parameters_sslpca')
  read(16,sslpcacontrols)
  close(16)

  argc=iargc()
  call getarg(0,programname)
  if (argc < 5) then
     call usage(programname)
  else
     call getarg(1,top_filename)
     call getarg(2,xtc_filename)
     call getarg(3,out_filename)
     call getarg(4,eig_filename)
     call getarg(5,smt_filename)
  end if

  allocate(mass(natom), nres(natom), resid(natom), cgnr(natom), charg(natom))
  
  ! get topological data
  open(19,file=top_filename,status='old')
  call get_top (19, natom, nres, atype, resid, res, atom, cgnr, charg, mass)
  close(19)

  deallocate(cgnr, resid, charg)

  allocate(coord(3,natom,nframe), bfcoord(3,natom,nframe))
  
  ! open xtc file
  call xtc % init(xtc_filename)

  call xtc % read

  i = 1
  do while (xtc % STAT == 0)
     ! read xtc file
     coord(:,:,i) = xtc % pos(:,:)
     i = i + 1
     call xtc % read
  end do
  
  call xtc % close

  allocate(com(3, nframe))

  ! remove the translation and rotation
  call bftrj(natom, nframe, coord, mass, com, bfcoord)

  deallocate(com)
  
  deallocate(coord)

  ndim = natom * 3
  
  allocate(x(ndim,nframe))

  do i = 1, nframe
     x(1:natom,i)           = bfcoord(1,1:natom,i)
     x(natom+1:natom*2,i)   = bfcoord(2,1:natom,i)
     x(natom*2+1:natom*3,i) = bfcoord(3,1:natom,i)
  enddo
  
  allocate(avx(ndim))
  allocate(S(ndim,ndim))
  allocate(A(ndim,ndim), iA(ndim,ndim))

  call ave_timeseries(x,avx,nframe,ndim)   ! avx <- x
  call vcv_timeseries(x,nframe,ndim,avx,S) ! S <- x

  open(19,file=smt_filename,status='new')
  write(19,'("The variance-covarucance matrix = ")')
  do i=1,ndim,1
     do j=1,ndim,1
        write(19,'(F8.3,1X)',advance='no'), S(i,j)
     end do
     write(19,'(/)',advance='no')
  end do
    
  deallocate(avx)

  call block_coordinate_descent_method(A,iA,S,rou,ndim) ! A, iA, <- S, rou, m

  write (19,'("The presicion matrix =")')
  do i=1,ndim
     do j=1,ndim
        write (19,'(F8.3,1X)',advance='no'),A(i,j)
     end do
     write (19,'(/)',advance='no')
  end do
  write (19,'("The inverse matrix of the presicion matrix =")')
  do i=1,ndim
     do j=1,ndim
        write (19,'(F8.3,1X)',advance='no'),iA(i,j)
     end do
     write (19,'(/)',advance='no')
  end do

  write (19,'("The Sparse structure of the time series is (the non-zero value of the matrix)")')
  do i=1,ndim,1
     do j=i+1,ndim,1
        if ( A(i,j) /= 0.0d0 ) then
           write(19,'(I3,1X,"-",1X,I3)'), i,j
        end if
     end do
  end do
  close(19)
  
  deallocate(S)
  deallocate(A)

  allocate(ev(ndim))
  
  ! calculate eigen value and eigen vector
  call eigvec(natom, iA, ev)

  sumev = sum(ev)
  
  cont = 0.0D0
  do i = 0, 3*natom-1
     cont = cont + ev(3*natom-i)/sumev
     write (*, '(I4,1X,E20.10,1X,F8.4)') , i, ev(3*natom-i), cont*100
  end do

  open(19,file=eig_filename,status='new')
  do i = 1, 3*natom
     write (19, '(I4,1X)',advance='no'), i
     write (19, '(F10.6,1X,F10.6,1X)',advance='no'), iA(pc1,i), iA(pc2,i)
     write(19,'(/)',advance='no')
  end do
  close(19)
  
  allocate(pc(2,nframe))

  pc = 0.0d0

  open(20,file=out_filename,status='new')
  do i = 1, nframe
     do j = 1, natom
        do k = 1, 3
           l = natom*(k-1)+j
           pc(1,i) = pc(1,i) + iA(pc1,l) * bfcoord(k,j,i)
           pc(2,i) = pc(2,i) + iA(pc2,l) * bfcoord(k,j,i)
        end do
     end do
     write (20, '(F8.4,1X,F8.4,1X)', advance='no'), pc(1,i), pc(2,i)
     write(20,'(/)',advance='no')
  end do
  close(20)
  
  deallocate(bfcoord)
  deallocate(ev, iA)
  deallocate(pc)

contains

  subroutine usage(programname)
    implicit none
    character(10) :: programname

    print *,"Usage: ",programname," topfilename xtcfilename outfilename eigfilename smtfilename"
    stop
  end subroutine usage





end program g_sslpca
