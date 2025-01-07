! This program calculates 2D-Free Energy Landscape spanned by Phi- Psi space
! from umbrella sampling simulation data. (2015-07-27)
program wham
  implicit none
  real :: pi, b, tol, minh, mins, maxh, maxs, wid_binh, wid_bins, dim,  C, d, mind, phsi(2)
  real, allocatable :: ef(:), efo(:), ew(:,:,:), bV(:), z(:,:,:), z0(:,:), rouU(:,:), rouB(:,:,:), ewps(:), phi(:), psi(:)
  integer :: n_sim, i, j, k, l, m, h, s, num_binh, num_bins
  integer, allocatable :: n(:)

  real, parameter :: kB = 1.3807488e-23
  real, parameter :: NA = 6.0221413e+23

  character INP_filename*50, Z_filename*50
  character FEL_filename*50

  character(10) :: programname
  integer :: argc

  NAMELIST /whamcontrols/ tol, num_binh, num_bins, n_sim, b, minh, mins, maxh, maxs
  open(17,file='parameters_wham')
  read(17,whamcontrols)
  close(17)

  b = kB * NA / b

  wid_binh = (maxh - minh)/num_binh
  wid_bins = (maxs - mins)/num_bins

  allocate(ef(n_sim))
  allocate(efo(n_sim))
  allocate(bV(n_sim))
  allocate(z(2,n_sim,10))
  allocate(z0(2,n_sim))
  allocate(rouU(num_binh,num_bins))
  allocate(rouB(num_binh,num_bins,n_sim))

  allocate(n(n_sim))

  pi = acos(-1.0d0)

  argc=iargc()
  call getarg(0,programname)
  if (argc < 2) then
     call usage(programname)
  else
     call getarg(1,INP_filename)
     call getarg(2,FEL_filename)
  end if

  ! 1. read (phi, psi) data
  open(19,file=INP_filename)
  do i = 1, n_sim
     read(19,'(A50,1X,2F6.3,1x,F6.3)'), Z_filename, z0(1,i), z0(2,i), bV(i)
     bV(i) = b * bV(i)
     open(21,file=Z_filename)
     n(i) = read_z(21,phi,psi,pi)
!     n(i) = read_z(21,z(1,i,:),z(2,i,:))
     z(1,i,:) = phi
     z(2,i,:) = psi
     close(21)
  end do
  close(19)

  ! 2. calculate exp(-beta x Wi(rj,k))
  do i = 1, n_sim
     allocate(ew(n_sim,n_sim,n(i)))
  end do

  do i = 1, n_sim
     do j = 1, n_sim
        do k = 1, n(j)
           ew(i,j,k) = W_func(z(:,j,l),z0(:,i),bV(i),pi)
           ew(i,j,k) = exp(-1.0*ew(i,j,k))
        end do
     end do
  end do

  ! 3. etaration for ef k
  ef = 1.0e0
  m = 0
  do while (mind > tol)
     efo = ef
     do k = 1, n_sim
        do i = 1, n_sim
           do l = 1, n(i)
              dim = 0.0e0
              do m = 1, n_sim
                 dim = dim + n(j) * ew(i,j,l) * ef(j)
              end do
              ef(k) = ef(k) + n(i) * ew(i,i,l) / dim
           end do
        end do
        ef(k) = 1.0e0 / ef(k)
     end do

     mind=abs(efo(1) - ef(1))
     do i = 2, n_sim
        d = abs(efo(i) - ef(i))
        if (mind > d) then
           mind = d
        end if
     end do

     m = m +1
     print *, m,'th mind=',mind
  end do

  ! 4. make histgram for biased simulations
  rouB = 0.0e0
  do i=1,n_sim
     do j=j,n(i)
        h = ceiling((z(1,i,j) - maxh) / wid_binh)
        s = ceiling((z(2,i,j) - maxs) / wid_bins)
        rouB(h,s,i) = rouB(h,s,i) + 1 / (n(i) * wid_binh * wid_bins)
     end do
  end do

  ! 5. create unbiased histgram
  rouU = 0.0e0
  do h = 1, n_sim
     do s = 1, n_sim
        phsi(1) = minh + wid_binh * (h - 0.5e0)
        phsi(2) = mins + wid_bins * (s - 0.5e0)

        do i = 1, n_sim
           ewps(i) = W_func(phsi,z0(:,i),bV(i),pi)
        end do

        rouU(h,s) = 0.0e0
        do i = 1, n_sim
           dim = 0.0e0
           do j = 1, n_sim
              dim = dim + n(j) * ewps(j) * ef(j)
           end do
           rouU(h,s) = rouU(h,s) + n(i) * rouB(h,s,i) / dim
        end do
     end do
  end do

  C = 0.0e0
  do h=1,num_binh
     do s=1,num_bins
        C = C + rouU(h,s)
     end do
  end do

  rouU(h,s) = rouU(h,s) / C

  ! 6. write (phi, psi, pmf) for output file
  open(23,file=FEL_filename,position='append')
  do h=1,num_binh
     do s=1,num_bins
        write(23,'(F6.3,1x,F6.3,1x,F6.3)'), &
             minh + wid_binh * (h - 0.5e0),mins + wid_bins * (s - 0.5e0), &
             -1.0*log(rouU(h,s))
     end do
  end do
  close(23)

contains

  ! read (phi,psi) data
  integer function read_z(unitnum,phi,psi,pi)
    use reallocate, only: reallocate_real

    implicit none
    integer,intent(in) :: unitnum
    real,intent(in) :: pi
    real,allocatable,intent(inout) :: phi(:), psi(:)

    integer i

    i = 1
    do 
       read(unitnum,'(F8.3,1X,F8.3)',end=999), phi(i),psi(i)
       phi(i) = phi(i) / 180.0 * pi
       psi(i) = psi(i) / 180.0 * pi
       i = i + 1
       call reallocate_real(phi,i)
       call reallocate_real(psi,i)
    end do
999 close(17)

    read_z = i
  end function read_z

  ! calculate restraint function W
  real function W_func(z,z0,bV,pi)
    implicit none
    real,intent(inout) ::  z(2),z0(2),bV,pi

    integer :: i
    real :: dZ,tpi

    tpi = 2.0*pi
    W_func = 0.0e0
    do i = 1,2
       z(i) = amod(z(i),tpi)
       z0(i) = amod(z0(i),tpi)
       dZ = z(i)-z0(i)
       dz = amod(dz,tpi)

       W_func = W_func + 0.5*bV*dZ*dZ
    end do
  end function W_func
    
  ! show usage
  subroutine usage(programname)
    implicit none
    character(10) :: programname

    print *,"Usage: ",programname,"inputfilename  outputfilename"
    stop
  end subroutine usage
end program wham
