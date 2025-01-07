! This program makes the list of atoms
! which belong to the solvent accesible surface area (SASA).
program sasa
    use gro
    implicit none
    type (gro_atom), allocatable, dimension(:) :: GA
    real(kind=8) pi 
    real(kind=8) l ! distance bet. atom i and j
    real(kind=8), allocatable :: pos(:,:)   ! positions of atom
    real(kind=8), allocatable :: r_vdw(:), r_i(:) ! vdW radius of atoms, i-th sphere radius
    real(kind=8) :: r_wat,                   ! vdW radius of wat.
    real(kind=8) :: phi, theta, dphi, dtheta, x, y, z
    integer, allocatable :: list_sasa(:)    ! list of atoms belongs to SASA
    integer, allocatable :: list_near(:,:)  ! list of atoms which is naer atom i
    integer :: nAtom                        ! # of atoms
    integer :: i, j

    character gro_filename*50

    character(10) :: argv, programname
    integer :: argc, iargc

    ! 0 reading parameters
    NAMELIST/controls/dphi, dtheta, r_wat
    open(17,file='parameters_sasa')
    read(17,nml=controls)
    close(17)

    pi = acos(-1.0)

    argc=iargc()
    call getarg(0,programname)
    if (argc < 2) then
       call usage(programname)
    else
       call getarg(1,gro_filename)
    end if

    call read_GRO (17, GA, nAtom, box_dim)

    allocate(pos(3, nMolA))
    allocate(r(nhist))

    ! make list of atom which is near atom i
    do i = 1, nMolA
       k = 1
       do j = 1, nMolA
          if ( i .neq j ) then
             xr = pos(:,i) - pos(:,j)
             !wrap distance (periodic boundary condition)
             xr = abs(xr - box_dim * nint(xr / box_dim))
             if ( xr .gt. ) then
                k = k + 1
                call reallocate(list_near,k)
                list_near(i,k)=j
             end if
          end if
       end do
    end do

    do i = 1, nAtom
       NUM = list_near(i,1)
       do j =1, NUM
          do t=1, numTheta
             do p=1, numPsi
                
             end do
          end do
       end do
    end do

contains
  subroutine usage(programname)
    implicit none
    character(10),intent(in) :: programname

    print *,"Usage: ",programname," inputfilename"
    stop
  end subroutine usage

end program sasa
