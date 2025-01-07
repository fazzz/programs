! This program calculates the com of mol A - com of mol B rdf
! in NVT ensemble with periodic cubic boundaries
program rdf
    use xdr, only: xtcfile
    implicit none
    type(xtcfile) :: xtc
    real(kind=8) pi 
    real(kind=8) dr !nm
    real(kind=8), allocatable :: posA(:,:), posB(:,:) ! positions of mol A/B
    real(kind=8), allocatable :: comA(:,:), comB(:,:) ! positions of com of mol A/B
    real(kind=8), allocatable :: massA(:), massB(:)
    real(kind=8), allocatable :: g(:), r(:) ! rdf
    real(kind=8) :: xr(3), dxr  !dxr: distance between two coms
    real(kind=8) :: r_max, box_dim, dv, rho
    integer :: nMolA, nMolB, nAtomMolA, nAtomMolB
    integer :: nhist, ng, ig, npos, nposA, nposB, i, j

    character xtc_filename*50
    character rdf_filename*50
    character itpofA_filename*50, itpofB_filename*50

    character(10) :: argv, programname
    integer :: argc, iargc

    ! 0 reading parameters
    NAMELIST/controls/nMolA, nMolB, nAtomMolA, nAtomMolB, dr
    open(17,file='parameters_rdf')
    read(17,nml=controls)
    close(17)

    pi = acos(-1.0)

    argc=iargc()
    call getarg(0,programname)
    if (argc < 4) then
       call usage(programname)
    else
       call getarg(1,xtc_filename)
       call getarg(2,rdf_filename)
       call getarg(3,itpofA_filename)
       call getarg(4,itpofB_filename)
    end if

    call xtc % init(xtc_filename)

    write(*,*),"nMolA ","nMolB"
    write(*,*),nMolA,nMolB
    write(*,*),nAtomMolA,nAtomMolB

    allocate(massA(nAtomMolA))
    allocate(massB(nAtomMolB))

    ! 1 reading mol data
    open(18,file=itpofA_filename,status='old')
!    write(*,*),nAtomMolA
    call get_top (18, nAtomMolA, massA)
    close(18)
!    write(*,*),massA

    open(19,file=itpofB_filename,status='old')
    call get_top (19, nAtomMolB, massB)
    close(19)

!    write(*,*),nMolA,nMolB
!    write(*,*),nAtomMolA,nAtomMolB

    ! number of total atoms 
    npos = nAtomMolA * nMolA + nAtomMolB * nMolB
    nposA = nAtomMolA * nMolA
    nposB = nAtomMolB * nMolB
    allocate(posA(3, nposA))
    allocate(posB(3, nposB))
    allocate(comA(3, nMolA))
    allocate(comB(3, nMolB))

    call xtc % read
    
    ! box information cannot be obtained until at least one read call 
    box_dim = xtc % box(1,1)
    r_max = box_dim / 2d0
    nhist = ceiling(r_max / dr)
    allocate(g(nhist))
    allocate(r(nhist))
    g = 0d0
    ng = 0
    r = [((i - 0.5) * dr, i = 1, nhist)]  ! set r-scales as the middle points (Fortran comprehension list) 

    do while ( xtc % STAT == 0 )
        ng = ng + 1
        write (*,*), ng
        posA = xtc % pos(:, 1:nposA)       !get the position of atoms of mol A
        posB = xtc % pos(:, nposA+1:npos)  !get the position of atoms of mol B
!        call pos2com(posA,nMolA,nAtomMolA,massA,comA)
!        call pos2com(posB,nMolB,nAtomMolB,massB,comB)

        comA = posA(:,1:nposA:nAtomMolA)
        comB = posB(:,1:nposB:nAtomMolB)

        do i = 1, nMolA
          do j = 1, nMolB
             xr = comA(:,i) - comB(:,j)             
             xr = abs(xr - box_dim * nint(xr / box_dim)) !wrap distance (periodic boundary condition)
             dxr = sqrt(sum(xr**2))
             ig = ceiling(dxr / dr)
             if (ig <= nhist) then
                if (ig == 0) then
                   ig = 1  !index of g(:) begins from 1
                end if
                g(ig) = g(ig) + 1
             end if
          end do
       end do
       call xtc % read
    end do

    ! 5. Close the file
    call xtc % close

    ! normalize rdf
    rho = dble(npos) / box_dim**3  !number density
    do i = 1, nhist
       dv = (4d0 / 3d0) * pi * (i**3 - (i-1)**3) * dr**3
       g(i) = g(i) / (ng * (nMolA + nMolB) * dv * rho)
    end do

    ! output results
    open(19,file=rdf_filename,position='append')
    do i = 1, nhist
      write(19,'(2f20.8)') r(i), g(i)
   end do
   close(19)

   deallocate(posA)
   deallocate(posB)
   deallocate(comA)
   deallocate(comB)

   deallocate(g)
   deallocate(r)

contains
  subroutine usage(programname)
    implicit none
    character(10),intent(in) :: programname

    print *,"Usage: ",programname," inputfilename outputfilename itpofAfilename itpofBfilename"
    stop
  end subroutine usage

  subroutine get_top( unitnum, NA, mass )
    implicit none

    integer,intent(in) :: unitnum
    integer,intent(in) :: NA
    real(kind=8),allocatable,dimension(:),intent(inout) :: mass

    integer :: nr, resid, cgnr
    real(kind=8) :: charg
    character(5) :: atype, res, atom

    integer i

    do i=1,NA
       read(unitnum,'(I6,A5,I6,1x,A5,1x,A5,I5,4x,F9.6,5x,F8.5)'),&
            nr, atype, resid, res, atom, cgnr, charg, mass(i)
    enddo

  end subroutine get_top

  subroutine pos2com( pos, nMol, nAtomMol, mass, com)
    implicit none

    real(kind=8), allocatable,intent(in) :: pos(:,:), mass(:)
    real(kind=8), allocatable,intent(inout) :: com(:,:)
    integer,intent(in) :: nMol, nAtomMol

    real(kind=8), allocatable :: position(:,:), compos(:)
    integer i

    allocate(position(3,nAtomMol))
    allocate(compos(3))

    do i = 1, nMol
       position = pos(:,1+(i-1)*nAtomMol:i*nAtomMol)
       call CenterOfMass(position,nAtomMol,mass,compos)
       com(:,i) = compos(:)
    end do

    deallocate(position)
    deallocate(compos)

  end subroutine pos2com

  subroutine CenterOfMass(pos, nAtom, mass, com)
    implicit none

    real(kind=8), allocatable, intent(in) :: pos(:,:), mass(:)
    integer, intent(in) :: nAtom
    real(kind=8), allocatable, intent(inout) :: com(:)

    real(kind=8) summass
    integer i

    summass = 0.0
    com  = 0.0

    do i=1,nAtom
       summass = summass + mass(i)
       com = com + pos(:,i)*mass(i)
    end do

    com = com / summass

  end subroutine CenterOfMass

end program rdf
