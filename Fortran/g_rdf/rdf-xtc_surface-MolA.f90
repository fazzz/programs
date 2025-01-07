! This program calculates the atom on surface of protein - atom a of MolA rdf
! in NPT ensemble with periodic cubic boundaries
program rdf
    use xdr, only: xtcfile
    implicit none
    type(xtcfile) :: xtc
    real(kind=8) pi 
    real(kind=8) dr !nm
    real(kind=8), allocatable :: posA(:,:), posPRO(:,:) ! positions of atom a/b of mol A/B
    real(kind=8), allocatable :: g(:), r(:) ! rdf
    real(kind=8) :: xr(3), dxr  !dxr: distance between two atoms
    real(kind=8) :: r_max, box_dim, dv, rho
    integer :: nMolA,  nAtomMolA, nAtomPRO
    integer :: nthbmolA, nthbatomamolA, nthemolA, nthbmolA_Sys
    integer :: nthbmolPRO, nthbatomamolPRO, nthemolPRO, nthbmolPRO_Sys
    integer :: nhist, ng, ig, i, j

    character xtc_filename*50
    character rdf_filename*50
    character itpofA_filename*50, itpofB_filename*50, itpofProt_filename*50

    character(10) :: argv, programname
    integer :: argc, iargc

    ! 0 reading parameters
    NAMELIST/controls/nMolA, nAtomMolA, nAtomPRO, nthbmolA, nthbatomamolA, nthbmolPRO, nthbatomamolPRO,  dr
    open(17,file='parameters_rdf')
    read(17,nml=controls)
    close(17)

    pi = acos(-1.0)

    argc=iargc()
    call getarg(0,programname)
    if (argc < 2) then
       call usage(programname)
    else
       call getarg(1,xtc_filename)
       call getarg(2,rdf_filename)
    end if

    call xtc % init(xtc_filename)

    ! number of total atoms 
    nthbmolA_Sys = nthbmolA + nthbatomamolA-1
    nthbmolB_Sys = nthbmolB + nthbatomamolB-1
    nthemolA = nthbmolA + nAtomMolA * nMolA - 1
    nthemolB = nthbmolB + nAtomMolB * nMolB - 1
    allocate(posA(3, nMolA))
    allocate(posB(3, nMolB))

    call xtc % read
    
    ! box information cannot be obtained until at least one read call 
    box_dim = xtc % box(1,1)
    r_max = box_dim / 2d0
    nhist = ceiling(r_max / dr)
    allocate(g(nhist))
    allocate(r(nhist))
    g = 0d0
    rho = 0d0
    ng = 0
    ! set r-scales as the middle points (Fortran comprehension list) 
    r = [((i - 0.5) * dr, i = 1, nhist)]  

    write (*,*) nthbmolA,nthbatomamolA,nthemolA,nAtomMolA
    write (*,*) nthbmolB,nthbatomamolB,nthemolB,nAtomMolB

    do while ( xtc % STAT == 0 )
        ng = ng + 1
        write (*,*), ng

        box_dim = xtc % box(1,1)
        rho = rho + dble(nMolA) / box_dim**3  !number density

        do i = 1, nMolA
           !get the position of atoms of mol A
           posA = xtc % pos(:,nthbmolA_Sys:nthemolA:nAtomMolA)
          do j = 1, nMolB
             !get the position of atoms of mol B
             posB = xtc % pos(:,nthbmolB_Sys:nthemolB:nAtomMolB)
             xr = posA(:,i) - posB(:,j)             
             !wrap distance (periodic boundary condition)
             xr = abs(xr - box_dim * nint(xr / box_dim))
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
    do i = 1, nhist
       dv = 4d0 * pi * (i*dr) **2 * dr
       g(i) = g(i) / (ng * nMolB * dv ) * rho
    end do

    ! output results
    open(19,file=rdf_filename,position='append')
    do i = 1, nhist
      write(19,'(2f20.8)') r(i), g(i)
   end do
   close(19)

   deallocate(posA)
   deallocate(posB)

   deallocate(g)
   deallocate(r)

contains
  subroutine usage(programname)
    implicit none
    character(10),intent(in) :: programname

    print *,"Usage: ",programname," inputfilename outputfilename"
    stop
  end subroutine usage

end program rdf
