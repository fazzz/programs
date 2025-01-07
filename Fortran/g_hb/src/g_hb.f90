! This program calculates the distribution function of the length of 
! the (heavy) atoms of the ligand and solvent
! in NPT ensemble with periodic cubic boundaries

program g_hb
    use xdr, only: xtcfile
    use gettop
    implicit none

    integer  i, j, k, l, m
    
    type(xtcfile) :: xtc
    real(kind=8) pi 
    real(kind=8) dr ! bin distance (nm)
    real(kind=8), allocatable :: posLig(:,:)  ! positions of atoms of ligand molecule
    real(kind=8), allocatable :: posSlv(:,:)  ! positions of atoms of solvent molecule
    real(kind=8), allocatable :: g(:,:),  r(:)   ! rdf g(idx for bin:idx for lig atom)
    real(kind=8) :: xr(3), dxr 
    real(kind=8) :: r_max, box_dim, dv, rho

    integer :: nMolSlv ! # of molecules of ligand and solvent in the system
    integer :: nAtmLig, nAtmSlv ! # of atoms consists of ligand and solvent molecule

    integer :: nHvyAtomLig
    integer :: nHvyAtomSlvMol, nHvyAtomSlvSys

    integer :: nAtmSlvSys ! # of total atoms of solvent molecule in the system

    integer :: IndxIniAtmLig, IndxFinAtmLig ! index for ini/fin atom for solute in the system
    integer :: IndxIniAtmSlv, IndxFinAtmSlv ! index for ini/fin atom for solvent in the system

    integer :: nhist, ng, ig

    ! parameter for atoms
    integer,allocatable,dimension(:) :: nresLig, residLig, cgnrLig
    integer,allocatable,dimension(:) :: nresSlv, residSlv, cgnrSlv
    double precision,allocatable,dimension(:) :: chargLig, massLig
    double precision,allocatable,dimension(:) :: chargSlv, massSlv
    character(5) :: atype
    character(5) :: resLig, resSlv
    character(5), allocatable,dimension(:) :: atomLig, atomSlv
    
    character xtc_filename*50
    character rdf_filename*50
    character itpofLig_filename*50, itpofSolv_filename*50

    character(10) :: argv, programname
    integer :: argc, iargc

    ! 0 reading parameters
    NAMELIST/controls/nMolSlv, nAtmLig, nAtmSlv, IndxIniAtmLig, IndxIniAtmSlv, dr
    open(17,file='parameters_df')
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
       call getarg(3,itpofLig_filename)
       call getarg(4,itpofSolv_filename)
    end if

    allocate(atomLig(nAtmLig), massLig(nAtmLig), nresLig(nAtmLig), residLig(nAtmLig), &
             cgnrLig(nAtmLig), chargLig(nAtmLig))

    open(20,file=itpofLig_filename,status='old')
    call get_top (20, nAtmLig, nresLig, atype, residLig, resLig, atomLig, &
         cgnrLig, chargLig, massLig)
    close(20)

    allocate(atomSlv(nAtmSlv), massSlv(nAtmSlv), nresSlv(nAtmSlv), residSlv(nAtmSlv), &
         cgnrSlv(nAtmSlv), chargSlv(nAtmSlv))

    open(21,file=itpofSolv_filename,status='old')
    call get_top (21, nAtmSlv, nresSlv, atype, residSlv, resSlv, atomSlv, &
         cgnrSlv, chargSlv, massSlv)
    close(21)

    deallocate(nresLig, residLig, cgnrLig, chargLig)
    deallocate(nresSlv, residSlv, cgnrSlv, chargSlv)

    call xtc % init(xtc_filename)

    nAtmSlvSys = nAtmSlv*nMolSlv

    IndxFinAtmLig = IndxIniAtmLig + nAtmLig - 1
    IndxFinAtmSlv = IndxIniAtmSlv + nAtmSlvSys - 1
    
    call xtc % read
    
    ! box information cannot be obtained until at least one read call 
    box_dim = xtc % box(1,1)
    r_max = box_dim / 2d0
    nhist = ceiling(r_max / dr)
    allocate(g(nhist,nAtmLig))
    allocate(r(nhist))
    g = 0d0
    ng = 0
    ! set r-scales as the middle points (Fortran comprehension list) 
    r = [((i - 0.5) * dr, i = 1, nhist)]  

    nHvyAtomLig = 0
    do i = 1, nAtmLig
       if ((index(atomLig(i),'H') /= 1)) then
          nHvyAtomLig = nHvyAtomLig + 1
       end if
    end do

    nHvyAtomSlvMol = 0
    do i = 1, nAtmSlv
       if ((index(atomSlv(i),'H') /= 1)) then
          nHvyAtomSlvMol = nHvyAtomSlvMol + 1
       end if
    end do
    nHvyAtomSlvSys = nHvyAtomSlvMol * nMolSlv

    allocate(posLig(3,nAtmLig))
    allocate(posSlv(3,nAtmSlvSys))

!    write(*,*) nHvyAtomLig, nHvyAtomSlvSys
!    write(*,*) IndxIniAtmLig, IndxFinAtmLig
!    write(*,*) IndxIniAtmSlv, IndxFinAtmSlv
    
    do while ( xtc % STAT == 0 )
       ng = ng + 1
       !       write (*,*), ng

       posLig = xtc % pos(:,IndxIniAtmLig:IndxFinAtmLig)
       posSlv = xtc % pos(:,IndxIniAtmSlv:IndxFinAtmSlv)

!       write(*,*) "YES"
!       write(*,*) "YES"
       
       do i = 1, nAtmLig
          if ((index(atomLig(i),'H') /= 1)) then
!             write(*,*) "YES"
             do j = 1, nMolSlv
                do l = 1, nAtmSlv 
                   if ((index(atomSlv(l),'H') /= 1)) then
!                      write(*,*) "YES"
                      xr = posLig(:,i) - posSlv(:,l+(j-1)*nAtmSlv)
                      !wrap distance (periodic boundary condition)
                      xr = abs(xr - box_dim * nint(xr / box_dim))
                      dxr = sqrt(sum(xr**2))

!                      write(*,*) dxr

                      ig = ceiling(dxr / dr) ! index of g(:)
                      if (ig <= nhist) then
                         if (ig == 0) then
                            ig = 1  !index of g(:) begins from 1
                         end if
                         g(ig,i) = g(ig,i) + 1
                      end if

                   end if
                end do
             end do
          end if
       end do
       call xtc % read
    end do

    ! 5. Close the file
    call xtc % close
       
    ! normalize rdf
    rho = dble(nMolSlv) / box_dim**3  !number density
    do i = 1, nhist
       do j = 1, nAtmLig
          dv = (4d0 / 3d0) * pi * (i**3 - (i-1)**3) * dr**3
          g(i,j) = g(i,j) / (ng  * dv  * rho)
       end do
    end do

    ! output results
    open(19,file=rdf_filename,position='append')
    do i = 1, nhist
       
       write(19,'(f20.8,3x)',advance='no') r(i)
       do j = 1, nHvyAtomLig
          write(19,'(f20.8,1x)',advance='no') g(i,j)
       end do
       write(19,*)
       
       !       write(19,'(2f20.8)') r(i), g(i,1)
    end do
    close(19)

    deallocate(atomLig)
    deallocate(atomSlv)
    
    deallocate(posLig)
    deallocate(posSlv)

    deallocate(g)
    deallocate(r)

    deallocate(massLig)
    deallocate(massSlv)
    
  contains
  subroutine usage(programname)
    implicit none
    character(10),intent(in) :: programname

    print *,"Usage: ",programname," xtcfilename(in) rdffilename(out) itpligfilename(in) itpslvfilename(in)"
    stop
  end subroutine usage
end program g_hb


