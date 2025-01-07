! This program calculates the distribution function of the length of 
! the most neiboring (heavy) atoms of the solute (protein) and solvent
! in NPT ensemble with periodic cubic boundaries

program g_df
    use xdr, only: xtcfile
    use gettop
    implicit none

    type(xtcfile) :: xtc
    real(kind=8) pi 
    real(kind=8) dr ! bin distance (nm)
    real(kind=8), allocatable :: posSlu(:,:)  ! positions of atoms of solute molecule
    real(kind=8), allocatable :: posSlv(:,:)  ! positions of atoms of solvent molecule
    real(kind=8), allocatable :: g(:), gminp(:), ghat(:), gpcom(:), r(:)   ! rdf
    real(kind=8) :: xr(3), dxr, min_dxr ! min_dxr: dist. bet. two atoms of solute and solvent
    real(kind=8) :: r_max, box_dim, dv, rho, rho2

    integer :: nMolSlu, nMolSlv ! # of molecules of solute and solvent in the system
    integer :: nAtmSlu, nAtmSlv ! # of atoms consists of solute and solvent molecule

    integer :: nAtmSluSys ! # of total atoms of solute molecule in the system
    integer :: nAtmSlvSys ! # of total atoms of solvent molecule in the system
    
    integer :: IndxIniAtmSlu, IndxFinAtmSlu ! index for ini/fin atom for solute in the system
    integer :: IndxIniAtmSlv, IndxFinAtmSlv ! index for ini/fin atom for solvent in the system

    integer :: nhist, ng, ig, i, j, k, l, m

    integer,allocatable,dimension(:) :: nresSlu, residSlu, cgnrSlu
    integer,allocatable,dimension(:) :: nresSlv, residSlv, cgnrSlv
    double precision,allocatable,dimension(:) :: chargSlu, massSlu
    double precision,allocatable,dimension(:) :: chargSlv, massSlv
    character(5) :: atype
    character(5) :: resSlu, resSlv
    character(5), allocatable,dimension(:) :: atomSlu, atomSlv

    double precision :: com(3)
    
    integer :: nHvyAtomSluMol, nHvyAtomSluSys
    integer :: nHvyAtomSlvMol, nHvyAtomSlvSys
    
    character xtc_filename*50
    character rdf_filename*50
    character itpofSolute_filename*50, itpofSolvent_filename*50

    character(10) :: argv, programname
    integer :: argc, iargc

    ! 0 reading parameters
    NAMELIST/controls/nMolSlu, nMolSlv, nAtmSlu, nAtmSlv, IndxIniAtmSlu, IndxIniAtmSlv, dr
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
       call getarg(3,itpofSolute_filename)
       call getarg(4,itpofSolvent_filename)
    end if

    allocate(atomSlu(nAtmSlu))
    allocate(massSlu(nAtmSlu), nresSlu(nAtmSlu), residSlu(nAtmSlu), &
         cgnrSlu(nAtmSlu), chargSlu(nAtmSlu))

    open(20,file=itpofSolute_filename,status='old')
    call get_top (20, nAtmSlu, nresSlu, atype, residSLu, resSlu, atomSlu, &
         cgnrSlu, chargSlu, massSlu)
    close(20)

    allocate(atomSlv(nAtmSlv))
    allocate(massSlv(nAtmSlv), nresSlv(nAtmSlv), residSlv(nAtmSlv), &
         cgnrSlv(nAtmSlv), chargSlv(nAtmSlv))

    open(21,file=itpofSolvent_filename,status='old')
    call get_top (21, nAtmSlv, nresSlv, atype, residSlv, resSlv, atomSlv, &
         cgnrSlv, chargSlv, massSlv)
    close(21)

    deallocate(nresSlu, residSlu, cgnrSlu, chargSlu)
    deallocate(nresSlv, residSlv, cgnrSlv, chargSlv)

    call xtc % init(xtc_filename)

    nAtmSluSys = nAtmSlu*nMolSlu
    nAtmSlvSys = nAtmSlv*nMolSlv

    IndxFinAtmSlu = IndxIniAtmSlu + nAtmSluSys - 1
    IndxFinAtmSlv = IndxIniAtmSlv + nAtmSlvSys - 1
    
    allocate(posSlu(3,nAtmSluSys))
    allocate(posSlv(3,nAtmSlvSys))

    call xtc % read
    
    ! box information cannot be obtained until at least one read call 
    box_dim = xtc % box(1,1)
    r_max = box_dim / 2d0
    nhist = ceiling(r_max / dr)
    allocate(g(nhist))
    allocate(gminp(nhist))
    allocate(ghat(nhist))
    allocate(gpcom(nhist))
    allocate(r(nhist))
    g = 0d0
    gminp = 0d0
    ghat = 0d0
    gpcom = 0d0
    ng = 0
    ! set r-scales as the middle points (Fortran comprehension list) 
    r = [((i - 0.5) * dr, i = 1, nhist)]  

    nHvyAtomSlvMol = 0
    do i = 1, nAtmSlv 
       if ((index(atomSlv(i),'H') /= 1)) then
          nHvyAtomSlvMol = nHvyAtomSlvMol + 1
       end if
    end do
    nHvyAtomSlvSys = nHvyAtomSlvMol * nMolSlv

    nHvyAtomSluMol = 0
    do i = 1, nAtmSlu
       if ((index(atomSlu(i),'H') /= 1)) then
          nHvyAtomSluMol = nHvyAtomSluMol + 1
       end if
    end do
    nHvyAtomSluSys = nHvyAtomSluMol * nMolSlu
    
    do while ( xtc % STAT == 0 )
       ng = ng + 1
       write (*,*), ng

       posSlu = xtc % pos(:,IndxIniAtmSlu:IndxFinAtmSlu)
       posSlv = xtc % pos(:,IndxIniAtmSlv:IndxFinAtmSlv)

       call get_com(nAtmSlu, posSlu, massSlu, com, atomSlu)

       write(*,*), com
       
!       m = 0
       do j = 1, nMolSlv
          do l = 1, nAtmSlv 
             if ((index(atomSlv(l),'H') /= 1)) then

                ! write(*,*), atomSlv(l)

!                m = m + 1
                
                xr = com(:) - posSlv(:,l+(j-1)*nAtmSlv)
                !wrap distance (periodic boundary condition)
                xr = abs(xr - box_dim * nint(xr / box_dim))
                dxr = sqrt(sum(xr**2))
                   
                ig = ceiling(dxr / dr)
                if (ig <= nhist) then
                   if (ig == 0) then
                      ig = 1  !index of g(:) begins from 1
                   end if
                   gpcom(ig) = gpcom(ig) + 1
                end if
             end if
          end do
       end do
!       write(*,*), m
    
       do i = 1, nMolSlu
          do j = 1, nMolSlv
             min_dxr = 100.0d0
             do k = 1, nAtmSlu
                do l = 1, nAtmSlv
                   if ((index(atomSlu(k),'H') /= 1) &
                        .and. (index(atomSlv(l),'H') /= 1)) then
                      xr = posSlu(:,k+(i-1)*nAtmSlu) - posSlv(:,l+(j-1)*nAtmSlv)
                      !wrap distance (periodic boundary condition)
                      xr = abs(xr - box_dim * nint(xr / box_dim))
                      dxr = sqrt(sum(xr**2))

                      if (dxr < min_dxr ) then
                         min_dxr = dxr
                      end if

                      ig = ceiling(dxr / dr)
                      if (ig <= nhist) then
                         if (ig == 0) then
                            ig = 1  !index of g(:) begins from 1
                         end if
                         ghat(ig) = ghat(ig) + 1
                      end if

                   end if
                end do
             end do

             ig = ceiling(min_dxr / dr)
             if (ig <= nhist) then
                if (ig == 0) then
                   ig = 1  !index of g(:) begins from 1
                end if
                gminp(ig) = gminp(ig) + 1
             end if

          end do
       end do
    
       call xtc % read
     end do

    ! 5. Close the file
    call xtc % close

    ! normalize rdf
    rho = dble(nMolSlv) / box_dim**3  !number density
    rho2 = dble(nHvyAtomSlvSys) / box_dim**3  !number density
    do i = 1, nhist
       dv = (4d0 / 3d0) * pi * (i**3 - (i-1)**3) * dr**3
       ghat(i)  =  ghat(i) / (ng * dble(nHvyAtomSlvMol) * dble(nHvyAtomSluSys) * dv *rho)
       gminp(i) =  gminp(i) / (ng * dble(nMolSlu) * dble(nMolSlv) ) !  * dv *rho
       gpcom(i) =  gpcom(i) / (ng * dv  * rho2)
    end do

    ! output results
    open(19,file=rdf_filename,position='append')
    do i = 1, nhist
       write(19,'(4f20.8)') r(i), ghat(i), gminp(i), gpcom(i)
    end do
    close(19)

    deallocate(atomSlu)
    deallocate(atomSlv)
    
    deallocate(posSlu)
    deallocate(posSlv)

    deallocate(g)
    deallocate(gminp)
    deallocate(ghat)
    deallocate(gpcom)
    deallocate(r)

    deallocate(massSlu)
    deallocate(massSlv)
    
contains
  subroutine usage(programname)
    implicit none
    character(10),intent(in) :: programname

    print *,"Usage: ",programname," inpfilename outfilename itpslufilename itpslvfilename"
    stop
  end subroutine usage

  subroutine get_com(natom, coord, mass, com, atom)
    implicit none

    integer, intent(in) :: natom
    double precision, intent(in)  :: coord(3, natom), mass(natom)
    character(5), allocatable,dimension(:), intent(in) :: atom
    double precision, intent(out) :: com(3)
    double precision :: summass

    integer :: i

    summass = 0d0
    do i = 1, natom
       if ((index(atom(i),'H') /= 1)) then
          summass = summass + mass(i)
       end if
    end do

    com = 0d0
    do i = 1, natom
       if ((index(atom(i),'H') /= 1)) then
          com(1) = com(1) + coord(1,i)*mass(i)
          com(2) = com(2) + coord(2,i)*mass(i)
          com(3) = com(3) + coord(3,i)*mass(i)
       end if
    end do
    com(1) = com(1) / summass
    com(2) = com(2) / summass
    com(3) = com(3) / summass

  end subroutine get_com
  
end program g_df
