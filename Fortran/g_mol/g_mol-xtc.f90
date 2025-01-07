! This program calculates the trajectry of mol (density)
! in NPT ensemble with periodic cubic boundaries
program g_mol
    use xdr, only: xtcfile
    implicit none
    type(xtcfile) :: xtc
    real :: mol, V, box_dim 
    real :: ave_mol
    real :: Anum = 6.02214129e23
    real :: nm3tol = 1.0e-24
    real :: nm3toM
    integer :: nMol
    integer :: i, j

    character xtc_filename*50
    character mol_filename*50

    character(10) :: argv, programname
    integer :: argc, iargc

    ! 0 reading parameters
    NAMELIST/controls/nMol
    open(17,file='parameters_g_mol')
    read(17,nml=controls)
    close(17)

    argc=iargc()
    call getarg(0,programname)
    if (argc < 2) then
       call usage(programname)
    else
       call getarg(1,xtc_filename)
       call getarg(2,mol_filename)
    end if

    call xtc % init(xtc_filename)

    write(*,*),"nMol "
    write(*,*),nMol

    open(20,file=mol_filename,position='append')

    nm3toM = 1.0/(Anum*nm3tol)

    write(*,'(F12.8)'),nm3toM

    call xtc % read
    
    i = 0
    do while ( xtc % STAT == 0 )
       box_dim = xtc % box(1,1)
       V = box_dim**3  !number density
       mol = nMol / V * nm3toM
       ave_mol = ave_mol + mol
       i = i + 1
       write(20,'(I4XF12.8XF12.8XF12.8)'),i,box_dim, V, mol

       call xtc % read
    end do

    ave_mol = ave_mol / i

    i = i + 1
    write(20,'(I4X2F12.8)'),i, V, mol

    ! 5. Close the file
    call xtc % close

   close(20)

contains
  subroutine usage(programname)
    implicit none
    character(10),intent(in) :: programname

    print *,"Usage: ",programname," inputfilename outputfilename"
    stop
  end subroutine usage
end program g_mol
