! This program convert the trajectry file (xtc)
! for desired order of components
program g_trjchang
    use xdr, only: xtcfile
    implicit none
    type(xtcfile) :: xtc, xtc_out
    integer :: i, j, k
    integer :: nAtom1, nAtom2
    integer :: nMol1, nMol2
    real, allocatable,dimension(:,:) ::  pos

    character xtc_filename*50
    character xtc_out_filename*50

    character(10) :: argv, programname
    integer :: argc, iargc

    ! 0 reading parameters
    NAMELIST/controls/nAtom1, nAtom2, nMol1, nMol2
    open(17,file='parameters_g_trjchang')
    read(17,nml=controls)
    close(17)

    ! 1 getting file names
    argc=iargc()
    call getarg(0,programname)
    if (argc < 2) then
       call usage(programname)
    else
       call getarg(1,xtc_filename)
       call getarg(2,xtc_out_filename)
    end if

    call xtc % init(xtc_filename)
    call xtc_out % init(xtc_out_filename,"w")

    call xtc % read

    allocate(pos(3, xtc % natoms))

    do while ( xtc % STAT == 0 )
       do i = 1, nMol1*nAtom1
          pos(:,nMol2*nAtom2+i) = xtc % pos(:,i)
       end do

       do i = 1, nMol2*nAtom2
          pos(:,i) = xtc % pos(:,nMol1*nAtom1+i)
       end do

       call xtc_out % write(xtc % natoms, xtc % step,  &
            xtc % time, xtc % box, pos, xtc % prec)

       call xtc % read
    end do

    deallocate(pos)
    
    ! 5. Close the file
    call xtc % close
    call xtc_out % close
    
contains
  subroutine usage(programname)
    implicit none
    character(10),intent(in) :: programname

    print *,"Usage: ",programname," inputfilename outputfilename"
    stop
  end subroutine usage
end program g_trjchang
