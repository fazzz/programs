! This program convert the trajectry file (xtc)
! for desired order of components
program g_showxtc
    use xdr, only: xtcfile
    implicit none
    type(xtcfile) :: xtc

    integer i
    
    character xtc_filename*50

    character(10) :: argv, programname
    integer :: argc, iargc

    ! 1 getting file names
    argc=iargc()
    call getarg(0,programname)
    if (argc < 1) then
       call usage(programname)
    else
       call getarg(1,xtc_filename)
    end if

    ! 2. Init xtc-file
    call xtc % init(xtc_filename)
    call xtc % read
    ! 3. write xtc-file
    write(*,'(a,f12.6,a,i0)') " Precision: ", xtc % prec, "  No. Atoms: ", xtc % NATOMS

    i = 0
    do while ( xtc % STAT == 0 )
       call xtc % read
       i = i + 1
    end do
    i = i - 1
    
    ! 4. Count steps
!    write(*,'(a,f12.6,a,i0)') " Time (ps): ", xtc % time, "  Step: ", xtc % STEP
!    write(*,'(a,i0)') "  Step: ", xtc % STEP
     write(*,'(a,i0)') "  Frames: ", i
    
    ! 4. Close the file
    call xtc % close
    
contains
  subroutine usage(programname)
    implicit none
    character(10),intent(in) :: programname

    print *,"Usage: ",programname," xtcfilename"
    stop
  end subroutine usage
end program g_showxtc
