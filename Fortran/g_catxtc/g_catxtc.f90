! This program convert the trajectry file (xtc)
! for desired order of components
program g_catxtc
    use xdr, only: xtcfile
    implicit none
    type(xtcfile) :: xtc

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

    do while ( xtc % STAT == 0 )
       ! 3. write xtc-file
       write(*,'(a,f12.6,a,i0)') " Time (ps): ", xtc % time, "  Step: ", xtc % STEP
       write(*,'(a,f12.6,a,i0)') " Precision: ", xtc % prec, "  No. Atoms: ", xtc % NATOMS
       write(*,'(3f9.3)') xtc % pos

       call xtc % read
    end do
    
    ! 4. Close the file
    call xtc % close
    
contains
  subroutine usage(programname)
    implicit none
    character(10),intent(in) :: programname

    print *,"Usage: ",programname," xtcfilename"
    stop
  end subroutine usage
end program g_catxtc
