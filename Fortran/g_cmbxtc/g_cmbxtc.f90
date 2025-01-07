! This program combines the two trajectry files (xtc)
! to one file.
program g_cmbxtc
    use xdr, only: xtcfile
    implicit none
    type(xtcfile) :: xtca, xtcb
    type(xtcfile) :: xtcab

    integer :: step_a, step
    real :: time_a, time
    
    character xtca_filename*50, xtcb_filename*50, xtcab_filename*50

    character(10) :: argv, programname
    integer :: argc, iargc

    argc=iargc()
    call getarg(0,programname)
    if (argc < 3) then
       call usage(programname)
    else
       call getarg(1,xtca_filename)
       call getarg(2,xtcb_filename)
       call getarg(3,xtcab_filename)
    end if

    call xtca % init(xtca_filename)
    call xtca % read

    call xtcab % init(xtcab_filename, 'w')

    do while ( xtca % STAT == 0 )
       step = xtca % step
       time = xtca % time
       call xtcab % write(xtca % natoms, xtca % step, xtca % time, xtca % box, xtca % pos, xtca % prec)
       
       call xtca % read
    end do
    
    call xtcb % init(xtcb_filename)
    call xtcb % read

    do while ( xtcb % STAT == 0 )
       step = step_a + xtca % step
       time = time_a + xtca % time
       call xtcab % write(xtcb % natoms, step, time, xtcb % box, xtcb % pos, xtcb % prec)

       call xtcb % read
    end do
    
    call xtca % close
    call xtcb % close
    call xtcab % close
    
contains
  subroutine usage(programname)
    implicit none
    character(10),intent(in) :: programname

    print *,"Usage: ",programname," xtcfilenamea xtcfilenameb xtcfilenameab"
    stop
  end subroutine usage
end program g_cmbxtc










