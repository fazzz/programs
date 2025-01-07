!  XDR Fortran Interface XTC Example Program with Wrapper
!  2014 (c) James W. Barnett <jbarnet4@tulane.edu>
!  https://github.com/wesbarnett/
!

program example

    ! 1. Use the xdr interface
    use xdr, only: xtcfile

    implicit none

    ! 2. Declare a variable of type xtcfile
    type(xtcfile) :: xtcf

    character xtc_filename*50

    character(10) :: argv, programname
    integer :: argc, iargc
    
    argc=iargc()
    call getarg(0,programname)
    if (argc < 1) then
       call usage(programname)
    else
       call getarg(1,xtc_filename)
    end if
    
    ! 3. Initialize it with the names of xtc files you want to read in and write out
    call xtcf % init(xtc_filename)

    ! 4. Read in each configuration. Everything is stored in the xtcfile type (precision, time,
    !    step, no of atoms, positions, etc.). Look in the xtc module for more details.
    !    You can save the positions in the loop for your calculations in another array, or 
    !    do your calculations after each read.

    call xtcf % read

    write(*,'(3f9.3)') xtcf % pos

    ! 5. Close the file
    call xtcf % close

contains
    subroutine usage(programname)
    implicit none
    character(10),intent(in) :: programname

    print *,"Usage: ",programname," xtcfilename"
    stop
  end subroutine usage

end program example
