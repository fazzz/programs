#include "copyright.i"
!*******************************************************************************
!
! Module: get_cmdline_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module get_cmdline_mod

  implicit none

  private       add_suffix

contains

#ifdef USE_PXFGETARG
!*******************************************************************************
!
! Subroutine:  pmemd_getarg
!
! Description: <TBS>
!
!*******************************************************************************

subroutine pmemd_getarg(i, c)

  implicit none

! Formal arguments:

  character(*)  :: c
  integer       :: i

! Local variables:

  integer       :: ilen
  integer       :: ierror

! Note that the pxfgetarg subroutine is broken in ifc 7.1; use getarg instead
! for any compiler that supports it!

  call pxfgetarg(i, c, ilen, ierror)

  if (ierror .ne. 0) then
    write(mdout, '(a)') 'ERROR: Command line parsing failed in getarg!'
    call mexit(6, 1)
  end if

  return

end subroutine pmemd_getarg
#else
! The preferred command line subroutine is getarg.
#define pmemd_getarg    getarg
#endif /* USE_PXFGETARG */

!*******************************************************************************
!
! Subroutine:  get_cmdline
!
! Description: <TBS>
!              
! OUTPUT: (to common)
!
! Author: George Seibel
!
!*******************************************************************************

subroutine get_cmdline

  use file_io_dat_mod
  use pmemd_lib_mod

  implicit none

! External system function declaration:

#ifdef USE_PXFGETARG
  integer       :: ipxfargc
#else
  integer       :: iargc
#endif

! Local variables:

  character(80) :: arg  ! temp buffer for each of the whitespace-
                        ! delimited command line words.
  integer       :: iarg ! arg pointer, final number of arguments.
  integer       :: indx

! Flags that determine if/when the -suffix gets used:

  logical       :: mdout_specified
  logical       :: mdinfo_specified
  logical       :: mdcrd_specified
  logical       :: mdvel_specified
  logical       :: mden_specified
  logical       :: restrt_specified
  logical       :: logfile_specified
  logical       :: outfile_suffix_specified

  integer       :: suffix_len

  mdout_specified = .false.
  mdinfo_specified = .false.
  mdcrd_specified = .false.
  mdvel_specified = .false.
  mden_specified = .false.
  restrt_specified = .false.
  logfile_specified = .false.
  outfile_suffix_specified = .false.

! Default file names:

  mdin_name   = 'mdin'
  mdout_name  = 'mdout'
  inpcrd_name = 'inpcrd'
  prmtop_name = 'prmtop'
  restrt_name = 'restrt'
  refc_name   = 'refc'
  mdvel_name  = 'mdvel'
  mden_name   = 'mden'
  mdcrd_name  = 'mdcrd'
  mdinfo_name = 'mdinfo'
  logfile_name = 'logfile'

! Default status of output: new

  owrite = 'N'

! Get com line arguments:

  iarg = 0
#ifdef USE_PXFGETARG
  indx = ipxfargc()
#else
  indx = iargc()
#endif

  if (indx .eq. 0) goto 20

  10 continue

  iarg = iarg + 1

  call pmemd_getarg(iarg, arg)

  if (arg .eq. '-O') then
    owrite = 'U'
  else if (arg .eq. '-i') then
    iarg = iarg + 1
    call pmemd_getarg(iarg, mdin_name)
  else if (arg .eq. '-o') then
    iarg = iarg + 1
    call pmemd_getarg(iarg, mdout_name)
    mdout_specified = .true.
  else if (arg .eq. '-p') then
    iarg = iarg + 1
    call pmemd_getarg(iarg, prmtop_name)
  else if (arg .eq. '-c') then
    iarg = iarg + 1
    call pmemd_getarg(iarg, inpcrd_name)
  else if (arg .eq. '-r') then
    iarg = iarg + 1
    call pmemd_getarg(iarg, restrt_name)
    restrt_specified = .true.
  else if (arg .eq. '-ref' .or. arg .eq.'-z') then
    iarg = iarg + 1
    call pmemd_getarg(iarg, refc_name)
  else if (arg .eq. '-e') then
    iarg = iarg + 1
    call pmemd_getarg(iarg, mden_name)
    mden_specified = .true.
  else if (arg .eq. '-v') then
    iarg = iarg + 1
    call pmemd_getarg(iarg, mdvel_name)
    mdvel_specified = .true.
  else if (arg .eq. '-x' .or. arg .eq.'-t') then
    iarg = iarg + 1
    call pmemd_getarg(iarg, mdcrd_name)
    mdcrd_specified = .true.
  else if (arg .eq. '-inf') then
    iarg = iarg + 1
    call pmemd_getarg(iarg, mdinfo_name)
    mdinfo_specified = .true.
  else if (arg .eq. '-l') then
    iarg = iarg + 1
    call pmemd_getarg(iarg, logfile_name)
    logfile_specified = .true.
  else if (arg .eq. '-suffix') then
    iarg = iarg + 1
    call pmemd_getarg(iarg, outfile_suffix)
    outfile_suffix_specified = .true.
  else if (arg .eq. '-help') then
    write(mdout, 9000)
    call mexit(6, 0)

#ifdef MPI
  else if (arg .eq. '-p4pg') then
    iarg = iarg + 1
  else if (arg .eq. '-p4wd') then
    iarg = iarg + 1
  else if (arg .eq. '-np') then
    iarg = iarg + 1
  else if (arg .eq. '-mpedbg') then
    continue
  else if (arg .eq. '-dbx') then
    continue
  else if (arg .eq. '-gdb') then
    continue
#endif
  else
    if (arg .eq. ' ') goto 20
    write(mdout, '(/,2x,a,a)') 'unknown flag: ', arg
    write(mdout, 9000)
    call mexit(6, 1)
  end if

  if (iarg .lt. indx) goto 10

  20 continue

  if (outfile_suffix_specified) then

    suffix_len = len_trim(outfile_suffix)

    if (.not. mdout_specified) &
      call add_suffix(mdout_name, outfile_suffix, suffix_len)

    if (.not. mdinfo_specified) &
      call add_suffix(mdinfo_name, outfile_suffix, suffix_len)

    if (.not. mdcrd_specified) &
      call add_suffix(mdcrd_name, outfile_suffix, suffix_len)

    if (.not. mdvel_specified) &
      call add_suffix(mdvel_name, outfile_suffix, suffix_len)

    if (.not. mden_specified) &
      call add_suffix(mden_name, outfile_suffix, suffix_len)

    if (.not. restrt_specified) &
      call add_suffix(restrt_name, outfile_suffix, suffix_len)

    if (.not. logfile_specified) &
      call add_suffix(logfile_name, outfile_suffix, suffix_len)

  end if

  return

9000 format(/, 2x, &
      'usage: pmemd  [-O] -i mdin -o mdout -p prmtop -c inpcrd -r restrt', &
      /14x, '[-ref refc -x mdcrd -v mdvel -e mden -inf mdinfo -l logfile]', &
      /14x, '[-suffix output_files_suffix]',/)

end subroutine get_cmdline

!*******************************************************************************
!
! Subroutine:  add_suffix
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine add_suffix(file_name, suffix, suffix_len)

  use file_io_dat_mod
  use pmemd_lib_mod

  implicit none

! Formal arguments:

  character(max_fn_len), intent(in out) :: file_name
  character(max_fn_len), intent(in)     :: suffix
  integer, intent(in)                   :: suffix_len

  if (suffix_len + len_trim(file_name) .lt. max_fn_len) then
    file_name = trim(file_name) // '.' // trim(suffix)
  else
    write(mdout, '(a)') 'ERROR: Filename with suffix too long!'
    call mexit(6, 1)
  end if

  return

end subroutine add_suffix

end module get_cmdline_mod
