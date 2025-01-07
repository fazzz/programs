program bestfit_ty

use xdr, only: xtcfile
use bestfit
	implicit none
	type(xtcfile) :: xtcbf
	type(xtcfile) :: xtc_out
	integer :: h, i, j, k, l, n, frame, slvnum, refslt_natom, bfslt_natom
	integer, dimension(:), allocatable :: num
	character , dimension(:), allocatable :: eltp1
	double precision :: center(3)
	double precision, allocatable, dimension(:) :: stmass, charge, ljep, ljsi
	double precision, allocatable, dimension(:,:) :: outcoord, bfcoord, sumcoord, meancoord
	double precision, allocatable, dimension(:,:,:) :: refcoord, cencoord
	real, allocatable, dimension(:,:) :: pos
	character bffilename*50, outfilename*50, SltInfo*20
	real, parameter :: massH = 1.00794          ! atomic weight (hydrogen)
	real, parameter :: massC = 12.0107          ! atomic weight (carbon)
	real, parameter :: massO = 15.9994          ! atomic weight (oxygen)
	real, parameter :: massN = 14.00674         ! atomic weight (nitrogen)
	real, parameter :: massS = 32.066           ! atomic weight (sulfur)
	real, parameter :: massP = 30.973761        ! atomic weight (phosphorus)
	real, parameter :: massB = 10.811           ! atomic weight (boron)
	real, parameter :: massK = 39.0983          ! atomic weight (potassium)
	real, parameter :: massF = 18.9984032       ! atomic weight (fluorine)
	real, parameter :: massI = 126.90447        ! atomic weight (iodine)

call calcmass(refslt_natom, stmass)
call inputfile(bffilename, frame, refslt_natom, refcoord)
call calcmean(refslt_natom, frame, refcoord, stmass, center, cencoord)
call calcfit(bffilename, refslt_natom, stmass, meancoord, bfcoord, cencoord)

contains
subroutine calcmass(refslt_natom, stmass)
   implicit none
   integer :: i, j, refslt_natom
   integer, dimension(:), allocatable :: num
   character , dimension(:), allocatable :: eltp1
   double precision, allocatable, dimension(:) :: stmass, charge, ljep, ljsi
   character SltInfo*20

      open(19, file = 'SltInfo', status = 'old')
         j = 1
         do
            read (19, *, end = 999)
            j = j + 1
         enddo
      999 refslt_natom = j - 1
      close(19)
      open(21, file = 'SltInfo', status = 'old')
         allocate( stmass(refslt_natom), num(refslt_natom), eltp1(refslt_natom), &
charge(refslt_natom), ljep(refslt_natom), ljsi(refslt_natom) )
         do i = 1, refslt_natom
            read(21, *), num(i), eltp1(i), charge(i), ljep(i), ljsi(i)
            if(eltp1(i) == 'H') stmass(i) = massH
            if(eltp1(i) == 'C') stmass(i) = massC
            if(eltp1(i) == 'O') stmass(i) = massO
            if(eltp1(i) == 'N') stmass(i) = massN
            if(eltp1(i) == 'S') stmass(i) = massS
            if(eltp1(i) == 'P') stmass(i) = massP
            if(eltp1(i) == 'B') stmass(i) = massB
            if(eltp1(i) == 'K') stmass(i) = massK
            if(eltp1(i) == 'F') stmass(i) = massF
            if(eltp1(i) == 'I') stmass(i) = massI
         enddo
      close(21)
end subroutine calcmass

subroutine inputfile(bffilename, frame, refslt_natom, refcoord)
	implicit none
	integer :: i, frame, slvnum, refslt_natom, bfslt_natom
	double precision, allocatable, dimension(:) :: stmass
	double precision, allocatable, dimension(:,:,:) :: refcoord
	character bffilename*50

	print *, "The number of frame is"
	read *, frame
	print *, "The number of water molecule is"
	read *, slvnum
	print *, "Input reffile name is"
	read *, bffilename
   call xtcbf % init( bffilename )
   call xtcbf % read
      bfslt_natom = xtcbf % NATOMS - slvnum * 3
   if(refslt_natom /= bfslt_natom) then
      print *, "Error:Inputfile is not correct"
      stop
   endif
!	do while ( xtcbf % STAT == 0 )
!		call xtcbf % read
!		frame = xtcbf % STEP
!	enddo
	call xtcbf % close
	call xtcbf % init( bffilename )
	call xtcbf % read
   allocate( refcoord(3, refslt_natom, frame) )
	do i = 1, frame
      refcoord(:,:,i) = xtcbf % pos(:, 1:refslt_natom)
		call xtcbf % read
	enddo
	call xtcbf % close
end subroutine inputfile

subroutine calcmean(refslt_natom, frame, refcoord, stmass, center, cencoord)
	implicit none
	integer :: h, i, j, k, l, n, frame, refslt_natom
	double precision :: center(3)
   double precision, allocatable, dimension(:) :: stmass
	double precision, allocatable, dimension(:,:) :: sumcoord, meancoord, outcoord, bfcoord
   double precision, allocatable, dimension(:,:,:) :: refcoord, cencoord
	call center_mass(refslt_natom, frame, refcoord, stmass, center)
	allocate( meancoord(3, refslt_natom), cencoord(3, refslt_natom, frame), &
sumcoord(3, refslt_natom), outcoord(3, refslt_natom), bfcoord(3, refslt_natom) )
	do j = 1, frame
		do i = 1, refslt_natom
			cencoord(1,i,j) = refcoord(1,i,j) - center(1)
			cencoord(2,i,j) = refcoord(2,i,j) - center(2)
			cencoord(3,i,j) = refcoord(3,i,j) - center(3)
		enddo
	enddo
	do k = 1, refslt_natom
			sumcoord(1,k) = sum(cencoord(1,k,:))
			sumcoord(2,k) = sum(cencoord(2,k,:))
			sumcoord(3,k) = sum(cencoord(3,k,:))
	enddo
			meancoord(:,:) = sumcoord(1:3,1:refslt_natom) / frame
	deallocate( cencoord )
	allocate( cencoord(3, refslt_natom, frame) )
	do n = 1, frame
			bfcoord(:,:) = refcoord(1:3,1:refslt_natom,n)
			call fit(refslt_natom, meancoord, bfcoord, stmass, outcoord)
			cencoord(:,:,n) = outcoord(1:3,1:refslt_natom)
	enddo

	do h = 1, 3
		do l = 1, refslt_natom
         sumcoord(1,l) = sum(cencoord(1,l,:))
         sumcoord(2,l) = sum(cencoord(2,l,:))
         sumcoord(3,l) = sum(cencoord(3,l,:))
		enddo
         meancoord(:,:) = sumcoord(1:3,1:refslt_natom) / frame
		deallocate( outcoord, cencoord )
		allocate( outcoord(3, refslt_natom), cencoord(3, refslt_natom, frame) )
		call fit(refslt_natom, meancoord, cencoord, stmass, outcoord)
		cencoord(:,:,h) = outcoord(1:3,1:refslt_natom)
	enddo
end subroutine calcmean

subroutine calcfit(bffilename, refslt_natom, stmass, meancoord, bfcoord, cencoord)
	implicit none
   integer :: i, j, refslt_natom
	double precision, allocatable, dimension(:) :: stmass
   double precision, allocatable, dimension(:,:) :: bfcoord, outcoord, meancoord
	double precision, allocatable, dimension(:,:,:) :: refcoord, cencoord
	real, allocatable, dimension(:,:) :: pos
   character bffilename*50, outfilename*50

	print *, "Outputfile name is"
	read *, outfilename
   call xtcbf % init( bffilename )
   call xtc_out % init(outfilename, 'w')
   call xtcbf % read
		do i = 1, frame
		allocate( pos(3, refslt_natom) )
			pos(:,:) = real( cencoord(1:3, 1:refslt_natom, i) )
			call xtc_out % write(refslt_natom, xtcbf % step, xtcbf % time, &
xtcbf % box, pos, xtcbf % prec)
		deallocate( pos )
		enddo
	call xtc_out % close
	call xtcbf % close
end subroutine calcfit

subroutine center_mass(refslt_natom, frame, refcoord, stmass, center)
    implicit none
    integer, intent(in) :: refslt_natom
    double precision, intent(in) :: refcoord(3, refslt_natom, frame), stmass(refslt_natom)
    double precision, intent(out) :: center(3)
    double precision :: summass
    integer :: i, j, frame

    summass = sum(stmass)
    if(summass == 0) write(*,*),"error"
	do j = 1, frame
    do i = 1, 3
       center(i) = dot_product(refcoord(i,:,j), stmass(:)) / summass
    end do
	enddo
end subroutine center_mass

end program bestfit_ty
