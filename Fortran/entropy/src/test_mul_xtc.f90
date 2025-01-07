program bestfit_ty

use xdr, only: xtcfile
use bestfit
	implicit none
	type(xtcfile) :: xtcbf
	type(xtcfile) :: xtc_out
	integer :: i, j, slvnum, refslt_natom, bfslt_natom, natom
	integer, dimension(:), allocatable :: num
	character , dimension(:), allocatable :: eltp1
	double precision, allocatable, dimension(:) :: stmass, charge, ljep, ljsi
	double precision, allocatable, dimension(:,:) :: refcoord, outcoord, bfcoord
	real, allocatable, dimension(:,:) :: pos
	character reffilename*50, bffilename*50, outfilename*50, proname*20, SltInfo*20
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

call inputfile(refslt_natom, slvnum, proname, natom, refcoord)
call calcmass(refslt_natom, stmass)
call calcfit(refslt_natom, slvnum, stmass, refcoord, bfcoord, outcoord)
!call output(bffilename, refslt_natom, outcoord)

contains
subroutine inputfile(refslt_natom, slvnum, proname, natom, refcoord)
	implicit none
	integer :: i, slvnum, refslt_natom, natom
	double precision, allocatable, dimension(:,:) :: refcoord
	character reffilename*50, proname*20
	print *, "The number of water molecule is"
	read *, slvnum
	print *, "Input reffile name is"
	read *, reffilename
	open(18, file = reffilename, status = 'old')
	read(18, '(20A)'), proname
	read(18, *), natom
	refslt_natom = natom - slvnum * 3
	allocate( refcoord(3, refslt_natom) )
   	do i = 1, refslt_natom
   	read(18, '(20X, 3F8.3)'), refcoord(1, i), refcoord(2, i), refcoord(3, i)
   	enddo
	close(18)
end subroutine inputfile

subroutine calcfit(refslt_natom, slvnum, stmass, refcoord, bfcoord, outcoord)
   integer :: slvnum, refslt_natom, bfslt_natom
   double precision, allocatable, dimension(:,:) :: refcoord, bfcoord, outcoord
	double precision, allocatable, dimension(:) :: stmass
	real, allocatable, dimension(:,:) :: pos
   character bffilename*50, outfilename*50
	print *, "Input bffile name is"
	read *, bffilename
	print *, "Outputfile name is"
	read *, outfilename
	call xtcbf % init( bffilename )
	call xtc_out % init(outfilename, 'w')
	call xtcbf % read
		bfslt_natom = xtcbf % NATOMS - slvnum * 3
		if(refslt_natom /= bfslt_natom) then
			print *, "Error: Inputfile is not correct"
			stop
		endif
		do while (xtcbf % STAT == 0)
				bfcoord = xtcbf % pos(:, 1:refslt_natom)
				allocate( outcoord(3, refslt_natom) )
			call fit(refslt_natom, refcoord, bfcoord, stmass, outcoord)
				allocate( pos(3, refslt_natom) )
				pos(:,:) = real( outcoord(1:3, 1:refslt_natom) )
			call xtc_out % write(refslt_natom, xtcbf % step, xtcbf % time, &
xtcbf % box, pos, xtcbf % prec)
			call xtcbf % read
				deallocate( outcoord )
				deallocate( pos )
		enddo
			call xtc_out % close
			call xtcbf % close
end subroutine calcfit

subroutine calcmass(refslt_natom, stmass)
	implicit none
	integer :: i, refslt_natom
	integer, dimension(:), allocatable :: num
	character , dimension(:), allocatable :: eltp1
	double precision, allocatable, dimension(:) :: stmass, charge, ljep, ljsi
	character SltInfo*20
		open(19, file = 'SltInfo', status = 'old')
		allocate( stmass(refslt_natom), num(refslt_natom), eltp1(refslt_natom), &
charge(refslt_natom), ljep(refslt_natom), ljsi(refslt_natom) )
			do i = 1, refslt_natom
				read(19, *), num(i), eltp1(i), charge(i), ljep(i), ljsi(i)
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
		close(19)
end subroutine calcmass

end program bestfit_ty
