program bestfit_ty

use bestfit
implicit none
integer :: i, slvnum, refslt_natom, natom, natom2
integer, dimension(:), allocatable :: num
character , dimension(:), allocatable :: eltp1
double precision, allocatable, dimension(:) :: stmass, charge, ljep, ljsi
double precision, allocatable, dimension(:,:) :: refcoord, outcoord, bfcoord
character reffilename*50, bffilename*50, filename*50, proname*20, proname2*20, SltInfo*20
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

call inputfile(refslt_natom, proname, natom, refcoord, bfcoord)
call calcmass(refslt_natom, stmass)
allocate( outcoord(3, refslt_natom) )
call fit(refslt_natom, refcoord, bfcoord, stmass, outcoord)
call outputfile(proname, natom, refslt_natom, outcoord)

contains
subroutine inputfile(refslt_natom, proname, natom, refcoord, bfcoord)
implicit none
integer :: i, slvnum, refslt_natom, natom, natom2
double precision, allocatable, dimension(:,:) :: refcoord, bfcoord
character reffilename*50, bffilename*50, filename*50, proname*20, proname2*20
print *, "The number of water molecule is"
read *, slvnum

print *, "Input reffile name is"
read *, reffilename
open(17, file = reffilename, status = 'old')
read(17, '(20A)'), proname
read(17, *), natom
refslt_natom = natom - slvnum * 3
allocate( refcoord(3, refslt_natom) )
   do i = 1, refslt_natom
	read(17, '(20X, 3F8.3)'), refcoord(1, i), refcoord(2, i), refcoord(3, i)
   enddo
close(17)

print *, "Input bffile name is"
read *, bffilename
open(18, file = bffilename, status = 'old')
read(18, '(20A)'), proname2
read(18, *), natom2
allocate( bfcoord(3, refslt_natom) )
   do i = 1, refslt_natom
	read(18, '(20X, 3F8.3)'), bfcoord(1, i), bfcoord(2, i), bfcoord(3, i)
   enddo
close(18)
if(proname /= proname2 .or. natom /= natom2) then
   print *, "Error:Inputfile is not correct"
   stop
endif
end subroutine inputfile

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

subroutine outputfile(proname, natom, refslt_natom, outcoord)
implicit none
integer :: i, refslt_natom, natom
double precision, allocatable, dimension(:,:) :: outcoord
character proname*20
print *, "Outputfile name is"
   read *, filename
   open(20, file = filename, status = 'replace')
      write(20, '(20A)'), proname
      write(20, *), natom
      do i = 1, refslt_natom
		write(20,*), outcoord(1, i), outcoord(2, i), outcoord(3, i)
      enddo
   close(20)
end subroutine outputfile

end program bestfit_ty

