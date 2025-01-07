program bestfit_ty

use xdr, only: xtcfile
use bestfit
        implicit none
        type(xtcfile) :: xtcbf
        type(xtcfile) :: xtc_out
        integer :: i, frame, refslt_natom
        double precision :: ev(3), center(3), varm(3,3), corm(3,3), vec(3,3)
        double precision, allocatable, dimension(:) :: stmass
        double precision, allocatable, dimension(:,:) :: xd ! Array for inputfile data
        double precision, allocatable, dimension(:,:,:) :: refcoord, cencoord
        character bffilename*50
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

        call inputmass( refslt_natom, stmass )
        call inputfile( bffilename, frame, refslt_natom, refcoord )
        call calcmean( refslt_natom, frame, refcoord, stmass, center, cencoord )
        call outputfile( bffilename, refslt_natom, cencoord )
        allocate( xd(3, refslt_natom))
        do i = 1, frame
           xd(:,:) = cencoord(1:3, 1:refslt_natom, i)
           call varcor( refslt_natom, xd, varm, corm )
           call jacobi( varm, ev, vec )
           print *, ev(1:3)
        enddo

contains
!Input the number of solute and masses of all atoms
subroutine inputmass(refslt_natom, stmass)
        implicit none
        integer :: i, j, refslt_natom
        integer, dimension(:), allocatable :: num
        character , dimension(:), allocatable :: eltp
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
          allocate( stmass(refslt_natom), num(refslt_natom), eltp(refslt_natom),&
charge(refslt_natom), ljep(refslt_natom), ljsi(refslt_natom) )
          do i = 1, refslt_natom
             read(21, *), num(i), eltp(i), charge(i), ljep(i), ljsi(i)
             if(eltp(i) == 'H') stmass(i) = massH
             if(eltp(i) == 'C') stmass(i) = massC
             if(eltp(i) == 'O') stmass(i) = massO
             if(eltp(i) == 'N') stmass(i) = massN
             if(eltp(i) == 'S') stmass(i) = massS
             if(eltp(i) == 'P') stmass(i) = massP
             if(eltp(i) == 'B') stmass(i) = massB
             if(eltp(i) == 'K') stmass(i) = massK
             if(eltp(i) == 'F') stmass(i) = massF
             if(eltp(i) == 'I') stmass(i) = massI
          enddo
        close(21)
end subroutine inputmass

!Read inputfile data
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
           allocate( refcoord(3, refslt_natom, frame) )
        do i = 1, frame
           refcoord(:,:,i) = xtcbf % pos(:, 1:refslt_natom)
           call xtcbf % read
        enddo
        call xtcbf % close
end subroutine inputfile

!Calcurate mean structure
subroutine calcmean(refslt_natom, frame, refcoord, stmass, center, cencoord)
        implicit none
        integer :: i, j, frame, refslt_natom
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
        do i = 1, refslt_natom
                sumcoord(1,i) = sum(cencoord(1,i,:))
                sumcoord(2,i) = sum(cencoord(2,i,:))
                sumcoord(3,i) = sum(cencoord(3,i,:))
        enddo
                meancoord(:,:) = sumcoord(1:3,1:refslt_natom) / frame
        deallocate( cencoord )
        allocate( cencoord(3, refslt_natom, frame) )
        do i = 1, frame
                bfcoord(:,:) = refcoord(1:3,1:refslt_natom,i)
                call fit(refslt_natom, meancoord, bfcoord, stmass, outcoord)
                cencoord(:,:,i) = outcoord(1:3,1:refslt_natom)
        enddo

!収束loop*3
        do i = 1, 3
           do j = 1, refslt_natom
                sumcoord(1,j) = sum(cencoord(1,j,:))
                sumcoord(2,j) = sum(cencoord(2,j,:))
                sumcoord(3,j) = sum(cencoord(3,j,:))
           enddo
                meancoord(:,:) = sumcoord(1:3,1:refslt_natom) / frame
                deallocate( outcoord, cencoord )
                allocate( outcoord(3, refslt_natom), cencoord(3, refslt_natom, frame) )
                call fit(refslt_natom, meancoord, cencoord, stmass, outcoord)
                cencoord(:,:,i) = outcoord(1:3,1:refslt_natom)
        enddo
end subroutine calcmean

!Variance-covariance matrix
subroutine varcor( refslt_natom, xd, varm, corm )
        integer,intent(in) :: refslt_natom
        double precision,intent(in) :: xd(3, refslt_natom)
        double precision,intent(out) :: varm(3,3), corm(3,3)
        integer :: i, j, k
        double precision :: xm(3), wmat(3, refslt_natom)

        do j = 1, refslt_natom
           do i = 1, 3
              wmat(i,j) = xd(i,j)
           enddo
        end do
        do i = 1, 3
           xm(i) = 0.0D0
           do j = 1, refslt_natom
              xm(i) = xm(i) + wmat(i,j)
           enddo
           xm(i) = xm(i) / dble(refslt_natom)
        end do
        do i = 1, 3
           do j = 1, 3
              varm(i,j) = 0.0D0
              do k = 1, refslt_natom
                 varm(i,j) = varm(i,j) + (wmat(i,k) - xm(i))*(wmat(j,k) - xm(j))
              enddo
              varm(i,j) = varm(i,j) / dble(refslt_natom - 1)
           enddo
        enddo
        do i = 1, 3
           do j = 1, 3
              corm(i,j) = varm(i,j) / sqrt(varm(i,i)) / sqrt(varm(j,j))
           enddo
        enddo
end subroutine varcor

!Jacobi method for eigen value analysis
subroutine jacobi(array,ev,vec)
        double precision,intent(inout) :: array(3,3)
        double precision,intent(out) :: ev(3), vec(3,3)
        integer :: i, j, k, im
        double precision :: s, c, aoff, eps, t, aa, aj, ak, vj, vk, vmax, bb

        do i = 1, 3
           do j = 1, 3
              vec(i,j) = 0.0
           enddo
           vec(i,i) = 1.0
        enddo

        aoff = 0.0
        do j = 1, 2
           do k = j+1, 3
              aoff = aoff + array(j,k)*array(j,k)
           enddo
        enddo
        eps = aoff*2.0
        do j = 1, 3
           eps = eps + array(j,j)*array(j,j)
        enddo
        eps = eps*0.00000000005

        do while (aoff>eps)
            do j = 1, 2
               do k = j+1, 3
                  if (array(j,k) == 0.0) then
                        t = 0.0
                  else
                        aa = (array(k,k) - array(j,j)) / (2.0*array(j,k))
                  endif
                  if (aa>=0) then
                        t = 1.0 / (aa + sqrt(1.0 + aa*aa))
                  else
                        t = 1.0 / (aa - sqrt(1.0 + aa*aa))
                  endif
                  c = 1.0 / sqrt(1.0+t*t)
                  s = t*c
                  array(j,j) = array(j,j) - t*array(j,k)
                  array(k,k) = array(k,k) + t*array(j,k)
                  array(j,k) = 0.0
                  do i = 1, j-1
                        aj = array(i,j)
                        ak = array(i,k)
                        array(i,j) = aj*c - ak*s
                        array(i,k) = aj*s + ak*c
                  enddo
                  do i = j+1, k-1
                        aj = array(j,i)
                        ak = array(i,k)
                        array(j,i) = aj*c - ak*s
                        array(i,k) = aj*s + ak*c
                  enddo
                  do i = k+1, 3
                        aj = array(j,i)
                        ak = array(k,i)
                        array(j,i) = aj*c - ak*s
                        array(k,i) = aj*s + ak*c
                  enddo
                  do i = 1, 3
                        vj = vec(i,j)
                        vk = vec(i,k)
                        vec(i,j) = vj*c - vk*s
                        vec(i,k) = vj*s + vk*c
                  enddo
              enddo
          enddo
          aoff = 0.0
          do j = 1, 2
              do k = j+1, 3
                 aoff = aoff + array(j,k)*array(j,k)
              enddo
          enddo
        enddo

        do k = 1, 3
           ev(k) = array(k,k)
        enddo

        do k = 1, 2
           vmax = -1.0E+30
           do i = k, 3
              if (ev(i)>=vmax) then
                  im = i
                  vmax = ev(i)
              endif
           enddo
           aa = ev(k)
           ev(k) = ev(im)
           ev(im) = aa
           do i = 1, 3
              bb = vec(i,k)
              vec(i,k) = vec(i,im)
              vec(i,im) = bb
           enddo
        enddo
end subroutine jacobi

!Calculate entropy
!subroutine calcentropy
!        implicit none
!        
!z = sum()
!end subroutine calcentropy

!Output to xtcfile
subroutine outputfile(bffilename, refslt_natom, cencoord)
        implicit none
        integer :: i, refslt_natom
        double precision, allocatable, dimension(:,:,:) :: cencoord
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
                call xtc_out % write(refslt_natom, xtcbf % step, xtcbf % time,&
xtcbf % box, pos, xtcbf % prec)
                deallocate( pos )
        enddo
        call xtc_out % close
        call xtcbf % close
end subroutine outputfile

!Calculate center of mass
subroutine center_mass(refslt_natom, frame, refcoord, stmass, center)
        implicit none
        integer, intent(in) :: refslt_natom
        double precision, intent(in) :: refcoord(3, refslt_natom, frame), stmass(refslt_natom)
        double precision, intent(out) :: center(3)
        double precision :: summass
        integer :: i, j, frame

        summass = sum(stmass)
        if(summass == 0) write(*,*), "error: sum of mass = 0"
        do j = 1, frame
           do i = 1, 3
              center(i) = dot_product(refcoord(i,:,j), stmass(:)) / summass
           end do
        enddo
end subroutine center_mass

end program bestfit_ty
