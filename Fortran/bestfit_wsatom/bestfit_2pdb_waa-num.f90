! This program perform bestfit two pdb strcturues
! based on all atoms
! You can specify the atom range to be fited.

program bestfit_2pdb
  use bestfit
  use gettop
  use pdb

  implicit none

  integer i,j

  integer N, Nref, Nbsf ! num of atoms
  integer iniatom, finatom

  type (PDB_atom), allocatable, dimension(:) :: PAref, PA !, PAout

!  real, allocatable, dimension(:,:) :: refcoord, coord, outcoord
!  real, allocatable, dimension(:) :: massHA

  double precision, allocatable, dimension(:,:) :: refcoord, coordbsf, coord, outcoord
  double precision, allocatable, dimension(:) :: massbsf

  double precision rmsd, len
  
  integer,allocatable,dimension(:) :: nr, resid, cgnr

  double precision,allocatable,dimension(:) :: charg, mass
!  character,allocatable,dimension(:,:) :: atype, res, atom
  character(5) :: atype, res, atom

  character(50) top_filename, pdb_reffilename, pdb_infilename, pdb_outfilename

  character(2) atm_name_HEAD

  character(10) :: argv, programname
  integer :: argc, iargc

!  integer, allocatable, dimension(:) :: serial
!  character, allocatable, dimension(:,:) :: atm_name
!  character, allocatable, dimension(:) :: alte_l_i
!  character, allocatable, dimension(:,:) :: res_name
!  character, allocatable, dimension (:) :: chain_i
!  integer, allocatable,dimension(:) :: res_seq_num
!  character, allocatable, dimension(:) :: code_for_insertion_res
  
  NAMELIST/aarange/iniatom, finatom
  open(20,file='parameters_bsf')
  read(20,aarange)
  close(20)

  Nbsf = finatom - iniatom
  Nbsf = Nbsf + 1

  if (Nbsf <= 0) then
     write(*,'("ERROR: parameter iniatom must be larger than parameter finatom.")')
     stop
  end if
  
  argc=iargc()
  call getarg(0,programname)
  if (argc < 3) then
     call usage(programname)
  else
     call getarg(1,top_filename)
     call getarg(2,pdb_reffilename)
     call getarg(3,pdb_infilename)
!     call getarg(4,pdb_outfilename)
  end if

  open(17,file=pdb_reffilename,status='old')
  call read_PDB (17, PAref, Nref)
  close(17)

!  NHAref = 0
!  do i=1,Nref
!     atm_name_HEAD = PAref(i)%atm_name
!     if ( atm_name_HEAD  .ne. "H" ) then
!        NHAref = NHAref + 1
!     end if
!  end do

!  write(*,*), NHAref

  open(19,file=pdb_infilename,status='old')
  call read_PDB (19, PA, N)
  close(19)

!  allocate(PAout(N))
!  PAout = PA
  
!  NHA = 0
!  do i=1,N
!     atm_name_HEAD = PA(i)%atm_name
!     if ( atm_name_HEAD .ne. "H" ) then
!        NHA = NHA + 1
!     end if
!  end do

!  write(*,*), NHA

!  if ( NHAref .ne. NHA ) then
!    write (*,*), "# of atoms of two PDBs is not equal."
!    stop
!  end if

!  write(*,*), N

  allocate(mass(N))
  allocate(nr(N))
  allocate(resid(N))
  allocate(cgnr(N))
  allocate(charg(N))
!  allocate(atype(5,N))
!  allocate(res(5,N))
!  allocate(atom(5,N))

  open(19,file=top_filename,status='old')
  call get_top (19, N, nr, atype, resid, res, atom, cgnr, charg, mass)
  close(19)

!  write(*,*), nr(10)
!  write(*,*), atom(:,10)
!  write(*,*), mass(:)

!  allocate(massHA(NHA))

!  allocate(refcoord(3,NHA))
!  allocate(coord(3,NHA))
!  allocate(outcoord(3,NHA))

  allocate(massbsf(Nbsf))

  allocate(refcoord(3,Nbsf))
  allocate(coordbsf(3,Nbsf))
  allocate(coord(3,N))
  allocate(outcoord(3,N))

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! j = 1                                 !
  ! do i=1,N                              !
  !    atm_name_HEAD = PA(i)%atm_name     !
  !    if ( atm_name_HEAD .ne. "H" ) then !
  !       massHA(j)=mass(i)               !
  !       refcoord(1,j)=PAref(i)%x        !
  !       refcoord(2,j)=PAref(i)%y        !
  !       refcoord(3,j)=PAref(i)%z        !
  !       coord(1,j)=PA(i)%x              !
  !       coord(2,j)=PA(i)%y              !
  !       coord(3,j)=PA(i)%z              !
  !       j = j + 1                       !
  !    end if                             !
  ! end do                                !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  massbsf(:)=mass(iniatom:finatom)
  do i=iniatom,finatom,1 
     refcoord(1,i-iniatom+1)=PAref(i)%x
     refcoord(2,i-iniatom+1)=PAref(i)%y
     refcoord(3,i-iniatom+1)=PAref(i)%z     

     coordbsf(1,i-iniatom+1)=PA(i)%x
     coordbsf(2,i-iniatom+1)=PA(i)%y
     coordbsf(3,i-iniatom+1)=PA(i)%z     
  end do
  coord(1,:)=PA%x
  coord(2,:)=PA%y
  coord(3,:)=PA%z

!  write(*,*), massbsf(:,10)

!  write(*,*), refcoord(:,1)
!  write(*,*), refcoord(:,2)
!  write(*,*), refcoord(:,3)

!  write(*,*), coordbsf(:,1)
!  write(*,*), coordbsf(:,2)
!  write(*,*), coordbsf(:,3)

!  write(*,*), refcoord(:,1000)

!  write(*,*), coord(:,1000)

!  write(*,*), coord(:,1)
!  write(*,*), coord(:,2)
!  write(*,*), coord(:,3)

!  call fit(NHA, refcoord, coord, massHA, outcoord)
  call fit_a_rotate_b(Nbsf, refcoord, coordbsf, massbsf, N, coord, outcoord)
  
  rmsd = 0.0d0
  do i=iniatom,finatom,1
     len = 0.0d0
     do j=1,3,1
        len = len + (refcoord(j,i-iniatom+1) - outcoord(j,i))**2
     end do
     rmsd = rmsd + len
  end do
  rmsd = rmsd / Nbsf
  rmsd = sqrt(rmsd)
  
!  write(*,'("The Heavy-Atom RMSD(A) from ref. to fited. is")'
!  open(21,file="bft.log",position='append')
!  write(*,'("The All-Atom RMSD(A) from ref. to fited. is F8.3")'),rmsd
!  write(*,'(F8.3)'),rmsd
!  close(21)
  
!  write(*,*), outcoord(:,1)
!  write(*,*), outcoord(:,2)
!  write(*,*), outcoord(:,3)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! j = 0                                 !
  ! do i=1,N                              !
  !    atm_name_HEAD = PA(i)%atm_name     !
  !    if ( atm_name_HEAD .ne. "H" ) then !
  !       PA(i)%x=outcoord(1,j)           !
  !       PA(i)%y=outcoord(2,j)           !
  !       PA(i)%z=outcoord(3,j)           !
  !       j = j +1                        !
  !    end if                             !
  ! end do                                !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
!  PAout%x=outcoord(1,:)
!  PAout%y=outcoord(2,:)
!  PAout%z=outcoord(3,:)

!  do i=1,N                                                       !
!     write(*,'("ATOM  ",I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3)'),& !
!          PAref(i)%serial,PAref(i)%atm_name,&                          !
!          PAref(i)%alte_l_i,PAref(i)%res_name,&                        !
!          PAref(i)%chain_i,mod(PAref(i)%res_seq_num,9999),&            !
!          PAref(i)%code_for_insertion_res,&                         !
!          PAref(i)%x,PAref(i)%y,PAref(i)%z                                !
!  enddo                                                          !
  
 do i=1,N
    write(*,'("ATOM  ",I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3)'),&
          PA(i)%serial,PA(i)%atm_name,&
          PA(i)%alte_l_i,PA(i)%res_name,&
          PA(i)%chain_i,mod(PA(i)%res_seq_num,9999),&
          PA(i)%code_for_insertion_res,&
          outcoord(1,i),outcoord(2,i),outcoord(3,i)
  enddo

!  open(21,file=pdb_outfilename,position='append')
!  call write_PDB_woH (21, PA, N)
!  call write_PDB (21, PA, N)
!  close(21)

  deallocate(PAref)
  deallocate(PA)

  deallocate(mass)

  deallocate(nr)
!  deallocate(atype) 
!  deallocate(resid) 
!  deallocate(res) 
!  deallocate(atom) 
  deallocate(cgnr) 
  deallocate(charg) 

  deallocate(refcoord)
  deallocate(coord)
  deallocate(outcoord)

contains
  subroutine usage(programname)
    implicit none
    character(10) :: programname
!    print *,"Usage: ",programname,"topfilename PDB1(reffile) PDB2(input) PDB3(output)"
    print *,"Usage: ",programname,"topfilename PDB1(reffile) PDB2(input)"
    stop
  end subroutine usage

end program bestfit_2pdb
