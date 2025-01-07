! This program perform bestfit two pdb strcturues

program bestfit_2pdbwaa
  use bestfit
  use ndx

  implicit none

  integer i,j,k,l

  integer N, Nref, Nbsf
  character(6) nameRecord
  integer serial
  character(4) atm_name
  character alte_l_i
  character(3) res_name
  character chain_i
  integer res_seq_num
  character code_for_insertion_res
  double precision x, y, z

  double precision, allocatable, dimension(:,:) :: refcoord, coordbsf, coord, outcoord, A
  double precision, allocatable, dimension(:) :: massbsf

  character(1) atm_name_HEAD
  
  integer Nndx_ref, Nndx_in
  integer,allocatable, dimension(:) :: atom_list_ref, atom_list_in
  
  character(400) ndx_reffilename, ndx_infilename
  character(400) pdb_reffilename, pdb_infilename

  character(10) :: argv, programname
  integer :: argc, iargc

  argc=iargc()
  call getarg(0,programname)
  if (argc < 4) then
     call usage(programname)
  else
     call getarg(1,pdb_reffilename)
     call getarg(2,pdb_infilename)
     call getarg(3,ndx_reffilename)
     call getarg(4,ndx_infilename)
  end if

  allocate(atom_list_ref(1))
  open(20,file=ndx_reffilename,status='old')
  call read_ndx(20, Nndx_ref, atom_list_ref)
  close(20)

  allocate(atom_list_in(1))
  open(21,file=ndx_infilename,status='old')
  call read_ndx(21, Nndx_in, atom_list_in)
  close(21)

  if (Nndx_ref .ne. Nndx_in) then
     write(*,'("ERROR: Nndx_ref must equal to Nndx_in.")')
     stop
  end if

  allocate(massbsf(Nndx_in))
  allocate(refcoord(3,Nndx_in))
  allocate(coordbsf(3,Nndx_in))

  open(17,file=pdb_reffilename,status='old')
  i=1
  j=1
  do
     read(17,'(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,10X,A2,A2)',end=991),&
          nameRecord,&
          serial,atm_name,&
          alte_l_i,res_name,&
          chain_i,&
          res_seq_num,&
          code_for_insertion_res,&
          x,y,z
     if ( i .eq. atom_list_ref(j) ) then
        refcoord(1,j)=x
        refcoord(2,j)=y
        refcoord(3,j)=z

        atm_name_HEAD = atm_name(2:2)
        if ( atm_name_HEAD .eq. "C" ) then
           massbsf(j)=12.0
        else if ( atm_name_HEAD .eq. "O" ) then
           massbsf(j)=16.0
        else if ( atm_name_HEAD .eq. "N" ) then
           massbsf(j)=14.0
        else
           massbsf(j)=1.0
        end if

        j = j + 1
     end if
     i = i + 1
  end do
  991 N = i - 1
  close(17)

!  open(19,file=pdb_infilename,status='old')
!  i=0
!  do
!     read(19,*,end=999)
!     i = i + 1
!  end do
!999 N = i - 1
!  close(19)

!  write(*,*) N

  allocate(coord(3,1))
!  allocate(coord(3,N))

!  write(*,*) "YES"

  open(19,file=pdb_infilename,status='old')
  i=0
  j=1
  do ! i=1,N,1
     read(19,'(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,10X,A2,A2)',end=999),&
          nameRecord,&
          serial,atm_name,&
          alte_l_i,res_name,&
          chain_i,&
          res_seq_num,&
          code_for_insertion_res,&
          x,y,z

!     read(19,'(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,10X,A2,A2)'),&
!          nameRecord,&
!          serial,atm_name,&
!          alte_l_i,res_name,&
!          chain_i,&
!          res_seq_num,&
!          code_for_insertion_res,&
!          coord(1,i), coord(2,i), coord(3,i)

     if ( i .gt. 0 ) then
        allocate(A(3,i))
        A=coord

        deallocate(coord)
        allocate(coord(3,i+1))
        do k=1,i,1
           do l=1,3,1
              coord(l,k)=A(l,k)
           end do
        end do

        deallocate(A)
     end if

     coord(1,i+1)=x
     coord(2,i+1)=y
     coord(3,i+1)=z

     if ( (i .eq. atom_list_in(j)) .and. (j .le. Nndx_in)) then
!     if ( i+1 .eq. atom_list_in(j) ) then
        coordbsf(1,j)=x
        coordbsf(2,j)=y
        coordbsf(3,j)=z

!        write(*,*) j
        j = j + 1
     end if
!     write(*,*) i
     
     i = i + 1
  end do
  999 N = i - 1
  close(19)

!  write(*,*) "YES3"

  allocate(outcoord(3,N))

  deallocate(atom_list_ref)
  deallocate(atom_list_in)

!  write(*,*) "YES4"

  call fit_a_rotate_b(Nndx_in, refcoord, coordbsf, massbsf, N, coord, outcoord)

!  write(*,*) "YES5"

  deallocate(massbsf)
  deallocate(refcoord)
  deallocate(coordbsf)
  deallocate(coord)

!  write(*,*) "YES6"

  open(25,file=pdb_infilename,status='old')
  do i = 1,N,1
     read(25,'(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3)'),&
          nameRecord,&
          serial,atm_name,&
          alte_l_i,res_name,&
          chain_i,&
          res_seq_num,&
          code_for_insertion_res,&
          x,y,z

     write(*,'(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3)'),&
          nameRecord,&
          serial,atm_name,&
          alte_l_i,res_name,&
          chain_i,res_seq_num,&
          code_for_insertion_res,&
          outcoord(1,i),outcoord(2,i),outcoord(3,i)
  enddo

  close(25)
  deallocate(outcoord)

contains
  subroutine usage(programname)
    implicit none
    character(10) :: programname
    print *,"Usage: ",programname,"PDB1(reffile) PDB2(input) NDX1(reffilename) NDX2(input)"
    stop
  end subroutine usage

end program bestfit_2pdbwaa












