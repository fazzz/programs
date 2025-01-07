module pdb
  implicit none

  type PDB_atom
!     character(6) nameRecord
     integer serial
     character(4) atm_name
     character alte_l_i
     character(3) res_name
     character chain_i

     integer res_seq_num
     character code_for_insertion_res
     double precision x, y, z
!     real occupancy,temperature_factor
!     character(2) seq_i
!     character(2) element
  end type PDB_atom

! write_PDB <- To write output coordinates in Protein Data Bank (PDB) format.
! read_PDB  <- To read inout coordinate in PDB format.

contains

  subroutine write_PDB (unitnum, PA, N)

    implicit none
    type (PDB_atom),allocatable,dimension(:),intent(in) :: PA
    integer, intent(in) :: unitnum
    integer, intent(in) :: N

     character(6) nameRecord
!     character alte_l_i
!     character chain_i

!     integer res_seq_num
!     character code_for_insertion_res
!     real occupancy,temperature_factor
!     character(2) seq_i
!     character(2) element

     character(1) atm_name_HEAD

    integer i

    nameRecord="ATOM  "
    
    write(unitnum,'("MODEL")')
    do i = 1,N
       ! write(unitnum,'(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,10X,A2,A2)'),&
       !      PA(i)%nameRecord,PA(i)%serial,PA(i)%atm_name,PA(i)%alte_l_i,PA(i)%res_name,&
       !      PA(i)%chain_i,mod(PA(i)%res_seq_num,9999),PA(i)%code_for_insertion_res,&
       !      PA(i)%x,PA(i)%y,PA(i)%z,&
       !      PA(i)%occupancy,PA(i)%temperature_factor,&
       !      PA(i)%seq_i,PA(i)%element

!     atm_name_HEAD = PA(i)%atm_name
!     if ( atm_name_HEAD  .ne. "H" ) then
       write(unitnum,'(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3)'),&
            nameRecord,&
            PA(i)%serial,PA(i)%atm_name,&
            PA(i)%alte_l_i,PA(i)%res_name,&
            PA(i)%chain_i,mod(PA(i)%res_seq_num,9999),&
            PA(i)%code_for_insertion_res,&
            PA(i)%x,PA(i)%y,PA(i)%z
!     end if

    end do
    write(unitnum,'("ENDML")')
    
  end subroutine write_PDB

  subroutine read_PDB (unitnum, PA, N)
    implicit none
    integer, intent(in) :: unitnum
    type (PDB_atom),allocatable,dimension(:),intent(inout) :: PA
    integer, intent(out) :: N
    
    type (PDB_atom),allocatable,dimension(:) :: A

    character(6) nameRecord
!    character alte_l_i
!    character chain_i
    
!    integer res_seq_num
!    character code_for_insertion_res
!    real occupancy,temperature_factor
!    character(2) seq_i
!    character(2) element

    integer i, j

    nameRecord="ATOM  "

    i = 1
    allocate(PA(i))
    do
       ! read(unitnum,'(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,10X,A2,A2)',end=999),&
       !      PA(i)%nameRecord,&
       !      PA(i)%serial,PA(i)%atm_name,PA(i)%alte_l_i,PA(i)%res_name,PA(i)%chain_i,&
       !      PA(i)%res_seq_num,PA(i)%code_for_insertion_res,PA(i)%x,PA(i)%y,PA(i)%z,&
       !      PA(i)%occupancy,PA(i)%temperature_factor,&
       !      PA(i)%seq_i,PA(i)%element

       read(unitnum,'(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,10X,A2,A2)',end=999),&
            nameRecord,&
            PA(i)%serial,PA(i)%atm_name,&
            PA(i)%alte_l_i,PA(i)%res_name,&
            PA(i)%chain_i,&
            PA(i)%res_seq_num,&
            PA(i)%code_for_insertion_res,&
            PA(i)%x,PA(i)%y,PA(i)%z!,&
!            occupancy,temperature_factor,&
!            seq_i,element

       allocate(A(i))
       
!       print *,i
!       print *,"YES"
       
       do j = 1,i
!          A(j)%nameRecord = PA(j)%nameRecord
          A(j)%serial = PA(j)%serial
          A(j)%atm_name = PA(j)%atm_name
          A(j)%alte_l_i = PA(j)%alte_l_i
          A(j)%res_name = PA(j)%res_name
          A(j)%chain_i = PA(j)%chain_i
          A(j)%res_seq_num = PA(j)%res_seq_num
          A(j)%code_for_insertion_res = PA(j)%code_for_insertion_res
          A(j)%x = PA(j)%x
          A(j)%y = PA(j)%y
          A(j)%z = PA(j)%z
!          A(j)%occupancy = PA(j)%occupancy
!          A(j)%temperature_factor = PA(j)%temperature_factor
!          A(j)%seq_i = PA(j)%seq_i
!          A(j)%element = PA(j)%element
       end do

!       A = PA

       deallocate(PA)
       i = i + 1
       allocate(PA(i))

!       print *,i
!       print *,"YES YES"

       do j = 1,i-1
!          PA(j)%nameRecord = A(j)%nameRecord
          PA(j)%serial = A(j)%serial
          PA(j)%atm_name = A(j)%atm_name
          PA(j)%alte_l_i = A(j)%alte_l_i
          PA(j)%res_name = A(j)%res_name
          PA(j)%chain_i = A(j)%chain_i
          PA(j)%res_seq_num = A(j)%res_seq_num
          PA(j)%code_for_insertion_res = A(j)%code_for_insertion_res
          PA(j)%x = A(j)%x
          PA(j)%y = A(j)%y
          PA(j)%z = A(j)%z
!          PA(j)%occupancy = A(j)%occupancy
!          PA(j)%temperature_factor = A(j)%temperature_factor
!          PA(j)%seq_i = A(j)%seq_i
!          PA(j)%element = A(j)%element
       end do

!       print *,"YES YES YES"

       !       PA = A
       deallocate(A)
    end do
999 N = i - 1

  end subroutine read_PDB

end module PDB
