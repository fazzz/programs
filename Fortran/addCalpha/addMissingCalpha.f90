! This program adds missing residues (Calpha) to N-Terminal and C-Terminal.
program addMissingCalpha
  use pdb

  implicit none

  integer i, j, k
  
  real,dimension(3) :: vec_C2C_NT, vec_C2C_CT ! Vector of Calpha (N-1) to Calpha (N)
  real,allocatable,dimension(:,:) :: r_Ca_NT, r_Ca_NT ! Cartesian coordinate of Calpha (adding)
  type (PDB_atom), allocatable, dimension(:) :: PA, PAadded

  integer IndexIniCaNT, IndexFinCaNT  ! # of Calpha to add
  integer IndexIniCaCT, IndexFinCaCT  ! # of Calpha to add

  integer NaddedCT, NaddedNT ! # of res.

  integer N, Nadded
  
  character inputPDBfilename * 50, outputPDBfilename * 50

  NCa = 1
  
  NAMELIST/addMiss/IndexIniCaNT, IndexFinCaNT, IndexIniCaCT, IndexFinCaCT, NaddedCT, NaddedNT
  open(16,file='parameters_aMC')
  read(16,addMiss)
  close(16)

  argc=iargc()
  call getarg(0,programname)
  if (argc < 2) then
     call usage(programname)
  else
     call getarg(1,inputPDBfilename)
     call getarg(2,outputPDBfilename)
  end if

  open(17,file=pdb_filename,status='old')
  call read_PDB (17, PA, N)
  close(17)

  vec_C2C_NT(1) = PA(IndexFinCa)%x - PA(IndexIniCa)%x
  vec_C2C_NT(2) = PA(IndexFinCa)%y - PA(IndexIniCa)%y
  vec_C2C_NT(3) = PA(IndexFinCa)%z - PA(IndexIniCa)%z

  vec_C2C_NT(1) = PA(IndexFinCa)%x - PA(IndexIniCa)%x
  vec_C2C_NT(2) = PA(IndexFinCa)%y - PA(IndexIniCa)%y
  vec_C2C_NT(3) = PA(IndexFinCa)%z - PA(IndexIniCa)%z

  Nadded = N + NaddedNT + NaddedCT + 2

  allocate(PAadded(Nadded))
  
  j = 1
  do i =1,N,1
     if (i .ne. Na) then
        PAadded(j) = PA(i)
        j = j + 1
     else
        PAadded(j) = PA(IndexFinCa)%x + k * vec_C2C(1)
        PAadded(j) = PA(IndexFinCa)%x + k * vec_C2C(2)
        PAadded(j) = PA(IndexFinCa)%x + k * vec_C2C(3)
        j = j + 1
     end if
        
  end do

end program addMissingCalpha
