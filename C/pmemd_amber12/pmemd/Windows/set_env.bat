::Script to set environment variables and start Visual Studio
::Choose cuda and MPI settings (ignore if neither apply)

::Choose cuda version here
::  use V4_0 or V3_2

set CUDA_PATH=%CUDA_PATH_V4_0%

::Derived variables
	set CUDA_INC_PATH=%CUDA_PATH%\include
	set CUDA_BIN_PATH=%CUDA_PATH%\bin

::Choose mpi implementation here
:: use msmpi or mpich2
goto msmpi

:msmpi
	echo "Using Microsoft MPI"
	set MPI_INC32=%MSMPI_INC%
	set MPI_INC64=%MSMPI_INC%
	set MPI_LIB32=%MSMPI_LIB32%
	set MPI_LIB64=%MSMPI_LIB64%

	:: These are the libraries we need to link to
	set MPI_LIBS32="$(MPI_LIB32)\msmpi.lib" "$(MPI_LIB32)\msmpifmc.lib"
	set MPI_LIBS64="$(MPI_LIB64)\msmpi.lib" "$(MPI_LIB64)\msmpifmc.lib"
	goto continue

:mpich2
	echo "Using Mpich2"
	set MPI_INC32=C:\Program Files (x86)\MPICH2\include
	set MPI_INC64=C:\Program Files\MPICH2\include
	set MPI_LIB32=C:\Program Files (x86)\MPICH2\lib
	set MPI_LIB64=C:\Program Files\MPICH2\lib
	:: These are the libraries we need to link to
	set MPI_LIBS32="$(MPI_LIB32)\mpi.lib" "$(MPI_LIB32)\fmpich2.lib"
	set MPI_LIBS64="$(MPI_LIB64)\mpi.lib" "$(MPI_LIB64)\fmpich2.lib"
	goto continue

:continue
echo 'Done!'
echo 'Starting Visual Studio...'
".\pmemd.sln"