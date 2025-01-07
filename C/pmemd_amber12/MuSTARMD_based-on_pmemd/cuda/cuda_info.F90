#include "../copyright.i"

!*******************************************************************************
!
! cuda_info
!
! Description: 
!
! By Ross C. Walker (SDSC, 2009)
!
! Routines to print info regarding CUDA based GPU Calculations.
!              
!*******************************************************************************

!*******************************************************************************
!
! Subroutine: gpu_write_cuda_info
!
! Description:
!
! Writes information about the GPU(s) being used for the calculation
! to the output file.
!              
!*******************************************************************************

subroutine gpu_write_cuda_info(iout_unit, mytaskid, numtasks)

  implicit none

! Formal arguments:

  integer, intent(in)   :: iout_unit !Unit to write info to.
  integer, intent(in)   :: mytaskid  !MPI task ID.
  integer, intent(in)   :: numtasks  !Number of MPI tasks.

! Local variables:
  
  integer :: gpu_dev_count,gpu_dev_id,gpu_dev_mem
  integer :: gpu_num_multi, gpu_num_proc, name_len, i
  double precision :: gpu_core_freq
  character(len=80) :: gpu_name

  write(iout_unit,'(a)')        '|--------------------- INFORMATION ----------------------'
  write(iout_unit,'(a)')        '| GPU (CUDA) Version of PMEMD in use: NVIDIA GPU IN USE.'
  write(iout_unit,'(a)')        '|                     Version 12.0'
  write(iout_unit,'(a)')        '| '
  write(iout_unit,'(a)')        '|                      03/19/2012'
  write(iout_unit,'(a)')        '| '
  write(iout_unit,'(a)')        '| Implementation by:'
  write(iout_unit,'(a)')        '|                    Ross C. Walker     (SDSC)'
  write(iout_unit,'(a)')        '|                    Scott Le Grand     (nVIDIA)'
  write(iout_unit,'(a)')        '|                    Duncan Poole       (nVIDIA)'
  write(iout_unit,'(a)')        '| '
  write(iout_unit,'(a)')        '| CAUTION: The CUDA code is currently experimental.'
  write(iout_unit,'(a)')        '|          You use it at your own risk. Be sure to'
  write(iout_unit,'(a)')        '|          check ALL results carefully.'
  write(iout_unit,'(a)')        '| '
  write(iout_unit,'(a)')        '| Precision model in use:'
#ifdef use_SPSP
  write(iout_unit,'(a)')        '|      [SPSP] - All Single Precision (ex. Shake).'
#elif use_DPDP
  write(iout_unit,'(a)')        '|      [DPDP] - All Double Precision.'
#else
  write(iout_unit,'(a)')        '|      [SPDP] - Hybrid Single/Double Precision (Default).'
#endif
  write(iout_unit,'(a)')        '| '
  write(iout_unit,'(a)')        '|--------------------------------------------------------'
  write(iout_unit,'(a)')        ' '

! Get info about GPU(s) in use.
  write(iout_unit,'(a)')        '|------------------- GPU DEVICE INFO --------------------'
#ifdef MPI  
  do i = 0, numtasks - 1
    if (i .eq. 0) then
#endif
      call gpu_get_device_info(gpu_dev_count,gpu_dev_id,gpu_dev_mem,gpu_num_proc,gpu_core_freq,gpu_name,name_len)
#ifdef MPI
    else
      call gpu_get_slave_device_info(i, gpu_dev_count,gpu_dev_id,gpu_dev_mem,gpu_num_proc,gpu_core_freq,gpu_name,name_len)
    end if
#endif
    write(iout_unit,'(a)')      '|'
#ifdef MPI
    write(iout_unit,'(a,i6)')   '|                         Task ID: ',i
#endif
    write(iout_unit,'(a,i6)')   '|   CUDA Capable Devices Detected: ',gpu_dev_count
    write(iout_unit,'(a,i6)')   '|           CUDA Device ID in use: ',gpu_dev_id
    write(iout_unit,'(a,a)')    '|                CUDA Device Name: ',gpu_name(1:name_len)
    write(iout_unit,'(a,i6,a)') '|     CUDA Device Global Mem Size: ',gpu_dev_mem,' MB'
    write(iout_unit,'(a,i6)')   '| CUDA Device Num Multiprocessors: ',gpu_num_proc
    write(iout_unit,'(a,f6.2,a)') '|           CUDA Device Core Freq: ',gpu_core_freq,' GHz'
    write(iout_unit,'(a)')      '|'
#ifdef MPI
  end do
#endif  
  write(iout_unit,'(a)')        '|--------------------------------------------------------'
  write(iout_unit,'(a)')   ' '

  return

end subroutine gpu_write_cuda_info

subroutine gpu_write_memory_info(iout_unit)

  implicit none

! Formal arguments:

  integer, intent(in)   :: iout_unit !Unit to write info to.

! Local variables:
  integer :: gpumemory, cpumemory
  
  call gpu_get_memory_info(gpumemory, cpumemory)

  write(iout_unit, '(a)')    '| GPU memory information:'
  write(iout_unit,'(a,i9)')  '| KB of GPU memory in use: ', gpumemory
  write(iout_unit,'(a,i9/)') '| KB of CPU memory in use: ', cpumemory

  return

end subroutine gpu_write_memory_info

