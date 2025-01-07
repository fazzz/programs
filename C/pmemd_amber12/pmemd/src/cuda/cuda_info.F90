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

subroutine gpu_write_cuda_info(iout_unit, mytaskid, numtasks, using_pme_potential,iamd)

  implicit none

! Formal arguments:

  integer, intent(in)   :: iout_unit !Unit to write info to.
  integer, intent(in)   :: mytaskid  !MPI task ID.
  integer, intent(in)   :: numtasks  !Number of MPI tasks.
  integer, intent(in)   :: using_pme_potential  
  integer, intent(in)   :: iamd  

! Local variables:
  
  integer :: gpu_dev_count,gpu_dev_id,gpu_dev_mem
  integer :: gpu_num_multi, gpu_num_proc, name_len, i
  double precision :: gpu_core_freq
  character(len=80) :: gpu_name

  write(iout_unit,'(a)')        '|--------------------- INFORMATION ----------------------'
  write(iout_unit,'(a)')        '| GPU (CUDA) Version of PMEMD in use: NVIDIA GPU IN USE.'
  write(iout_unit,'(a)')        '|                     Version 12.2'
  write(iout_unit,'(a)')        '| '
  write(iout_unit,'(a)')        '|                      01/10/2013'
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
#ifdef use_SPFP
  write(iout_unit,'(a)')        '|      [SPFP] - Mixed Single/Double/Fixed Point Precision.'
  write(iout_unit,'(a)')        '|               (Default)'
#elif use_DPDP
  write(iout_unit,'(a)')        '|      [DPDP] - All Double Precision.'
#elif use_SPDP
  write(iout_unit,'(a)')        '|      [SPDP] - Hybrid Single/Double Precision.'
#else
  write(iout_unit,'(a)')        '|      ERROR ERROR ERROR - Unknown GPU Precision Model'
  call mexit(iout_unit, 1)
#endif
  write(iout_unit,'(a)')        '| '
  write(iout_unit,'(a)')        '|--------------------------------------------------------'
  write(iout_unit,'(a)')        ' '

!Write citation information for the GPU code to mdout.
  call write_cuda_citation(iout_unit, using_pme_potential,iamd)

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

subroutine write_cuda_citation(iout_unit, using_pme_potential,iamd)

!Write the citation information for the GPU code.
  
  integer, intent(in)   :: using_pme_potential
  integer, intent(in)   :: iamd

  write(iout_unit,'(a)')        '|----------------- CITATION INFORMATION -----------------'
  write(iout_unit,'(a)')        '|'
  write(iout_unit,'(a)')        '|    When publishing work that utilized the CUDA version'
  write(iout_unit,'(a)')        '|    of AMBER, please cite the following in addition to'
  write(iout_unit,'(a)')        '|    the regular AMBER citations:' 
  write(iout_unit,'(a)')        '|' 
if (using_pme_potential/=0) then
  write(iout_unit,'(a)')        '|  - Romelia Salomon-Ferrer; Andreas W. Goetz; Duncan'
  write(iout_unit,'(a)')        '|    Poole; Scott Le Grand; Ross C. Walker "Routine'
  write(iout_unit,'(a)')        '|    microsecond molecular dynamics simulations with'
  write(iout_unit,'(a)')        '|    AMBER - Part II: Particle Mesh Ewald", J. Chem.' 
  write(iout_unit,'(a)')        '|    Theory Comput., 2012, (In prep).'
  write(iout_unit,'(a)')        '|'
endif
  write(iout_unit,'(a)')        '|  - Andreas W. Goetz; Mark J. Williamson; Dong Xu;'
  write(iout_unit,'(a)')        '|    Duncan Poole; Scott Le Grand; Ross C. Walker'
  write(iout_unit,'(a)')        '|    "Routine microsecond molecular dynamics simulations'
  write(iout_unit,'(a)')        '|    with AMBER - Part I: Generalized Born", J. Chem.'
  write(iout_unit,'(a)')        '|    Theory Comput., 2012, 8 (5), pp1542-1555.'
#ifdef use_SPFP
  write(iout_unit,'(a)')        '|'
  write(iout_unit,'(a)')        '|  - Scott Le Grand; Andreas W. Goetz; Ross C. Walker'
  write(iout_unit,'(a)')        '|    "SPFP: Speed without compromise - a mixed precision'
  write(iout_unit,'(a)')        '|    model for GPU accelerated molecular dynamics'
  write(iout_unit,'(a)')        '|    simulations.", Comp. Phys. Comm., 2013, 184'
  write(iout_unit,'(a)')        '|    pp374-380, DOI: 10.1016/j.cpc.2012.09.022'
#endif
if (iamd/=0) then
  write(iout_unit,'(a)')        '|'
  write(iout_unit,'(a)')        '|    When publishing work that utilized the CUDA version'
  write(iout_unit,'(a)')        '|    of AMD, please cite the following in addition to'
  write(iout_unit,'(a)')        '|    the regular AMBER citations:' 
  write(iout_unit,'(a)')        '|' 
  write(iout_unit,'(a)')        '|  - Levi C. T. Pierce; Romelia Salomon-Ferrer; '
  write(iout_unit,'(a)')        '|    Cesar Augusto F de Oliveira; J. Andrew McCammon'
  write(iout_unit,'(a)')        '|    and Ross C. Walker "Routine access to milli-second '
  write(iout_unit,'(a)')        '|    time scales with accelerated molecular dynamics".' 
  write(iout_unit,'(a)')        '|    J. Chem. Theory Comput., 2012, 8(9), pp2997-3002.'
  write(iout_unit,'(a)')        '|    DOI: 10.1021/ct300284c.'
  write(iout_unit,'(a)')        '|'
endif
  write(iout_unit,'(a)')        '|'
  write(iout_unit,'(a)')        '|--------------------------------------------------------'
  write(iout_unit,'(a)')        ' '

  return

end subroutine write_cuda_citation

subroutine gpu_write_memory_info(iout_unit)

  implicit none

! Formal arguments:

  integer, intent(in)   :: iout_unit !Unit to write info to.

! Local variables:
  integer :: gpumemory, cpumemory
  
  call gpu_get_memory_info(gpumemory, cpumemory)

  write(iout_unit, '(a)')    '| GPU memory information (estimate):'
  write(iout_unit,'(a,i9)')  '| KB of GPU memory in use: ', gpumemory
  write(iout_unit,'(a,i9/)') '| KB of CPU memory in use: ', cpumemory

  return

end subroutine gpu_write_memory_info

