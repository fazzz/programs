#include "copyright.i"

!*******************************************************************************
!
! Module: amd_mod
!
! Description: 
!
! Module for controlling accelerated molecular dynamics calculations
!
! Written by Romelia Salomon, 5/2011
!              
!*******************************************************************************

module amd_mod

  use file_io_dat_mod
  use parallel_dat_mod
  use mdin_ctrl_dat_mod
  use dihedrals_mod

  implicit none

! Everything is private by default

! Constants
 
  double precision, parameter :: J_PER_CAL = 4.184d0 !  This is defined as the thermochemical calorie
  double precision, parameter :: JPKC = J_PER_CAL * 1000.0d0 !kilocalories per joule
  double precision, parameter :: BOLTZMANN = 1.380658d-23 !Boltzmann's constant in J/K
  double precision, parameter :: AVOGADRO = 6.0221367d+23 !Avogadro's number
  double precision, parameter :: KB = (BOLTZMANN * AVOGADRO) / JPKC !Boltzmann's constant in internal units

! Variables
!

!AMD
  double precision, save              :: totalenergy,tboostall,fwgt,totdih,tboost,dihsum
  integer, save                       :: num_amd_recs, num_amd_lag, amd_ntwx
  double precision, allocatable, save             :: amd_weights_and_energy(:)

contains

!*******************************************************************************
!
! Subroutine: amd_setup
!
! Description: Sets up the AMD run. Opens the log file, sets up exchange tables
!              etc.
!
!*******************************************************************************

subroutine amd_setup(ntwx)
  
  use pmemd_lib_mod

  implicit none

! Formal arguments:
  integer               :: ntwx

! Local variables:
  integer :: alloc_failed


  allocate(amd_weights_and_energy(6), stat = alloc_failed)

  if (alloc_failed .ne. 0) call setup_alloc_error

  fwgt = 1.0
  fwgtd = 1.0
  tboostall = 0.0
  tboost = 0.0
  num_amd_recs = 1
  num_amd_lag = 0
  amd_ntwx = ntwx

#ifdef CUDA
  call gpu_amd_setup(iamd, iamdlag, ntwx, EthreshP, alphaP, EthreshD, alphaD, temp0)
#endif
  return
end subroutine amd_setup

!*******************************************************************************
!
! Subroutine:  calculate_amd_dih_weights
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine calculate_amd_dih_weights(atm_cnt,totdih_ene,frc,crd)

  implicit none

! Formal arguments:
  integer                       :: atm_cnt
  double precision              :: totdih_ene
  double precision              :: frc(3, atm_cnt)
  double precision              :: crd(3, atm_cnt)

! Calculate the boosting weight for amd

  totdih = 0.0d0
  tboost = 0.0d0
  fwgtd   = 1.0d0
  totdih = totdih_ene

  if(totdih.le.EthreshD)then
    if(num_amd_lag .eq. 0)then
    tboost = ((EthreshD - totdih)**2)/ &
             ((alphaD + (EthreshD - totdih))*temp0*KB)
    fwgtd = (alphaD**2)/((alphaD + EthreshD - totdih)**2)
    end if
  end if
  if (ntf .le. 5) then
    if (cit_nphih .gt. 0) then
      call get_dihed_forces_amd(cit_nphih, cit_dihed, crd, frc)
    end if
  end if

  if (ntf .le. 6) then
    if (cit_nphia .gt. 0) then
      call get_dihed_forces_amd(cit_nphia, cit_dihed(cit_nphih + 1), crd, frc)
    end if
  end if

  return

end subroutine calculate_amd_dih_weights

!*******************************************************************************
!
! Subroutine:  calculate_amd_total_weights
!
! Description: <TBS>
!              
!*******************************************************************************

#ifdef MPI
subroutine calculate_amd_total_weights(atm_cnt,tot_potenergy,totdih_ene,frc,crd,my_atm_lst)
#else
subroutine calculate_amd_total_weights(atm_cnt,tot_potenergy,totdih_ene,frc,crd)
#endif
  implicit none

! Formal arguments:
  integer                       :: atm_cnt
  double precision              :: tot_potenergy
  double precision              :: totdih_ene
  double precision              :: frc(3, atm_cnt)
  double precision              :: crd(3, atm_cnt)
#ifdef MPI
  integer                       :: my_atm_lst(*)
#endif

! Local variables:
  integer              :: atm_lst_idx,i

!AMD DUAL BOOST CALC START
   if(iamd.gt.0)then
     totalenergy = 0.0d0
     tboostall = 0.0d0
     fwgt   = 1.0d0
     totalenergy = tot_potenergy + (tboost*temp0*KB)
     if (((iamd == 1).or.(iamd == 3)) .and. (totalenergy.le.EthreshP)) then
       if (num_amd_lag .eq. 0) then
       tboostall = ((EthreshP - totalenergy)**2)/ &
                 ((alphaP + (EthreshP - totalenergy))*temp0*KB)
       fwgt = (alphaP**2)/((alphaP + EthreshP - totalenergy)**2)
#ifdef MPI
       do atm_lst_idx = 1, my_atm_cnt
         i = my_atm_lst(atm_lst_idx)
#else
       do i = 1, atm_cnt
#endif
         frc(:,i) = frc(:,i)*fwgt
       enddo 
       end if
     end if
     if (master.and. (num_amd_recs.eq.amd_ntwx)) then
       amd_weights_and_energy(1) = tot_potenergy 
       amd_weights_and_energy(2) = totdih_ene
       amd_weights_and_energy(3) = fwgt
       amd_weights_and_energy(4) = fwgtd
       amd_weights_and_energy(5) = tboostall
       amd_weights_and_energy(6) = tboost
     endif
     num_amd_recs = num_amd_recs +1
     if(num_amd_lag .eq. iamdlag)then
       num_amd_lag = 0
     else
       num_amd_lag = num_amd_lag +1
     endif 
   end if

  return

end subroutine calculate_amd_total_weights


!*******************************************************************************
!
! Subroutine:  calculate_amd_total_weights_gb
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine calculate_amd_total_weights_gb(atm_cnt,tot_potenergy,totdih_ene,frc,crd)

  use gb_parallel_mod

  implicit none

! Formal arguments:
  integer                       :: atm_cnt
  double precision              :: tot_potenergy
  double precision              :: totdih_ene
  double precision              :: frc(3, atm_cnt)
  double precision              :: crd(3, atm_cnt)

! Local variables:
  integer              :: atm_lst_idx,i

!AMD DUAL BOOST CALC START
   if(iamd.gt.0)then
     totalenergy = 0.0d0
     tboostall = 0.0d0
     fwgt   = 1.0d0
     totalenergy = tot_potenergy + (tboost*temp0*KB)
     if (((iamd == 1).or.(iamd == 3)) .and. (totalenergy.le.EthreshP)) then
       if (num_amd_lag .eq. 0) then
       tboostall = ((EthreshP - totalenergy)**2)/ &
                 ((alphaP + (EthreshP - totalenergy))*temp0*KB)
       fwgt = (alphaP**2)/((alphaP + EthreshP - totalenergy)**2)
#ifdef MPI
       call gb_amd_apply_weight(atm_cnt,frc,fwgt)
#else
       frc(:,:) = frc(:,:)*fwgt
#endif
       end if
     end if
     if (master.and. (num_amd_recs.eq.amd_ntwx)) then
       amd_weights_and_energy(1) = tot_potenergy 
       amd_weights_and_energy(2) = totdih_ene
       amd_weights_and_energy(3) = fwgt
       amd_weights_and_energy(4) = fwgtd
       amd_weights_and_energy(5) = tboostall
       amd_weights_and_energy(6) = tboost
     endif

     num_amd_recs = num_amd_recs +1
     if(num_amd_lag .eq. iamdlag)then
       num_amd_lag = 0
     else
       num_amd_lag = num_amd_lag +1
     endif 
   end if

  return

end subroutine calculate_amd_total_weights_gb


!*******************************************************************************
!
! Subroutine:  write_amd_weights
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine write_amd_weights(ntwx,total_nstep)

  use file_io_mod
!  use mdin_ctrl_dat_mod

  implicit none

! Formal arguments:

  integer               :: ntwx
  integer               :: total_nstep

! Local variables:

    write(amdlog,'(2x,2i10,6f22.12)') ntwx,total_nstep, &
     amd_weights_and_energy(1), amd_weights_and_energy(2), &
     amd_weights_and_energy(3), amd_weights_and_energy(4), &
     amd_weights_and_energy(5), amd_weights_and_energy(6)
  
  num_amd_recs = 1

  return

end subroutine write_amd_weights

end module amd_mod
