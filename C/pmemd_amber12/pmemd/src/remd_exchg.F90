#include "copyright.i"

!*******************************************************************************
!
! Module: remd_exchg_mod
!
! Description: 
!
! This module holds all of the subroutines that conduct exchanges.
! 
! This module needed to be introduced external to the REMD module because many
! of the more advanced exchange subroutines introduced dependencies for remd_mod
! that resulted in cyclic dependencies that prevented compiling. It's probably
! good to have them separate, anyway
!
!*******************************************************************************

module remd_exchg_mod

#ifdef MPI
  use parallel_dat_mod
  use random_mod, only : random_state, amrset_gen, amrand_gen
  use remd_mod

  implicit none

! Constants

  double precision, parameter :: ONEKB = 503.01d0 ! 1/Kb in internal units
  double precision, parameter :: TINY = 1.d-6

! Random number generator

  type(random_state), save, private :: remd_rand_gen

contains

!*******************************************************************************
!
! Subroutine: temperature_exchange
!
! Description: Performs the temperature replica exchange attempt. We don't need
!              any extra energies here.
!
!*******************************************************************************

subroutine temperature_exchange(atm_cnt, vel, my_pot_ene_tot, t_dim, &
                                actual_temperature, print_exch_data, mdloop)

  use mdin_ctrl_dat_mod, only : temp0
  use pmemd_lib_mod,     only : mexit

  implicit none
 
! Passed variables

  integer, intent(in)             :: atm_cnt, mdloop

  double precision, intent(in)    :: my_pot_ene_tot

  double precision, intent(inout) :: vel(3, atm_cnt)
  
  double precision, intent(in)    :: actual_temperature
  
  integer, intent(in)             :: t_dim ! which dimension T-exchanges are

  logical, intent(in)             :: print_exch_data

! Local variables

  type :: exchange_data ! facilitates gather for remlog
    sequence
    double precision :: scaling
    double precision :: real_temp
    double precision :: pot_ene_tot
    double precision :: temp0
    double precision :: new_temp0
    double precision :: struct_num
  end type exchange_data

  integer, parameter :: SIZE_EXCHANGE_DATA = 6 ! for mpi_gather

  type (exchange_data) :: exch_buffer(numgroups) ! gather buffer

  type (exchange_data) :: my_exch_data ! stores data for this replica

  double precision  :: neighbor_pot_ene
  double precision  :: neighbor_temperature
  double precision  :: my_temperature
  double precision  :: delta
  double precision  :: metrop
  double precision  :: random_value
  double precision  :: success_ratio

  integer           :: my_idx
  integer           :: neighbor_idx
  integer           :: neighbor_rank
  integer           :: replica_indices(numgroups)
  integer           :: i
  integer           :: success_array(numgroups)
  integer           :: success_buf(numgroups)

  logical           :: success
  logical           :: i_do_exchg

  integer, dimension(mpi_status_size) :: stat_array ! needed for mpi_recv

! Explanation of local variables:
!
! replica_indices: this is an array in which each replica's index occupies that
!                  replica's position in the array. It's used during mpi_gather
!                  by the master so it knows how to write out the remlog
!
! delta: The calculated DELTA value used as part of the detailed balance calc.
!
! metrop: The actual value compared to a random_value to decide success
!
! success: Did this exchange attempt succeed?
!
! success_array: buffer for keeping track of which temperatures exchanged
!                successfully. Used in a mpi_reduce at the end
!
! success_buf: buffer for mpi_reduce call on success_array
!
! my_idx: My index in the temperature table (num_from_min + 1)
!
! neighbor_idx: Neighbor index in the temperature table
! 
! neighbor_rank: the master_rank of the replica with whom we need to be
!                doing ALL of our communication. Trying to do everything
!                with the 'partners' array is a bit too confusing.


! Set the variables that we know at the beginning

  my_exch_data%real_temp   = actual_temperature
  my_exch_data%temp0       = temp0
  my_exch_data%new_temp0   = temp0 ! changed later if exch. succeeds
  my_exch_data%struct_num  = -1    ! not implemented yet
  my_exch_data%pot_ene_tot = my_pot_ene_tot
  my_exch_data%scaling     = -1.d0
  my_idx                   = num_from_min + 1
  success_array(:)         = 0

  if (master) then

    ! Call amrand_gen here to generate a random number. Sander does this
    ! to pick out structures from a reservoir which pmemd doesn't implement
    ! (yet). To keep the random #s synched, though, and to ensure proper
    ! regression testing, this must be done to keep the random #s synched
    ! b/w sander and pmemd

    call amrand_gen(remd_rand_gen, random_value)

#ifdef VERBOSE_REMD
    write(mdout, '(26("="),a,26("="))') 'REMD EXCHANGE CALCULATION'
    write(mdout, '(a,i10,a,i1)') 'Exch= ', exch_num(t_dim), ' RREMD= ', 0
#endif

    ! We alternate exchanges: temp up followed by temp down
    if (even_exchanges(t_dim)) then
      i_do_exchg = even_replica(t_dim)
    else
      i_do_exchg = .not. even_replica(t_dim)
    end if

    even_exchanges(t_dim) = .not. even_exchanges(t_dim)

    ! Find out what rank our neighbor is (-1 for indexing from 0), and
    ! what its index is in the temperature table. If we are controlling
    ! the exchange, then our partner will be our 2nd partner (above). 
    ! Otherwise it will be our 1st partner (below)
    
    if (i_do_exchg) then

      neighbor_rank = partners(t_dim, 2) - 1
      neighbor_temperature = temperatures(partners(t_dim, 2))

      if (my_idx .eq. numgroups) then
        neighbor_idx = 1
      else
        neighbor_idx = my_idx + 1
      end if

    else

      neighbor_rank = partners(t_dim, 1) - 1
      neighbor_temperature = temperatures(partners(t_dim, 1))

      if (my_idx .eq. 1) then
        neighbor_idx = numgroups
      else
        neighbor_idx = my_idx - 1
      end if

    end if

    ! Trade potential energies

    call mpi_sendrecv(my_pot_ene_tot, 1, mpi_double_precision, &
                      neighbor_rank, remd_tag, &
                      neighbor_pot_ene, 1, mpi_double_precision, &
                      neighbor_rank, remd_tag, pmemd_master_comm, &
                      stat_array, err_code_mpi)

#ifdef VERBOSE_REMD
    write(mdout, '(a,f6.2,2(a7,i2),a7,f10.2)') 'Replica          Temp= ', &
      my_exch_data%temp0, ' Indx= ', my_idx, ' Rep#= ', master_rank+1, &
      ' EPot= ', my_pot_ene_tot
#endif

    ! Calculate the exchange probability. Note that if delta is negative, then
    ! metrop will be > 1, and it will always succeed (so we don't need an extra
    ! conditional)
    
    if (i_do_exchg) then
#ifdef VERBOSE_REMD
      write(mdout, '(a,f6.2,2(a7,i2),a7,f10.2)') 'Partner          Temp= ', &
        neighbor_temperature, ' Indx= ', neighbor_idx, ' Rep#= ', &
        neighbor_rank + 1, ' EPot= ', neighbor_pot_ene
#endif

      delta = (my_pot_ene_tot - neighbor_pot_ene) * &
              (my_exch_data%temp0 - neighbor_temperature) * ONEKB / &
              (my_exch_data%temp0 * neighbor_temperature)

      metrop = exp(-delta)

      call amrand_gen(remd_rand_gen, random_value)

      success = random_value .lt. metrop

      ! Let my neighbor know about our success

      call mpi_send(success, 1, mpi_logical, neighbor_rank, &
                    remd_tag, pmemd_master_comm, err_code_mpi)

      if (success) then
        success_array(my_idx) = 1
        my_exch_data%new_temp0 = neighbor_temperature
        my_exch_data%scaling = sqrt(my_exch_data%new_temp0 / my_exch_data%temp0)
      end if

#ifdef VERBOSE_REMD
      write(mdout,'(a8,E16.6,a8,E16.6,a12,f10.2)') &
               "Metrop= ",metrop," delta= ",delta," o_scaling= ", &
               1 / my_exch_data%scaling
#endif

    else

#ifdef VERBOSE_REMD
      ! Even if we don't control exchange, still print REMD info to mdout

      write(mdout, '(a,f6.2,2(a7,i2),a7,f10.2)') 'Partner          Temp= ', &
        neighbor_temperature, ' Indx= ', neighbor_idx, ' Rep#= ', &
        neighbor_rank + 1, ' EPot= ', neighbor_pot_ene
      write(mdout, '(a)') 'Not controlling exchange.'
#endif

      call amrand_gen(remd_rand_gen, random_value) ! to stay synchronized

      ! Get the message from the exchanging replica about our success
      call mpi_recv(success, 1, mpi_logical, neighbor_rank, &
                    remd_tag, pmemd_master_comm, stat_array, err_code_mpi)

      ! We scale velocities only if we succeed
      if (success) then
        my_exch_data%new_temp0 = neighbor_temperature
        my_exch_data%scaling = sqrt(my_exch_data%new_temp0 / my_exch_data%temp0)
      end if

    end if

#ifdef VERBOSE_REMD
    write(mdout,'(a8,E16.6,a12,f10.2,a10,L1)') &
      'Rand=   ', random_value, ' MyScaling= ', my_exch_data%scaling, &
      ' Success= ', success
    
    write(mdout,'(24("="),a,24("="))') "END REMD EXCHANGE CALCULATION"
#endif

  end if ! master

  ! We'll broadcast all of the my_exch_data here once. From that, we can tell
  ! if the exchange succeeded based on whether or not the temperatures have
  ! changed. Then make sure we update temp0 so the simulation reflects that.

  call mpi_bcast(my_exch_data, SIZE_EXCHANGE_DATA, mpi_double_precision, &
                 0, pmemd_comm, err_code_mpi)

  temp0 = my_exch_data%new_temp0

  if (.not. master) &
    success = abs(my_exch_data%temp0 - my_exch_data%new_temp0) .gt. TINY

#ifdef VERBOSE_REMD
  if (master) &
    write(mdout, '(2(a,f6.2))') &
      'REMD: checking to see if bath T has changed: ', &
      my_exch_data%temp0, '->', my_exch_data%new_temp0
#endif

  if (success) &
    call rescale_velocities(atm_cnt, vel, my_exch_data%scaling)

  ! We want to keep track of our exchange successes, so gather them
  ! from all of the masters here, and increment them in exchange_successes
  ! We do this whether or not we print so we can keep track of exchange
  ! success ratios. For the index and exchange data, only gather those to
  ! the master_master if we actually plan to print. This should be much
  ! more efficient for high exchange attempts.

  ! Then recollect the temperatures as some may have changed, and reset our
  ! partners.

  if (master) then

    if (print_exch_data) then
      call mpi_gather(my_exch_data, SIZE_EXCHANGE_DATA, mpi_double_precision, &
                      exch_buffer, SIZE_EXCHANGE_DATA, mpi_double_precision, &
                      0, pmemd_master_comm, err_code_mpi)
    
      call mpi_gather(my_idx, 1, mpi_integer, &
                      replica_indices, 1, mpi_integer, &
                      0, pmemd_master_comm, err_code_mpi)
    end if

    call mpi_reduce(success_array, success_buf, numgroups, mpi_integer, &
                    mpi_sum, 0, pmemd_master_comm, err_code_mpi)

    call collect_temps

    call set_temp_partners(t_dim)

    if (master_master) &
      exchange_successes(t_dim,:) = exchange_successes(t_dim,:) + success_buf(:)

  end if

  ! Write the data to the remlog

  if (print_exch_data .and. master_master) then
    write(remlog, '(a,i8)') '# exchange ', sum(exch_num, 1)
    do i = 1, numgroups
      success_ratio = float(exchange_successes(t_dim, replica_indices(i))) / &
                      float(exch_num(t_dim)) * 2
      write(remlog, '(i2,6f10.2,i8)') &
            i, &
            exch_buffer(i)%scaling, &
            exch_buffer(i)%real_temp, &
            exch_buffer(i)%pot_ene_tot, &
            exch_buffer(i)%temp0, &
            exch_buffer(i)%new_temp0, &
            success_ratio, &
            int(exch_buffer(i)%struct_num)
    end do
  end if

  ! Increment exchange counter

  if (master) exch_num(t_dim) = exch_num(t_dim) + 1

  ! Let's let everyone catch up, I don't know if this is necessary...

  remd_modwt = .true.

  return

end subroutine temperature_exchange

!*******************************************************************************
!
! Subroutine: hamiltonian_exchange
!
! Description: Performs the Hamiltonian replica exchange attempt. Additional
!              energies need to be computed here.
!
!*******************************************************************************

subroutine hamiltonian_exchange(atm_cnt, my_atm_lst, crd, vel, frc, my_pot_ene_tot, &
                                gbl_img_atm_map, gbl_atm_img_map, h_dim)

  use gb_force_mod
  use mdin_ctrl_dat_mod
  use pme_force_mod

  implicit none

! Passed variables

  integer          :: atm_cnt
  integer          :: my_atm_lst(*)
  integer          :: gbl_img_atm_map(*)
  integer          :: gbl_atm_img_map(*)
  double precision :: crd(3,atm_cnt)
  double precision :: frc(3,atm_cnt)
  double precision :: vel(3,atm_cnt)
  integer          :: h_dim
  double precision :: my_pot_ene_tot

! Local variables

  type :: ene_temp ! Useful type to limit communication
    sequence
    double precision :: neighbor_rep ! neighbor replica # (for remlog printing)
    double precision :: energy_1 ! energy with MY coordinates
    double precision :: energy_2 ! energy with THEIR coordinates
    double precision :: left_fe  ! free energy for exchanging to the left
    double precision :: right_fe ! free energy for exchanging to the right
    double precision :: temperature
  end type ene_temp
  integer, parameter :: SIZE_ENE_TEMP = 6

  double precision :: my_pot_ene_tot_2   ! MY pot ene with THEIR coordinates
  double precision :: delta 
  double precision :: metrop
  double precision :: random_value
  double precision :: success_ratio
  double precision :: virial(3)          ! For pme_force_call
  double precision :: ekcmt(3)           ! For pme_force_call
  double precision :: pme_err_est        ! For pme_force_call
  double precision :: vel_scale_coff     ! For scaling velocities at different temperatures
  logical          :: i_do_exchg
  logical          :: rank_success(numgroups) ! log of all replica's success
  logical          :: success
  integer          :: success_array(numgroups) ! to keep track of successes
  integer          :: success_buf(numgroups) ! mpi_reduce buffer for successes
  character        :: success_char ! to print T or F if we succeeded in remlog
  integer          :: neighbor_rank
  integer          :: istat(mpi_status_size)
  integer          :: ier, i
  integer          :: rep

  type(gb_pot_ene_rec)  :: my_new_gb_ene  ! gb energy record for THEIR coords
  type(pme_pot_ene_rec) :: my_new_pme_ene ! pme energy record for THEIR coords
  type(ene_temp)   :: my_ene_temp
  type(ene_temp)   :: neighbor_ene_temp
  type(ene_temp)   :: ene_temp_buffer(numgroups) ! to collect/print REMD data


  ! Get coordinates and forces from GPU
#ifdef CUDA
  call gpu_download_crd(crd)
  call gpu_download_frc(frc)
#endif

  ! Set up some stuff for pme_force call:

  virial(:) = 0.d0
  ekcmt(:) = 0.d0


  ! First determine if we're controlling the exchange or not

  if (master) then

    ! Initialize some data

    success_array(:) = 0
    my_ene_temp%energy_1 = my_pot_ene_tot
    my_ene_temp%temperature = temp0

#ifdef VERBOSE_REMD
    write(mdout, '(24("="),a,24("="))') ' H-REMD EXCHANGE CALCULATION '
    write(mdout, '(a,i10,a,i1)') 'Exch= ', exch_num, ' RREMD= ', 0
#endif
    
    ! Alternate exchanges

    if (even_exchanges(h_dim)) then
      i_do_exchg = even_replica(h_dim)
    else
      i_do_exchg = .not. even_replica(h_dim)
    end if

    even_exchanges(h_dim) = .not. even_exchanges(h_dim)
    
    ! Set partner rank 
    if (i_do_exchg) then
      my_ene_temp%neighbor_rep = partners(h_dim, 2)
      neighbor_rank = partners(h_dim, 2) - 1
    else
      my_ene_temp%neighbor_rep = partners(h_dim, 1)
      neighbor_rank = partners(h_dim, 1) - 1
    end if

    ! Exchange coordinates with your neighbor

    call mpi_sendrecv(crd, atm_cnt * 3, mpi_double_precision, &
                      neighbor_rank, remd_tag, &
                      crd_temp, atm_cnt * 3, mpi_double_precision, &
                      neighbor_rank, remd_tag, &
                      pmemd_master_comm, istat, err_code_mpi)

  end if ! master

  ! Now broadcast the temporary coordinates to the rest of pmemd_comm. I don't
  ! know if this is necessary or not (i.e. is it done in the load balancing?)

  call mpi_bcast(crd_temp, atm_cnt * 3, mpi_double_precision, &
                 0, pmemd_comm, err_code_mpi)

#ifdef CUDA
  ! Update the crd and frc if we are using the GPU
  ! JMS: Is this step necessary??
  call gpu_upload_crd(crd_temp)
  call gpu_upload_frc(frc_temp)
#endif

  ! Now it's time to get the energy of the other structure:

  if (using_gb_potential) then

    call gb_force(atm_cnt, crd_temp, frc_temp, my_new_gb_ene, irespa, .true.)

    my_ene_temp%energy_2 = my_new_gb_ene%total

  else if (using_pme_potential) then
    
    ! First .true. is new_list, next is need_pot_enes, .false. is need_virials

    call pme_force(atm_cnt, crd_temp, frc_temp, gbl_img_atm_map, gbl_atm_img_map, &
                   my_atm_lst, .true., .true., .false., my_new_pme_ene, &
                   virial, ekcmt, pme_err_est)
    
    my_ene_temp%energy_2 = my_new_pme_ene%total

  end if


#ifdef CUDA
  ! Update the results from the GPU
  ! JMS: Is this necessary??
  call gpu_download_crd(crd_temp)
  call gpu_download_frc(frc_temp)
#endif

  ! Now trade that potential energy with your partner

  if (master) then

    call mpi_sendrecv(my_ene_temp, SIZE_ENE_TEMP, mpi_double_precision, &
                      neighbor_rank, remd_tag, &
                      neighbor_ene_temp, SIZE_ENE_TEMP, mpi_double_precision, &
                      neighbor_rank, remd_tag, &
                      pmemd_master_comm, istat, err_code_mpi)

  ! Now it's time to calculate the exchange criteria. This detailed balance was
  ! derived under the assumption that *if* we have different T's, then they are
  ! NOT exchanged. Swapping temperatures requires a slightly different equation
  ! and then requires you to update temp0 at the end with an mpi_sendrecv call.

    if (i_do_exchg) then

      delta = -ONEKB / temp0 * (my_ene_temp%energy_2 - my_ene_temp%energy_1) - &
              ONEKB / neighbor_ene_temp%temperature * &
              (neighbor_ene_temp%energy_2 - neighbor_ene_temp%energy_1)

      metrop = exp(delta)

      call amrand_gen(remd_rand_gen, random_value)

      success = random_value .lt. metrop

      ! Let my neighbor know about our success
      call mpi_send(success, 1, mpi_logical, neighbor_rank, &
                    remd_tag, pmemd_master_comm, err_code_mpi)

      if (success) then
#ifdef CUDA
        call gpu_download_vel(vel)
#endif
        success_array(master_rank+1) = 1
        call mpi_sendrecv_replace(vel, atm_cnt * 3, mpi_double_precision, &
                   neighbor_rank, remd_tag, &
                   neighbor_rank, remd_tag, &
                   pmemd_master_comm, istat, err_code_mpi)
        vel_scale_coff = sqrt(temp0 / neighbor_ene_temp%temperature)
      end if

#ifdef VERBOSE_REMD
      write(mdout, '(a,f16.6)') 'My Eptot_1:       ', my_ene_temp%energy_1
      write(mdout, '(a,f16.6)') 'My Eptot_2:       ', my_ene_temp%energy_2
      write(mdout, '(a,f16.6)') 'Neighbor Eptot_1: ', neighbor_ene_temp%energy_1
      write(mdout, '(a,f16.6)') 'Neighbor Eptot_2: ', neighbor_ene_temp%energy_2
      write(mdout, '(3(a,f16.6))') 'Delta= ', delta, ' Metrop= ', metrop, &
                                   ' Random #= ', random_value
      if (success) then
        write(mdout, '(a)') 'Exchange Succeeded!'
      else
        write(mdout, '(a)') 'Exchange Failed!'
      end if
#endif

    else
      
      call amrand_gen(remd_rand_gen, random_value) ! keep in sync

      call mpi_recv(success, 1, mpi_logical, neighbor_rank, &
                    remd_tag, pmemd_master_comm, istat, err_code_mpi)
      if (success) then
#ifdef CUDA
        call gpu_download_vel(vel)
#endif
        call mpi_sendrecv_replace(vel, atm_cnt * 3, mpi_double_precision, &
                      neighbor_rank, remd_tag, &
                      neighbor_rank, remd_tag, &
                      pmemd_master_comm, istat, err_code_mpi)
        vel_scale_coff = sqrt(temp0 / neighbor_ene_temp%temperature)
      end if

    end if ! i_do_exchg

  end if ! master
  
  ! Now broadcast our success for all to hear! This is the ONLY part that is not
  ! for masters only

  call mpi_bcast(success, 1, mpi_logical, 0, pmemd_comm, err_code_mpi)

#ifdef CUDA
  ! If the exchange succeeded, then the coordinates and forces on the gpu
  ! are already correct and we don't need to do anything. Otherwise, we
  ! need to send back the originals.
  if (.not. success) then
    ! Upload the orginal coordinates again
    call gpu_upload_crd(crd)
    call gpu_upload_frc(frc)
  else
    call mpi_bcast(vel,3 * atm_cnt, mpi_double_precision, 0, pmemd_comm,err_code_mpi)
    call mpi_bcast(vel_scale_coff,1 , mpi_double_precision,0,pmemd_comm,err_code_mpi)
    call gpu_upload_vel(vel)
    call rescale_velocities(atm_cnt, vel, vel_scale_coff)
  end if
#else
  ! Update our coordinates if we succeeded. Note that the exchange probability
  ! equation used currently satisfies detailed balance only for NOT exchanging
  ! temperatures. We also update our force array, since we won't have to run
  ! our force call on the first pass through our next MD loop in runmd
  if (success) then
    crd(1, :) = crd_temp(1, :)
    crd(2, :) = crd_temp(2, :)
    crd(3, :) = crd_temp(3, :)
    frc(1, :) = frc_temp(1, :)
    frc(2, :) = frc_temp(2, :)
    frc(3, :) = frc_temp(3, :)
    call mpi_bcast(vel,3 * atm_cnt , mpi_double_precision, 0, pmemd_comm,err_code_mpi)
    call mpi_bcast(vel_scale_coff,1 , mpi_double_precision,0,pmemd_comm,err_code_mpi)
    call rescale_velocities(atm_cnt, vel, vel_scale_coff)
  end if
#endif /* CUDA */

  ! Calculate the free energies for jumping right and jumping left. This follows
  ! free energy perturbation (FEP). The equation for FEP is
  !
  !                          /        -(Eb - Ea)     \
  ! D G_(a->b) = -Kb * T ln | < exp (-------------) > |
  !                          \          Kb * T       /  A
  !
  ! Where the exponential is averaged. Thus, we have to keep a running total of
  ! the exponential Eb - Ea term (total_right_fe/total_left_fe), which can be
  ! combined into a total D G for this step via the above equation. These are
  ! stored in the ene_temp type to minimize communications. They can be 
  ! calculated regardless of whether we succeeded or not!
  ! For the layout of the rem.log file: left = up, right = down.
  ! Then calculate the 'averaged' free energy right after, keeping in mind that
  ! right/left alternate. Therefore, we introduce a counter for the number of
  ! left-moves and number of right-moves to make sure we divide by the right
  ! number of samples.

  if (master) then
    if (i_do_exchg) then

      num_right_exchg = num_right_exchg + 1

      total_right_fe = total_right_fe + &
                   exp((my_ene_temp%energy_1 - neighbor_ene_temp%energy_2) * &
                   ONEKB / my_ene_temp%temperature)

      my_ene_temp%right_fe = my_ene_temp%temperature / ONEKB * &
                             log(total_right_fe / num_right_exchg)
   
      my_ene_temp%left_fe = 0.d0

    else

      num_left_exchg = num_left_exchg + 1
      
      total_left_fe = total_left_fe + &
                   exp((my_ene_temp%energy_1 - neighbor_ene_temp%energy_2) * &
                   ONEKB / my_ene_temp%temperature)

      my_ene_temp%left_fe = my_ene_temp%temperature / ONEKB * &
                            log(total_left_fe / num_left_exchg)

      my_ene_temp%right_fe = 0.d0

    end if

    ! Collect data and write it out to the rem.log

    call mpi_reduce(success_array, success_buf, numgroups, mpi_integer, &
                    mpi_sum, 0, pmemd_master_comm, err_code_mpi)
    
    call mpi_gather(my_ene_temp, SIZE_ENE_TEMP, mpi_double_precision, &
                    ene_temp_buffer, SIZE_ENE_TEMP, mpi_double_precision, &
                    0, pmemd_master_comm, err_code_mpi)

    call mpi_gather(success, 1, mpi_logical, rank_success, 1, mpi_logical, &
                    0, pmemd_master_comm, err_code_mpi)
    if (master_master) then

      exchange_successes(h_dim,:) = exchange_successes(h_dim,:) + success_buf(:)

      write(remlog, '(a,i8)') '# exchange ', sum(exch_num, 1)

      do rep = 1, numgroups

        ! Did we succeed or not? We need to print a character
        if (rank_success(rep)) then
          success_char = 'T'
        else
          success_char = 'F'
        end if

        ! Calculate the success ratio:

        success_ratio = float(exchange_successes(h_dim, rep)) / &
                        float(exch_num(h_dim)) * 2

        ! Print info to the remlog

        write(remlog, '(2i6,5f10.2,4x,a,2x,f10.2)') rep, &
              int(ene_temp_buffer(rep)%neighbor_rep), &
              ene_temp_buffer(rep)%temperature, &
              ene_temp_buffer(rep)%energy_1, &
              ene_temp_buffer(rep)%energy_2, &
              ene_temp_buffer(rep)%left_fe, &
              ene_temp_buffer(rep)%right_fe, &
              success_char, &
              success_ratio

      end do

    end if ! master_master

#ifdef VERBOSE_REMD
    write(mdout, '(26("="),a,26("="))') 'END H-REMD CALCULATION'
#endif

    exch_num(h_dim) = exch_num(h_dim) + 1

  end if ! master

  remd_modwt = .true.

  return

end subroutine hamiltonian_exchange

!*******************************************************************************
!
! Subroutine: setup_remd_randgen
!
! Description: 
!
! Initializes the REMD random number stream. This *should* go in remd.fpp
! (subroutine remd_setup), but the Intel compilers fail when it's taken from a
! different module. I think this is caused because the size of random_state is
! not a multiple of its largest member (causes a warning), so ifort confuses the
! types. This seems to me to be a compiler bug, but such it is.
!
!*******************************************************************************

subroutine setup_remd_randgen

  use mdin_ctrl_dat_mod, only : ig

  implicit none

  integer rand_seed

  rand_seed = mod(repnum-1, numgroups) / 2 * 2

  call amrset_gen(remd_rand_gen, ig + rand_seed)

  return

end subroutine setup_remd_randgen
#endif /* MPI */

end module remd_exchg_mod
