#ifndef DIRFRC_NOVEC

#ifdef BUILD_PAIRS_CALC_EFV
#ifdef BUILD_PAIRS_CALC_VEC_2CUT
!*******************************************************************************
!
! Internal Subroutine:  pairs_calc_efv_2cut
!
! Description:  Direct force computation on one atom.  Different implementations
!               are covered by conditional compilation in this file, but all
!               versions here should only be invoked if DIRFRC_NOVEC is
!               undefined
!*******************************************************************************

#ifdef DIRFRC_EFS
subroutine pairs_calc_efv_2cut(img_frc, img_crd, img_qterm, ef_tbl, eed_cub, &
                               ico, ipairs_sublst, img_iac, cn1, cn2, x_tran)
#else
subroutine pairs_calc_efv_2cut(img_frc, img_crd, img_qterm, eed_cub, &
                               ico, ipairs_sublst, img_iac, cn1, cn2, x_tran)
#endif

#else

!*******************************************************************************
!
! Internal Subroutine:  pairs_calc_efv
!
! Description:  Direct force computation on one atom.  Different implementations
!               are covered by conditional compilation in this file, but all
!               versions here should only be invoked if DIRFRC_NOVEC is
!               undefined
!*******************************************************************************

#ifdef DIRFRC_EFS
subroutine pairs_calc_efv(img_frc, img_crd, img_qterm, ef_tbl, eed_cub, &
                          ico, ipairs_sublst, img_iac, cn1, cn2, x_tran)
#else
subroutine pairs_calc_efv(img_frc, img_crd, img_qterm, eed_cub, &
                          ico, ipairs_sublst, img_iac, cn1, cn2, x_tran)
#endif

#endif /* BUILD_PAIRS_CALC_VEC_2CUT */
#endif /* BUILD_PAIRS_CALC_EFV */

#ifdef BUILD_PAIRS_CALC_FV
#ifdef BUILD_PAIRS_CALC_VEC_2CUT
!*******************************************************************************
!
! Internal Subroutine:  pairs_calc_fv_2cut
!
! Description:  Direct force computation on one atom.  Different implementations
!               are covered by conditional compilation in this file, but all
!               versions here should only be invoked if DIRFRC_NOVEC is
!               undefined
!*******************************************************************************

#ifdef DIRFRC_EFS
subroutine pairs_calc_fv_2cut(img_frc, img_crd, img_qterm, f_tbl, eed_cub, &
                              ico, ipairs_sublst, img_iac, cn1, cn2, x_tran)
#else
subroutine pairs_calc_fv_2cut(img_frc, img_crd, img_qterm, eed_cub, &
                              ico, ipairs_sublst, img_iac, cn1, cn2, x_tran)
#endif

#else

!*******************************************************************************
!
! Internal Subroutine:  pairs_calc_fv
!
! Description:  Direct force computation on one atom.  Different implementations
!               are covered by conditional compilation in this file, but all
!               versions here should only be invoked if DIRFRC_NOVEC is
!               undefined
!*******************************************************************************

#ifdef DIRFRC_EFS
subroutine pairs_calc_fv(img_frc, img_crd, img_qterm, f_tbl, eed_cub, &
                         ico, ipairs_sublst, img_iac, cn1, cn2, x_tran)
#else
subroutine pairs_calc_fv(img_frc, img_crd, img_qterm, eed_cub, &
                         ico, ipairs_sublst, img_iac, cn1, cn2, x_tran)
#endif

#endif /* BUILD_PAIRS_CALC_VEC_2CUT */
#endif /* BUILD_PAIRS_CALC_FV */

#ifdef BUILD_PAIRS_CALC_F
#ifdef BUILD_PAIRS_CALC_VEC_2CUT
!*******************************************************************************
!
! Internal Subroutine:  pairs_calc_f_2cut
!
! Description:  Direct force computation on one atom.  Different implementations
!               are covered by conditional compilation in this file, but all
!               versions here should only be invoked if DIRFRC_NOVEC is
!               undefined
!*******************************************************************************

#ifdef DIRFRC_EFS
subroutine pairs_calc_f_2cut(img_frc, img_crd, img_qterm, f_tbl, eed_cub, &
                             ico, ipairs_sublst, img_iac, cn1, cn2, x_tran)
#else
subroutine pairs_calc_f_2cut(img_frc, img_crd, img_qterm, eed_cub, &
                             ico, ipairs_sublst, img_iac, cn1, cn2, x_tran)
#endif

#else

!*******************************************************************************
!
! Internal Subroutine:  pairs_calc_f
!
! Description:  Direct force computation on one atom.  Different implementations
!               are covered by conditional compilation in this file, but all
!               versions here should only be invoked if DIRFRC_NOVEC is
!               undefined
!*******************************************************************************

#ifdef DIRFRC_EFS
subroutine pairs_calc_f(img_frc, img_crd, img_qterm, f_tbl, eed_cub, &
                        ico, ipairs_sublst, img_iac, cn1, cn2, x_tran)
#else
subroutine pairs_calc_f(img_frc, img_crd, img_qterm, eed_cub, &
                        ico, ipairs_sublst, img_iac, cn1, cn2, x_tran)
#endif

#endif /* BUILD_PAIRS_CALC_VEC_2CUT */
#endif /* BUILD_PAIRS_CALC_F */

  implicit none

! Formal arguments:

  double precision              :: img_frc(3, *)
  double precision              :: img_crd(3, *)
  double precision              :: img_qterm(*)
#ifdef DIRFRC_EFS
#ifdef NEED_ENE
  double precision, intent(in)  :: ef_tbl(*)
#else
  double precision, intent(in)  :: f_tbl(*)
#endif
#endif
  double precision, intent(in)  :: eed_cub(*)
  integer, intent(in)           :: ico(*)
  integer, intent(in)           :: ipairs_sublst(*)
  integer, intent(in)           :: img_iac(*)
  double precision, intent(in)  :: cn1(*), cn2(*)
  double precision, intent(in)  :: x_tran(1:3, 0:17)

! Local variables:

  integer, parameter            :: mask27 = Z"07FFFFFF"
  double precision, parameter   :: half = 1.d0/2.d0
  double precision, parameter   :: third = 1.d0/3.d0
  integer, parameter            :: vec_size = 128

  double precision      :: dumx, dumy, dumz

  double precision      :: cgi
  double precision      :: cgi_cgj
  double precision      :: b0, b1
  double precision      :: df
  double precision      :: dfx, dfy, dfz
  double precision      :: f6, r6, f12
#ifdef HAS_10_12
  double precision      :: f10, r10
#endif
#ifdef DIRFRC_EFS
  double precision      :: du, du2, du3 ! 'u' is another name for delr2
  double precision      :: del_efs
  double precision      :: dens_efs
  double precision      :: lowest_efs_u
#endif /* DIRFRC_EFS */

  ! Variables used with erfc switch table; name are historical:

  double precision      :: switch
  double precision      :: d_switch_dx
  double precision      :: x, dx, e3dx, e4dx2
  double precision      :: delr, delrinv

  integer               :: iaci
  integer               :: ic
  integer               :: ind
  integer               :: img_j
  integer               :: itran
  integer               :: sublst_head, sublst_max
  integer               :: vec_idx, vec_max
#ifdef SLOW_INDIRECTVEC
  integer               :: enc_img, vec_cnt
#else
  integer               :: nxt_idx, nxt_cnt
#endif

  double precision      :: delr2, delr2inv
#ifdef SLOW_INDIRECTVEC
  double precision      :: delx, dely, delz
#else
  integer               :: nxt(vec_size)
#endif

  ! Arrays that seem to cut cache problems on some machines...

  integer               :: img_j_vec(vec_size)
  double precision      :: del_vec(3, vec_size)
  double precision      :: delr2_vec(vec_size)

#ifdef DIRFRC_EFS
  dens_efs = efs_tbl_dens
  del_efs = 1.d0 / dens_efs
  lowest_efs_u = lowest_efs_delr2
#endif /* DIRFRC_EFS */

! First loop over the ee evaluation-only pairs:

  dumx = 0.d0
  dumy = 0.d0
  dumz = 0.d0

  cgi = img_qterm(img_i)

  sublst_head = 0
  sublst_max = ee_eval_cnt

  do while (sublst_head .lt. sublst_max)

    vec_max = min(sublst_max - sublst_head, vec_size)

#ifdef SLOW_INDIRECTVEC

    vec_cnt = 0

#ifdef DIRFRC_COMTRANS
    if (common_tran .eq. 1) then

      do vec_idx = 1, vec_max
        img_j = ipairs_sublst(sublst_head + vec_idx)
        delx = img_crd(1, img_j) + x_tran(1, 13)
        dely = img_crd(2, img_j) + x_tran(2, 13)
        delz = img_crd(3, img_j) + x_tran(3, 13)
        delr2 = delx * delx + dely * dely + delz * delz
#ifdef BUILD_PAIRS_CALC_VEC_2CUT
        if (delr2 .lt. es_cut2) then
#else
        if (delr2 .lt. max_nb_cut2) then
#endif
          vec_cnt = vec_cnt + 1
          img_j_vec(vec_cnt) = img_j
          del_vec(1, vec_cnt) = delx
          del_vec(2, vec_cnt) = dely
          del_vec(3, vec_cnt) = delz
          delr2_vec(vec_cnt) = delr2
        end if
      end do

    else
#endif /* DIRFRC_COMTRANS */

      do vec_idx = 1, vec_max
        enc_img = ipairs_sublst(sublst_head + vec_idx)
        img_j = iand(enc_img, mask27)
        itran = ishft(enc_img, -27)
        delx = img_crd(1, img_j) + x_tran(1, itran)
        dely = img_crd(2, img_j) + x_tran(2, itran)
        delz = img_crd(3, img_j) + x_tran(3, itran)
        delr2 = delx * delx + dely * dely + delz * delz
#ifdef BUILD_PAIRS_CALC_VEC_2CUT
        if (delr2 .lt. es_cut2) then
#else
        if (delr2 .lt. max_nb_cut2) then
#endif
          vec_cnt = vec_cnt + 1
          img_j_vec(vec_cnt) = img_j
          del_vec(1, vec_cnt) = delx
          del_vec(2, vec_cnt) = dely
          del_vec(3, vec_cnt) = delz
          delr2_vec(vec_cnt) = delr2
        end if
      end do

#ifdef DIRFRC_COMTRANS
    end if
#endif /* DIRFRC_COMTRANS */

#else

#ifdef DIRFRC_COMTRANS
    if (common_tran .eq. 1) then

      do vec_idx = 1, vec_max
        img_j = ipairs_sublst(sublst_head + vec_idx)
        img_j_vec(vec_idx) = img_j
        del_vec(1, vec_idx) = img_crd(1, img_j) + x_tran(1, 13)
        del_vec(2, vec_idx) = img_crd(2, img_j) + x_tran(2, 13)
        del_vec(3, vec_idx) = img_crd(3, img_j) + x_tran(3, 13)
        delr2_vec(vec_idx) = del_vec(1, vec_idx) * del_vec(1, vec_idx) + &
                             del_vec(2, vec_idx) * del_vec(2, vec_idx) + &
                             del_vec(3, vec_idx) * del_vec(3, vec_idx)
      end do

    else
#endif /* DIRFRC_COMTRANS */

      do vec_idx = 1, vec_max
        img_j = iand(ipairs_sublst(sublst_head + vec_idx), mask27)
        itran = ishft(ipairs_sublst(sublst_head + vec_idx), -27)
        img_j_vec(vec_idx) = img_j
        del_vec(1, vec_idx) = img_crd(1, img_j) + x_tran(1, itran)
        del_vec(2, vec_idx) = img_crd(2, img_j) + x_tran(2, itran)
        del_vec(3, vec_idx) = img_crd(3, img_j) + x_tran(3, itran)
        delr2_vec(vec_idx) = del_vec(1, vec_idx) * del_vec(1, vec_idx) + &
                             del_vec(2, vec_idx) * del_vec(2, vec_idx) + &
                             del_vec(3, vec_idx) * del_vec(3, vec_idx)
      end do

#ifdef DIRFRC_COMTRANS
    end if
#endif /* DIRFRC_COMTRANS */

    nxt_cnt = 0
    do vec_idx = 1, vec_max
#ifdef BUILD_PAIRS_CALC_VEC_2CUT
      if (delr2_vec(vec_idx) .lt. es_cut2) then
#else
      if (delr2_vec(vec_idx) .lt. max_nb_cut2) then
#endif
        nxt_cnt = nxt_cnt + 1
        nxt(nxt_cnt) = vec_idx
      end if 
    end do
#endif /* SLOW_INDIRECTVEC */

#ifdef SLOW_INDIRECTVEC
    do vec_idx = 1, vec_cnt
#else
    do nxt_idx = 1, nxt_cnt

      vec_idx = nxt(nxt_idx)
#endif /* SLOW_INDIRECTVEC */

      delr2 = delr2_vec(vec_idx)

      img_j = img_j_vec(vec_idx)
      cgi_cgj = cgi * img_qterm(img_j)

#ifdef DIRFRC_EFS
      if (delr2 .ge. lowest_efs_u) then

        ! Do the Coulomb part of the direct sum using efs: 

        ind = int(dens_efs * delr2)
        du = delr2 - dble(ind) * del_efs
        du2 = du * du
        du3 = du * du2

#ifdef NEED_ENE
        ind = ishft(ind, 3)             ! 8 * ind

        b0 = cgi_cgj * (ef_tbl(1 + ind) + du * ef_tbl(2 + ind) + &
             du2 * ef_tbl(3 + ind) + du3 * ef_tbl(4 + ind))

        df = cgi_cgj * (ef_tbl(5 + ind) + du * ef_tbl(6 + ind) + &
             du2 * ef_tbl(7 + ind) + du3 * ef_tbl(8 + ind))
#else
        ind = ishft(ind, 2)             ! 4 * ind

        df = cgi_cgj * (f_tbl(1 + ind) + du * f_tbl(2 + ind) + &
             du2 * f_tbl(3 + ind) + du3 * f_tbl(4 + ind))
#endif /* NEED_ENE */

#ifdef NEED_ENE
        eed_stk = eed_stk + b0
#endif /* NEED_ENE */

#ifdef NEED_VIR
        eedvir_stk = eedvir_stk - df * delr2
#endif /* NEED_VIR */

      else
#endif /* DIRFRC_EFS */
      
        ! Do the Coulomb part of the direct sum using erfc spline table: 

        delr = sqrt(delr2)
        delrinv = 1.d0 / delr

        x = dxdr * delr
        ind = int(eedtbdns_stk * x)
        dx = x - dble(ind) * del
        ind = ishft(ind, 2)             ! 4 * ind

        e3dx  = dx * eed_cub(3 + ind)
        e4dx2 = dx * dx *  eed_cub(4 + ind)

        switch = eed_cub(1 + ind) + &
                 dx * (eed_cub(2 + ind) + half * (e3dx + third * e4dx2))

        d_switch_dx = eed_cub(2 + ind) + e3dx + half * e4dx2

        b0 = cgi_cgj * delrinv * switch

        b1 = (b0 - cgi_cgj * d_switch_dx * dxdr)

#ifdef NEED_VIR
        eedvir_stk = eedvir_stk - b1
#endif /* NEED_VIR */
#ifdef NEED_ENE
        eed_stk = eed_stk + b0
#endif /* NEED_ENE */
        df = b1 * delrinv * delrinv

#ifdef DIRFRC_EFS
      end if
#endif /* DIRFRC_EFS */

      dfx = del_vec(1, vec_idx) * df
      dfy = del_vec(2, vec_idx) * df
      dfz = del_vec(3, vec_idx) * df

#ifdef NEED_VIR
      vxx = vxx - del_vec(1, vec_idx) * dfx
      vxy = vxy - del_vec(1, vec_idx) * dfy
      vxz = vxz - del_vec(1, vec_idx) * dfz
      vyy = vyy - del_vec(2, vec_idx) * dfy
      vyz = vyz - del_vec(2, vec_idx) * dfz
      vzz = vzz - del_vec(3, vec_idx) * dfz
#endif /* NEED_VIR */

      dumx = dumx + dfx
      dumy = dumy + dfy
      dumz = dumz + dfz

      img_frc(1, img_j) = img_frc(1, img_j) + dfx
      img_frc(2, img_j) = img_frc(2, img_j) + dfy
      img_frc(3, img_j) = img_frc(3, img_j) + dfz

    end do

    sublst_head = sublst_head + vec_max

  end do

  sublst_max = sublst_head + full_eval_cnt

  iaci = ntypes_stk * (img_iac(img_i) - 1)
  
  do while (sublst_head .lt. sublst_max)

    vec_max = min(sublst_max - sublst_head, vec_size)

#ifdef SLOW_INDIRECTVEC

    vec_cnt = 0

#ifdef DIRFRC_COMTRANS
    if (common_tran .eq. 1) then

      do vec_idx = 1, vec_max
        img_j = ipairs_sublst(sublst_head + vec_idx)
        delx = img_crd(1, img_j) + x_tran(1, 13)
        dely = img_crd(2, img_j) + x_tran(2, 13)
        delz = img_crd(3, img_j) + x_tran(3, 13)
        delr2 = delx * delx + dely * dely + delz * delz
        if (delr2 .lt. max_nb_cut2) then
          vec_cnt = vec_cnt + 1
          img_j_vec(vec_cnt) = img_j
          del_vec(1, vec_cnt) = delx
          del_vec(2, vec_cnt) = dely
          del_vec(3, vec_cnt) = delz
          delr2_vec(vec_cnt) = delr2
        end if
      end do

    else
#endif /* DIRFRC_COMTRANS */

      do vec_idx = 1, vec_max
        enc_img = ipairs_sublst(sublst_head + vec_idx)
        img_j = iand(enc_img, mask27)
        itran = ishft(enc_img, -27)
        delx = img_crd(1, img_j) + x_tran(1, itran)
        dely = img_crd(2, img_j) + x_tran(2, itran)
        delz = img_crd(3, img_j) + x_tran(3, itran)
        delr2 = delx * delx + dely * dely + delz * delz
        if (delr2 .lt. max_nb_cut2) then
          vec_cnt = vec_cnt + 1
          img_j_vec(vec_cnt) = img_j
          del_vec(1, vec_cnt) = delx
          del_vec(2, vec_cnt) = dely
          del_vec(3, vec_cnt) = delz
          delr2_vec(vec_cnt) = delr2
        end if
      end do

#ifdef DIRFRC_COMTRANS
    end if
#endif /* DIRFRC_COMTRANS */

#else

#ifdef DIRFRC_COMTRANS
    if (common_tran .eq. 1) then

      do vec_idx = 1, vec_max
        img_j = ipairs_sublst(sublst_head + vec_idx)
        img_j_vec(vec_idx) = img_j
        del_vec(1, vec_idx) = img_crd(1, img_j) + x_tran(1, 13)
        del_vec(2, vec_idx) = img_crd(2, img_j) + x_tran(2, 13)
        del_vec(3, vec_idx) = img_crd(3, img_j) + x_tran(3, 13)
        delr2_vec(vec_idx) = del_vec(1, vec_idx) * del_vec(1, vec_idx) + &
                             del_vec(2, vec_idx) * del_vec(2, vec_idx) + &
                             del_vec(3, vec_idx) * del_vec(3, vec_idx)
      end do

    else
#endif /* DIRFRC_COMTRANS */

      do vec_idx = 1, vec_max
        img_j = iand(ipairs_sublst(sublst_head + vec_idx), mask27)
        itran = ishft(ipairs_sublst(sublst_head + vec_idx), -27)
        img_j_vec(vec_idx) = img_j
        del_vec(1, vec_idx) = img_crd(1, img_j) + x_tran(1, itran)
        del_vec(2, vec_idx) = img_crd(2, img_j) + x_tran(2, itran)
        del_vec(3, vec_idx) = img_crd(3, img_j) + x_tran(3, itran)
        delr2_vec(vec_idx) = del_vec(1, vec_idx) * del_vec(1, vec_idx) + &
                             del_vec(2, vec_idx) * del_vec(2, vec_idx) + &
                             del_vec(3, vec_idx) * del_vec(3, vec_idx)
      end do

#ifdef DIRFRC_COMTRANS
    end if
#endif /* DIRFRC_COMTRANS */

    nxt_cnt = 0
    do vec_idx = 1, vec_max
      if (delr2_vec(vec_idx) .lt. max_nb_cut2) then
        nxt_cnt = nxt_cnt + 1
        nxt(nxt_cnt) = vec_idx
      end if 
    end do
#endif /* SLOW_INDIRECTVEC */

#ifdef SLOW_INDIRECTVEC
    do vec_idx = 1, vec_cnt
#else
    do nxt_idx = 1, nxt_cnt

      vec_idx = nxt(nxt_idx)
#endif /* SLOW_INDIRECTVEC */

      delr2 = delr2_vec(vec_idx)

      img_j = img_j_vec(vec_idx)
      ic = ico(iaci + img_iac(img_j))

#ifdef BUILD_PAIRS_CALC_VEC_2CUT
      if (delr2 .lt. es_cut2) then
#endif
      cgi_cgj = cgi * img_qterm(img_j)

#ifdef DIRFRC_EFS
      if (delr2 .ge. lowest_efs_u) then

        ! Do the Coulomb part of the direct sum using efs: 

        ind = int(dens_efs * delr2)
        du = delr2 - dble(ind) * del_efs
        du2 = du * du
        du3 = du * du2

#ifdef NEED_ENE
        ind = ishft(ind, 3)             ! 8 * ind

        b0 = cgi_cgj * (ef_tbl(1 + ind) + du * ef_tbl(2 + ind) + &
             du2 * ef_tbl(3 + ind) + du3 * ef_tbl(4 + ind))

        df = cgi_cgj * (ef_tbl(5 + ind) + du * ef_tbl(6 + ind) + &
             du2 * ef_tbl(7 + ind) + du3 * ef_tbl(8 + ind))
#else
        ind = ishft(ind, 2)             ! 4 * ind

        df = cgi_cgj * (f_tbl(1 + ind) + du * f_tbl(2 + ind) + &
             du2 * f_tbl(3 + ind) + du3 * f_tbl(4 + ind))
#endif /* NEED_ENE */

#ifdef NEED_ENE
        eed_stk = eed_stk + b0
#endif /* NEED_ENE */

#ifdef NEED_VIR
        eedvir_stk = eedvir_stk - df * delr2
#endif /* NEED_VIR */

        delr2inv = 1.d0 / delr2

      else
#endif /* DIRFRC_EFS */
      
        ! Do the Coulomb part of the direct sum using erfc spline table: 

        delr = sqrt(delr2)
        delrinv = 1.d0 / delr

        x = dxdr * delr
        ind = int(eedtbdns_stk * x)
        dx = x - dble(ind) * del
        ind = ishft(ind, 2)             ! 4 * ind

        e3dx  = dx * eed_cub(3 + ind)
        e4dx2 = dx * dx *  eed_cub(4 + ind)

        switch = eed_cub(1 + ind) + &
                 dx * (eed_cub(2 + ind) + half * (e3dx + third * e4dx2))

        d_switch_dx = eed_cub(2 + ind) + e3dx + half * e4dx2

        b0 = cgi_cgj * delrinv * switch

        b1 = (b0 - cgi_cgj * d_switch_dx * dxdr)

#ifdef NEED_VIR
        eedvir_stk = eedvir_stk - b1
#endif /* NEED_VIR */
#ifdef NEED_ENE
        eed_stk = eed_stk + b0
#endif /* NEED_ENE */
        delr2inv = delrinv * delrinv
        df = b1 * delr2inv

#ifdef DIRFRC_EFS
      end if
#endif /* DIRFRC_EFS */

#ifdef BUILD_PAIRS_CALC_VEC_2CUT
      else

        delr2inv = 1.d0 / delr2
        df = 0.d0

      end if
#endif /* BUILD_PAIRS_CALC_VEC_2CUT */

#ifdef HAS_10_12
      if (ic .gt. 0) then
#endif
        r6 = delr2inv * delr2inv * delr2inv
        f6 = cn2(ic) * r6
        f12 = cn1(ic) * (r6 * r6)
#ifdef NEED_ENE
        evdw_stk = evdw_stk + f12 - f6
#endif /* NEED_ENE */
        df = df + (12.d0 * f12 - 6.d0 * f6) * delr2inv
#ifdef HAS_10_12
      else if (ic .lt. 0) then
        ! This code allows 10-12 terms; in many (most?) (all?) cases, the
        ! only "nominal" 10-12 terms are on waters, where the asol and bsol
        ! parameters are always zero; hence we can skip this part.
        ic = - ic
        r10 = delr2inv * delr2inv * delr2inv * delr2inv * delr2inv
        f10 = gbl_bsol(ic) * r10
        f12 = gbl_asol(ic) * (r10 * delr2inv)
#ifdef NEED_ENE
        ehb_stk = ehb_stk + f12 - f10
#endif /* NEED_ENE */
        df = df + (12.d0 * f12 - 10.d0 * f10) * delr2inv
      end if
#endif /* HAS_10_12 */

      dfx = del_vec(1, vec_idx) * df
      dfy = del_vec(2, vec_idx) * df
      dfz = del_vec(3, vec_idx) * df

#ifdef NEED_VIR
      vxx = vxx - del_vec(1, vec_idx) * dfx
      vxy = vxy - del_vec(1, vec_idx) * dfy
      vxz = vxz - del_vec(1, vec_idx) * dfz
      vyy = vyy - del_vec(2, vec_idx) * dfy
      vyz = vyz - del_vec(2, vec_idx) * dfz
      vzz = vzz - del_vec(3, vec_idx) * dfz
#endif /* NEED_VIR */

      dumx = dumx + dfx
      dumy = dumy + dfy
      dumz = dumz + dfz

      img_frc(1, img_j) = img_frc(1, img_j) + dfx
      img_frc(2, img_j) = img_frc(2, img_j) + dfy
      img_frc(3, img_j) = img_frc(3, img_j) + dfz

    end do

    sublst_head = sublst_head + vec_max

  end do

  img_frc(1, img_i) = img_frc(1, img_i) - dumx
  img_frc(2, img_i) = img_frc(2, img_i) - dumy
  img_frc(3, img_i) = img_frc(3, img_i) - dumz

  return

#ifdef BUILD_PAIRS_CALC_EFV
#undef NEED_ENE
#undef NEED_VIR
#ifdef BUILD_PAIRS_CALC_VEC_2CUT
end subroutine pairs_calc_efv_2cut
#else
end subroutine pairs_calc_efv
#endif
#endif /* BUILD_PAIRS_CALC_EFV */

#ifdef BUILD_PAIRS_CALC_FV
#undef NEED_VIR
#ifdef BUILD_PAIRS_CALC_VEC_2CUT
end subroutine pairs_calc_fv_2cut
#else
end subroutine pairs_calc_fv
#endif
#endif /* BUILD_PAIRS_CALC_FV */

#ifdef BUILD_PAIRS_CALC_F
#ifdef BUILD_PAIRS_CALC_VEC_2CUT
end subroutine pairs_calc_f_2cut
#else
end subroutine pairs_calc_f
#endif
#endif /* BUILD_PAIRS_CALC_F */

#endif /* !DIRFRC_NOVEC */
