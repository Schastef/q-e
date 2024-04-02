!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
SUBROUTINE dynmat_us()
  !-----------------------------------------------------------------------
  !! This routine calculates the electronic term: \(\langle\psi|V"-eS"|\psi\rangle\)
  !! of the dynamical matrix. Eq. (B32) of PRB 64, 235118 (2001) is calculated
  !! here. Eqs. (B33) and (B34) in \(\texttt{addusdynmat}\).
  !
  USE kinds,                ONLY : DP
  USE constants,            ONLY : tpi
  USE ions_base,            ONLY : nat, ityp, ntyp => nsp, tau
  USE uspp,                 ONLY : nkb, vkb
  USE scf,                  ONLY : rho
  USE fft_base,             ONLY : dfftp
  USE fft_interfaces,       ONLY : fwfft
  USE buffers,              ONLY : get_buffer
  USE gvect,                ONLY : g, ngm, igtongl
  USE wvfct,                ONLY : npwx, nbnd, wg, et
  USE lsda_mod,             ONLY : lsda, current_spin, isk, nspin
  USE vlocal,               ONLY : vloc
  USE klist,                ONLY : xk, ngk, igk_k
  USE wavefunctions,        ONLY : evc
  USE cell_base,            ONLY : omega, tpiba2
  USE uspp_param,           ONLY : nh, nhm
  USE noncollin_module,     ONLY : noncolin, npol
#if defined(__CUDA)
  USE becmod,               ONLY : calbec, bec_type, allocate_bec_type, &
                                   deallocate_bec_type, beccopy,        &
                                   allocate_bec_type_acc, becupdate,    &
                                   deallocate_bec_type_acc
#else
  USE becmod,               ONLY : calbec, bec_type, allocate_bec_type, &
                                   deallocate_bec_type, beccopy
#endif
  USE modes,                ONLY : u
  USE dynmat,               ONLY : dyn
  USE phus,                 ONLY : alphap
  USE units_lr,             ONLY : iuwfc, lrwfc
  USE io_global,            ONLY : stdout
  USE mp_pools,             ONLY : my_pool_id, inter_pool_comm
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
  USE Coul_cut_2D,          ONLY : do_cutoff_2D 
  USE Coul_cut_2D_ph,       ONLY : cutoff_dynmat0 

  USE lrus,                 ONLY : becp1
  USE qpoint,               ONLY : nksq, ikks
  USE control_lr,           ONLY : nbnd_occ, lgamma
  USE uspp_init,            ONLY : init_us_2
  USE control_flags,        ONLY : offload_type

  IMPLICIT NONE
  
  INTEGER :: icart, jcart, na_icart, na_jcart, na, ng, nt, ik, &
       ig, is, ibnd, nu_i, nu_j, ijkb0, ikb, jkb, ih, jh, ikk, &
       js,  ijs, npw
  ! counters
  ! ikk: record position of wfc at k

  REAL(DP) :: gtau, fac, wgg
  ! the product G*\tau_s
  ! auxiliary variable
  ! the true weight of a K point

  COMPLEX(DP) :: work, dynwrk (3 * nat, 3 * nat), fact, tmpdynwrk
  ! work space
  TYPE (bec_type) :: gammap(3,3)
  COMPLEX(DP), ALLOCATABLE :: rhog (:), aux1 (:,:), work1 (:), &
               work2 (:), deff_nc(:,:,:,:)
  REAL(DP), ALLOCATABLE :: deff(:,:,:)
  ! fourier transform of rho
  ! the second derivative of the beta
  ! work space
#if defined(__CUDA)
  TYPE (bec_type) :: bectmp
#endif

  CALL start_clock ('dynmat_us')
  ALLOCATE (rhog  ( dfftp%nnr))
  ALLOCATE (work1 ( npwx))
  ALLOCATE (work2 ( npwx))
  ALLOCATE (aux1  ( npwx*npol , nbnd))
  IF (noncolin) THEN
     ALLOCATE (deff_nc( nhm, nhm, nat, nspin ))
  ELSE
     ALLOCATE (deff(nhm, nhm, nat ))
  END IF
#if defined(__CUDA)
  CALL allocate_bec_type_acc( nkb, nbnd, bectmp )
#endif
  DO icart=1,3
     DO jcart=1,3
        CALL allocate_bec_type(nkb,nbnd, gammap(icart,jcart))
     ENDDO
  ENDDO

  dynwrk (:,:) = (0.d0, 0.0d0)
  !
  !   We first compute the part of the dynamical matrix due to the local
  !   potential

  !   ... only the first pool does the calculation (no sum over k needed)
  IF ( my_pool_id /= 0 ) GOTO 100

  !
  rhog (:) = (0.d0, 0.d0)
  rhog (:) = CMPLX(rho%of_r(:, 1), 0.d0,kind=DP)

  CALL fwfft ('Rho', rhog, dfftp)
  CALL start_clock('dynus1')
  !
  ! there is a delta ss'
  !
  !$acc data copyin(dfftp, igtongl, ityp, rhog, tau, vloc) copy(dynwrk) present(g)
  !$acc data copyin(dfftp%nl) 
  !$acc parallel loop collapse(3) reduction(+:tmpdynwrk)
  DO na = 1, nat
     DO icart = 1, 3
        DO jcart = 1, 3
           na_icart = 3 * (na - 1) + icart
           na_jcart = 3 * (na - 1) + jcart
           tmpdynwrk = (0.0,0.0)
           !$acc loop vector reduction(+:tmpdynwrk)
           DO ng = 1, ngm
              gtau = tpi * (g (1, ng) * tau (1, na) + &
                            g (2, ng) * tau (2, na) + &
                            g (3, ng) * tau (3, na) )
              fac = omega * vloc (igtongl (ng), ityp (na) ) * tpiba2 * &
                   ( DBLE (rhog (dfftp%nl (ng) ) ) * COS (gtau) - &
                    AIMAG (rhog (dfftp%nl (ng) ) ) * SIN (gtau) )
              tmpdynwrk = tmpdynwrk - fac * g (icart, ng) * g (jcart, ng)
           END DO
           dynwrk (na_icart, na_jcart) = tmpdynwrk
        ENDDO
     ENDDO
  ENDDO
  !$acc end data
  !$acc end data
  IF (do_cutoff_2D) call cutoff_dynmat0(dynwrk, rhog)  

  CALL mp_sum (dynwrk, intra_bgrp_comm)
  CALL stop_clock('dynus1')
  !
  ! each pool contributes to next term
  !
100 CONTINUE
  !
  ! Here we compute  the nonlocal Ultra-soft contribution
  !
  CALL start_clock('dynusl')
  !$acc data present(g,igk_k) copyin(xk,evc) create(aux1)
  DO ik = 1, nksq
     ikk = ikks(ik)
     IF (lsda) current_spin = isk (ikk)
     npw = ngk(ikk)
     IF (nksq > 1) THEN
             CALL get_buffer (evc, lrwfc, iuwfc, ikk)
             !$acc update device(evc)
     ENDIF
     CALL start_clock('dynus2')
     CALL init_us_2 (npw, igk_k(1,ikk), xk (1, ikk), vkb, .true.)
     !
     !    We first prepare the gamma terms, which are the second derivatives
     !    becp terms.
     !
     DO icart = 1, 3
        DO jcart = 1, icart
           !$acc kernels
           aux1(:,:)=(0.d0,0.d0)
           !$acc end kernels
           !$acc parallel 
           !$acc loop collapse(2)
           DO ibnd = 1, nbnd
              DO ig = 1, npw
                 aux1 (ig, ibnd) = - evc (ig, ibnd) * tpiba2 * &
                      (xk (icart, ikk) + g (icart, igk_k(ig,ikk) ) ) * &
                      (xk (jcart, ikk) + g (jcart, igk_k(ig,ikk) ) )
              ENDDO
           END DO
           IF (noncolin) THEN
              !$acc loop collapse(2)
              DO ibnd = 1, nbnd
                 DO ig = 1, npw
                    aux1 (ig+npwx, ibnd) = - evc (ig+npwx, ibnd) * tpiba2 * &
                      (xk (icart, ikk) + g (icart, igk_k(ig,ikk) ) ) * &
                      (xk (jcart, ikk) + g (jcart, igk_k(ig,ikk) ) )
                 ENDDO
              ENDDO
           END IF
           !$acc end parallel
#if defined(__CUDA)
           CALL calbec ( offload_type, npw, vkb, aux1, bectmp )
           CALL becupdate( offload_type, gammap, icart, 3, jcart, 3, bectmp )
#else
           CALL calbec ( offload_type, npw, vkb, aux1, gammap(icart,jcart) )
#endif
           IF (jcart < icart) &
#if defined(__CUDA)
             CALL becupdate( offload_type, gammap, jcart, 3, icart, 3, bectmp )
#else
             CALL beccopy (gammap(icart,jcart),gammap(jcart,icart), nkb, nbnd)
#endif
        ENDDO
     ENDDO
     CALL stop_clock('dynus2')
     !
     !   And then compute the contribution from the US pseudopotential
     !   which is  similar to the KB one
     !
     DO ibnd = 1, nbnd_occ (ikk)
        wgg = wg (ibnd, ikk)
        IF (noncolin) THEN
           CALL compute_deff_nc(deff_nc,et(ibnd,ikk))
        ELSE
           CALL compute_deff(deff,et(ibnd,ikk))
        ENDIF
        ijkb0 = 0
        DO nt = 1, ntyp
           DO na = 1, nat
              IF (ityp (na) == nt) THEN
                 DO icart = 1, 3
                    na_icart = 3 * (na - 1) + icart
                    DO jcart = 1, 3
                       na_jcart = 3 * (na - 1) + jcart
                       DO ih = 1, nh (nt)
                          ikb = ijkb0 + ih
                          DO jh = 1, nh (nt)
                             jkb = ijkb0 + jh
                             IF (noncolin) THEN
                                ijs=0
                                DO is=1,npol
                                   DO js=1,npol
                                      ijs=ijs+1
                                      dynwrk(na_icart,na_jcart) = &
                                        dynwrk(na_icart,na_jcart) + &
                                             wgg* deff_nc(ih,jh,na,ijs) * &
                                  (CONJG(gammap(icart,jcart)%nc(ikb,is,ibnd))*&
                                     becp1(ik)%nc (jkb, js, ibnd) + &
                                     CONJG(becp1(ik)%nc(ikb, is, ibnd) ) * &
                                     gammap(icart,jcart)%nc(jkb, js, ibnd) + &
                                     CONJG(alphap(icart,ik)%nc(ikb,is,ibnd))* &
                                     alphap(jcart,ik)%nc(jkb, js, ibnd) + &
                                     CONJG(alphap(jcart,ik)%nc(ikb,is,ibnd))*&
                                     alphap(icart,ik)%nc(jkb, js, ibnd) )
                                   END DO
                                END DO
                             ELSE
                                dynwrk(na_icart,na_jcart) = &
                                  dynwrk(na_icart,na_jcart) + &
                                  deff (ih, jh, na)* wgg * &
                                  (CONJG(gammap(icart,jcart)%k(ikb,ibnd)) *&
                                   becp1(ik)%k (jkb, ibnd) + &
                                   CONJG (becp1(ik)%k (ikb, ibnd) ) * &
                                   gammap(icart,jcart)%k(jkb,ibnd) + &
                                   CONJG (alphap(icart,ik)%k(ikb, ibnd) ) * &
                                   alphap(jcart,ik)%k(jkb, ibnd) + &
                                   CONJG (alphap(jcart,ik)%k(ikb, ibnd) ) * &
                                   alphap(icart,ik)%k(jkb, ibnd) )
                             END IF
                          ENDDO
                       ENDDO
                    ENDDO
                 ENDDO
                 ijkb0 = ijkb0 + nh (nt)
              ENDIF
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !$acc end data
  !
  !   For true US pseudopotentials there is an additional term in the second
  !   derivative which is due to the change of the self consistent D part
  !   when the atom moves. We compute these terms in an additional routine
  !
  CALL stop_clock('dynusl')
  CALL addusdynmat (dynwrk)
  !
  CALL mp_sum ( dynwrk, inter_pool_comm )
  !
  !      do na = 1,nat
  !         do nb = 1,nat
  !           WRITE( stdout, '(2i3)') na,nb
  !            do icart = 1,3
  !              na_icart = 3*(na-1)+icart
  !               WRITE( stdout,'(6f13.8)')  &
  !                     (dynwrk(na_icart,3*(nb-1)+jcart), jcart=1,3)
  !            end do
  !         end do
  !      end do
  !      call stop_ph(.false.)
  !
  !  We rotate the dynamical matrix on the basis of patterns
  !

  CALL rotate_pattern_add(nat, u, dyn, dynwrk)

  IF (noncolin) THEN
     DEALLOCATE (deff_nc)
  ELSE
     DEALLOCATE (deff)
  END IF
#if defined(__CUDA)
  CALL deallocate_bec_type_acc(bectmp)
#endif
  DO icart=1,3
     DO jcart=1,3
        CALL deallocate_bec_type(gammap(icart,jcart))
     ENDDO
  ENDDO
  DEALLOCATE (aux1)
  DEALLOCATE (work2)
  DEALLOCATE (work1)
  DEALLOCATE (rhog)

  CALL stop_clock ('dynmat_us')
  RETURN
END SUBROUTINE dynmat_us
