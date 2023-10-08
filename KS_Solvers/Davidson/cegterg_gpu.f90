!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! NOTE (Ivan Carnimeo, May, 05th, 2022): 
!   cegterg and regterg have been ported to GPU with OpenACC, 
!   the previous CUF versions (cegterg_gpu and regterg_gpu) have been removed, 
!   and now cegterg and regterg are used for both CPU and GPU execution.
!   If you want to see the previous code checkout to commit: df3080b231c5daf52295c23501fbcaa9bfc4bfcc (on Thu Apr 21 06:18:02 2022 +0000)
!
#define ZERO ( 0.D0, 0.D0 )
#define ONE  ( 1.D0, 0.D0 )
!
#if ! defined(__CUDA)
! workaround for some old compilers that don't like CUDA fortran code
SUBROUTINE pcegterg_gpu( )
end SUBROUTINE pcegterg_gpu
#else
!----------------------------------------------------------------------------
!
!  Wrapper for subroutine with distributed matrixes (written by Carlo Cavazzoni)
!
!----------------------------------------------------------------------------
SUBROUTINE pcegterg_gpu(h_psi_gpu, s_psi_gpu, uspp, g_psi_gpu, &  
                        npw, npwx, nvec, nvecx, npol, evc_d, ethr, &
                        e_d, btype, notcnv, lrot, dav_iter , nhpsi )
  !----------------------------------------------------------------------------
  !
  ! ... iterative solution of the eigenvalue problem:
  !
  ! ... ( H - e S ) * evc = 0
  !
  ! ... where H is an hermitean operator, e is a real scalar,
  ! ... S is an uspp matrix, evc is a complex vector
  !
  USE util_param,       ONLY : DP, stdout
  USE mp_bands_util,    ONLY : intra_bgrp_comm, inter_bgrp_comm, root_bgrp_id, nbgrp, my_bgrp_id
  USE mp,               ONLY : mp_bcast, mp_root_sum, mp_sum, mp_barrier, &
                               mp_size, mp_type_free, mp_allgather
  USE device_memcpy_m,  ONLY : dev_memcpy, dev_memset
  !
  IMPLICIT NONE
  !
  include 'laxlib.fh'
  !
  INTEGER, INTENT(IN) :: npw, npwx, nvec, nvecx, npol
    ! dimension of the matrix to be diagonalized
    ! leading dimension of matrix evc, as declared in the calling pgm unit
    ! integer number of searched low-lying roots
    ! maximum dimension of the reduced basis set
    !    (the basis set is refreshed when its dimension would exceed nvecx)
    ! number of spin polarizations
  INTEGER, PARAMETER :: blocksize = 256
  INTEGER :: numblock
    ! chunking parameters
  COMPLEX(DP), DEVICE, INTENT(INOUT) :: evc_d(npwx*npol,nvec)
    !  evc   contains the  refined estimates of the eigenvectors
  REAL(DP), INTENT(IN) :: ethr
    ! energy threshold for convergence: root improvement is stopped,
    ! when two consecutive estimates of the root differ by less than ethr.
  LOGICAL, INTENT(IN) :: uspp
    ! if .FALSE. : S|psi> not needed
  INTEGER, INTENT(IN) :: btype(nvec)
    ! band type ( 1 = occupied, 0 = empty )
  LOGICAL, INTENT(IN) :: lrot
    ! .TRUE. if the wfc have already been rotated
  REAL(DP), DEVICE, INTENT(OUT) :: e_d(nvec)
    ! contains the estimated roots.
  INTEGER, INTENT(OUT) :: dav_iter, notcnv
    ! integer  number of iterations performed
    ! number of unconverged roots
  INTEGER, INTENT(OUT) :: nhpsi
    ! total number of indivitual hpsi
  !
  ! ... LOCAL variables
  !
  REAL(DP), ALLOCATABLE :: e(:)
  
  INTEGER, PARAMETER :: maxter = 20
    ! maximum number of iterations
  !
  INTEGER :: kter, nbase, np, kdim, kdmx, n, m, ipol, nb1, nbn
    ! counter on iterations
    ! dimension of the reduced basis
    ! counter on the reduced basis vectors
    ! do-loop counters
  INTEGER :: i, j, k, ierr
  REAL(DP), ALLOCATABLE :: ew(:)
  COMPLEX(DP), ALLOCATABLE :: hl(:,:), sl(:,:), vl(:,:), psi_w(:,:), spsi_w(:,:), hpsi_w(:,:) 
    ! Hamiltonian on the reduced basis
    ! S matrix on the reduced basis
    ! eigenvectors of the Hamiltonian
    ! eigenvalues of the reduced hamiltonian
    ! work space, contains psi
    ! the product of S and psi
    ! the product of H and psi
  LOGICAL, ALLOCATABLE :: conv(:)
    ! true if the root is converged
  REAL(DP) :: empty_ethr 
    ! threshold for empty bands
  INTEGER :: idesc(LAX_DESC_SIZE), idesc_old(LAX_DESC_SIZE)
  INTEGER, ALLOCATABLE :: irc_ip( : )
  INTEGER, ALLOCATABLE :: nrc_ip( : )
  INTEGER, ALLOCATABLE :: rank_ip( :, : )
    ! matrix distribution descriptors
  INTEGER :: nx
    ! maximum local block dimension
  LOGICAL :: la_proc
    ! flag to distinguish procs involved in linear algebra
  INTEGER, ALLOCATABLE :: notcnv_ip( : )
  INTEGER, ALLOCATABLE :: ic_notcnv( : )
  !
  INTEGER :: np_ortho(2), ortho_parent_comm
  LOGICAL :: do_distr_diag_inside_bgrp
  !
  REAL(DP), EXTERNAL :: myddot
  !
  EXTERNAL  h_psi_gpu, s_psi_gpu, g_psi_gpu
  INTEGER :: idx1, idx2
    ! h_psi(npwx,npw,nvec,psi,hpsi)
    !     calculates H|psi> 
    ! s_psi(npwx,npw,nvec,psi,spsi)
    !     calculates S|psi> (if needed)
    !     Vectors psi,hpsi,spsi are dimensioned (npwx,nvec)
    ! g_psi(npwx,npw,notcnv,psi,e)
    !    calculates (diag(h)-e)^-1 * psi, diagonal approx. to (h-e)^-1*psi
    !    the first nvec columns contain the trial eigenvectors
  !
  nhpsi = 0
  CALL start_clock( 'cegterg_gpu' )
  !
  CALL laxlib_getval( np_ortho = np_ortho, ortho_parent_comm = ortho_parent_comm, &
    do_distr_diag_inside_bgrp = do_distr_diag_inside_bgrp )
  !
  IF ( nvec > nvecx / 2 ) CALL errore( 'pcegterg', 'nvecx is too small', 1 )
  !
  ! ... threshold for empty bands
  !
  empty_ethr = MAX( ( ethr * 5.D0 ), 1.D-5 )
  !
  IF ( npol == 1 ) THEN
     !
     kdim = npw
     kdmx = npwx
     !
  ELSE
     !
     kdim = npwx*npol
     kdmx = npwx*npol
     !
  END IF
  !
  ! compute the number of chuncks
  numblock  = (npw+blocksize-1)/blocksize

  !
  ALLOCATE(  e( nvec ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' pcegterg ',' cannot allocate e (host) ', ABS(ierr) )
  !
  ! ... Initialize the matrix descriptor
  !
  ALLOCATE( ic_notcnv( np_ortho(2) ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' pcegterg ',' cannot allocate ic_notcnv ', ABS(ierr) )
  !
  ALLOCATE( notcnv_ip( np_ortho(2) ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' pcegterg ',' cannot allocate notcnv_ip ', ABS(ierr) )
  !
  ALLOCATE( irc_ip( np_ortho(1) ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' pcegterg ',' cannot allocate irc_ip ', ABS(ierr) )
  !
  ALLOCATE( nrc_ip( np_ortho(1) ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' pcegterg ',' cannot allocate nrc_ip ', ABS(ierr) )
  !
  ALLOCATE( rank_ip( np_ortho(1), np_ortho(2) ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' pcegterg ',' cannot allocate rank_ip ', ABS(ierr) )
  !
  CALL desc_init( nvec, nx, la_proc, idesc, rank_ip, irc_ip, nrc_ip )
  !
  IF( la_proc ) THEN
     !
     ! only procs involved in the diagonalization need to allocate local 
     ! matrix block.
     !
     ALLOCATE( vl( nx , nx ), STAT=ierr )
!$acc enter data create(vl) 
     IF( ierr /= 0 ) &
        CALL errore( ' pcegterg ',' cannot allocate vl_d ', ABS(ierr) )
     !
     ALLOCATE( sl( nx , nx ), STAT=ierr )
!$acc enter data create(sl)
     IF( ierr /= 0 ) &
        CALL errore( ' pcegterg ',' cannot allocate sl_d ', ABS(ierr) )
     !
     ALLOCATE( hl( nx , nx ) , STAT=ierr )
!$acc enter data create(hl) 
     IF( ierr /= 0 ) &
        CALL errore( ' pcegterg ',' cannot allocate hl_d ', ABS(ierr) )
     !
  ELSE
     !
     ALLOCATE( vl( 1 , 1 ), STAT=ierr )
!$acc enter data create(vl) 
     IF( ierr /= 0 ) &
        CALL errore( ' pcegterg ',' cannot allocate vl ', ABS(ierr) )
     !
     ALLOCATE( sl( 1 , 1 ), STAT=ierr )
     IF( ierr /= 0 ) &
        CALL errore( ' pcegterg ',' cannot allocate sl ', ABS(ierr) )
     !
     ALLOCATE( hl( 1 , 1 ), STAT=ierr )
     IF( ierr /= 0 ) &
        CALL errore( ' pcegterg ',' cannot allocate hl ', ABS(ierr) )
     !
  END IF
  !
  ALLOCATE( ew( nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' pcegterg ',' cannot allocate ew ', ABS(ierr) )
!$acc enter data create(ew) 
  !
  ALLOCATE( conv( nvec ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' pcegterg ',' cannot allocate conv ', ABS(ierr) )
  !
  notcnv = nvec
  nbase  = nvec
  conv   = .FALSE.
  !
  !CALL buffer%lock_buffer(psi_d, (/npwx*npol, nvecx/), ierr)
  ALLOCATE (psi_w(npwx*npol, nvecx)) 
!$acc enter data create(psi_w(1:npwx*npol, 1:nvecx)) 
  !
  !CALL buffer%lock_buffer(hpsi_d, (/npwx*npol, nvecx/), ierr)
  ALLOCATE(hpsi_w(npwx*npol, nvecx))
!$acc enter data create(hpsi_w(:npwx*npol, :nvecx)) 
  !
  !CALL buffer%lock_buffer(spsi_d, (/npwx*npol, nvecx/), ierr)
  ALLOCATE (spsi_w(npwx*npol, nvecx)) 
!$acc enter data create(spsi_w(1:npwx*npol,1:nvecx)) 
  ! 
!$acc host_data use_device(psi_w, spsi_w, hpsi_w) 
  CALL dev_memcpy(psi_w, evc_d, (/1, npwx*npol /), 1 , (/ 1, nvec /) )
  !
  ! ... hpsi contains h times the basis vectors
  !

  CALL h_psi_gpu( npwx, npw, nvec, psi_w, hpsi_w ) ; nhpsi = nhpsi + nvec
  !
  IF ( uspp ) CALL s_psi_gpu( npwx, npw, nvec, psi_w, spsi_w )
!$acc end host_data 
  !
  ! ... hl contains the projection of the hamiltonian onto the reduced
  ! ... space, vl contains the eigenvectors of hl. Remember hl, vl and sl
  ! ... are all distributed across processors, global replicated matrixes
  ! ... here are never allocated
  !
  CALL start_clock( 'cegterg_gpu:init' )
  CALL compute_distmat_gpu( hl, psi_w, hpsi_w )
  !
  IF ( uspp ) THEN
     !
     CALL compute_distmat_gpu( sl, psi_w, spsi_w )
     !
  ELSE
     !
     CALL compute_distmat_gpu( sl, psi_w, psi_w )
     !
  END IF
  CALL stop_clock( 'cegterg_gpu:init' )
  !
  IF ( lrot ) THEN
     !
     CALL set_e_from_h_gpu()
     e = e_d
     !
     CALL set_to_identity_gpu( vl, idesc )
     !
  ELSE
     !
     ! ... diagonalize the reduced hamiltonian
     !     Calling block parallel algorithm
     !
     CALL start_clock( 'cegterg_gpu:diag' )
     IF ( do_distr_diag_inside_bgrp ) THEN ! NB on output of pdiaghg ew and vl are the same across ortho_parent_comm
        ! only the first bgrp performs the diagonalization
#if defined (__SCALAPACK) 
!$acc update host(hl,sl,vl) 
        IF ( my_bgrp_id == root_bgrp_id ) CALL laxlib_pcdiaghg( nbase, hl, sl, nx, ew, vl, idesc) 
!$acc update device(vl) 
#else 
!$acc host_data use_device(vl,sl, hl) 
        IF( my_bgrp_id == root_bgrp_id ) CALL laxlib_pcdiaghg_gpu( nbase, hl, sl, nx, ew, vl, idesc, .false. )
!$acc end host_data
#endif 
!$acc host_data use_device(vl) 
        IF( nbgrp > 1 ) THEN ! results must be brodcast to the other band groups
           CALL mp_bcast( vl, root_bgrp_id, inter_bgrp_comm )
           CALL mp_bcast( ew, root_bgrp_id, inter_bgrp_comm )
        ENDIF
!$acc end host_data
     ELSE
#if defined (__SCALAPACK) 
!$acc update host(vl,sl,hl) 
        CALL laxlib_pcdiaghg( nbase, hl, sl, nx, ew, vl, idesc)
!$acc update device(vl) 
#else 
!$acc host_data use_device(vl,sl, hl) 
        CALL laxlib_pcdiaghg_gpu( nbase, hl, sl, nx, ew, vl, idesc, .false. )
!$acc end host_data
#endif 
     END IF
     CALL stop_clock( 'cegterg_gpu:diag' )
     !
     e(1:nvec) = ew(1:nvec)
     e_d(1:nvec) = ew(1:nvec)
     !
  END IF
  !
  ! ... iterate
  !
  iterate: DO kter = 1, maxter
     !
     dav_iter = kter
     !
     CALL start_clock( 'cegterg_gpu:update' )
     !
     CALL reorder_v()
     !
     nb1 = nbase + 1
     !
     ! ... expand the basis set with new basis vectors ( H - e*S )|psi> ...
     !
     !ew_d = ew ! NB: ew_d is needed by hpsi_dot_v_gpu
     !$acc update device(ew) 
     CALL hpsi_dot_v_gpu()
     !
     CALL stop_clock( 'cegterg_gpu:update' )
     !
     ! ... approximate inverse iteration
     !
!$acc host_data use_device(ew, psi_w) 
     CALL g_psi_gpu( npwx, npw, notcnv, npol, psi_w(1,nb1), ew(nb1) )
!$acc end host_data
     !
     ! ... "normalize" correction vectors psi(:,nb1:nbase+notcnv) in 
     ! ... order to improve numerical stability of subspace diagonalization 
     ! ... (cdiaghg) ew is used as work array :
     !
     ! ...         ew = <psi_i|psi_i>,  i = nbase + 1, nbase + notcnv
     !
!$acc kernels  
     DO n = 1, notcnv
        !
        nbn = nbase + n
        !
        IF ( npol == 1 ) THEN
           !
           ew(n) = dot_product( psi_w(1:npw,nbn), psi_w(1:npw,nbn))
           !
        ELSE
           !
           ew(n) = dot_product( psi_w(1:npw,nbn), psi_w(1:npw,nbn)) + &
                   dot_product( psi_w(npwx+1:npwx+npw,nbn), psi_w(npwx+1:npwx:npwx+npw,nbn))
           !
        END IF
        !
     END DO
!$acc end kernels 
     !
!$acc host_data use_device(ew) 
     CALL mp_sum( ew( 1:notcnv ), intra_bgrp_comm )
!$acc end host_data
     !

!$acc kernels present(ew,psi_w) 
     DO i = 1, notcnv
        idx2 = nbase+i
        psi_w(1:npol*npwx, idx2) = psi_w(1:npol*npwx,idx2)/SQRT( ew(i) )
     END DO
!$acc end kernels
     !
     ! ... here compute the hpsi and spsi of the new functions
     !
!$acc host_data use_device(psi_w, spsi_w, hpsi_w) 
     CALL h_psi_gpu( npwx, npw, notcnv, psi_w(1,nb1), hpsi_w(1,nb1) ) ; nhpsi = nhpsi + notcnv
     !
     IF ( uspp ) CALL s_psi_gpu( npwx, npw, notcnv, psi_w(1,nb1), spsi_w(1,nb1) )
!$acc end host_data
     !
     ! ... update the reduced hamiltonian
     !
     CALL start_clock( 'cegterg_gpu:overlap' )
     !
     ! we need to save the old descriptor in order to redistribute matrices 
     !
     idesc_old = idesc
     !
     ! ... RE-Initialize the matrix descriptor
     !
     CALL desc_init( nbase+notcnv, nx, la_proc, idesc, rank_ip, irc_ip, nrc_ip )
     !
     IF( la_proc ) THEN

        !  redistribute hl and sl (see dsqmred), since the dimension of the subspace has changed
        !

!FIXME data should be copied directly from hl device to vl host 
!$acc kernels 
        vl = hl 
!$acc end kernels 
!$acc update host(vl) 
!
!$acc exit data finalize delete(hl) 
        DEALLOCATE( hl)
        !
        ALLOCATE( hl( nx , nx ), STAT=ierr )
        IF( ierr /= 0 ) &
           CALL errore( ' pcegterg ',' cannot allocate hl ', ABS(ierr) )
!$acc enter data create(hl) 
        !
        CALL laxlib_zsqmred( nbase, vl, idesc_old(LAX_DESC_NRCX), idesc_old, nbase+notcnv, hl, nx, idesc )
!$acc update device(hl) 
!FIXME data should be copied directly from hl device to vl host  
!$acc kernels 
        vl = sl
!$acc end kernels 
!$acc update host(vl) 
!
!$acc exit data finalize delete(sl) 
        DEALLOCATE( sl) 
        ALLOCATE( sl( nx , nx ), STAT=ierr )
!$acc enter data create(sl) 
        IF( ierr /= 0 ) &
           CALL errore( ' pcegterg ',' cannot allocate sl ', ABS(ierr) )

        CALL laxlib_zsqmred( nbase, vl, idesc_old(LAX_DESC_NRCX), idesc_old, nbase+notcnv, sl, nx, idesc )
!$acc update device(sl) 
!$acc exit data finalize delete(vl) 
        DEALLOCATE( vl )
        ALLOCATE( vl( nx , nx ), STAT=ierr )
!$acc enter data create(vl) 
        IF( ierr /= 0 ) &
           CALL errore( ' pcegterg ',' cannot allocate vl ', ABS(ierr) )
     END IF
     !
     !
     CALL update_distmat_gpu( hl, psi_w, hpsi_w )
     !
     IF ( uspp ) THEN
        !
        CALL update_distmat_gpu( sl, psi_w, spsi_w )
        !
     ELSE
        !
        CALL update_distmat_gpu( sl, psi_w, psi_w )
        !
     END IF
     !
     CALL stop_clock( 'cegterg_gpu:overlap' )
     !
     nbase = nbase + notcnv
     !
     ! ... diagonalize the reduced hamiltonian
     !     Call block parallel algorithm
     !
     CALL start_clock( 'cegterg_gpu:diag' )
     IF ( do_distr_diag_inside_bgrp ) THEN ! NB on output of pdiaghg ew and vl are the same across ortho_parent_comm
        ! only the first bgrp performs the diagonalization
#if defined (__SCALAPACK) 
!$acc update host(hl,sl,vl) 
        IF( my_bgrp_id == root_bgrp_id ) CALL laxlib_pcdiaghg( nbase, hl, sl, nx, ew, vl, idesc) 
!$acc update device(vl) 
#else 
!$acc host_data use_device(hl, sl, vl) 
        IF( my_bgrp_id == root_bgrp_id ) CALL laxlib_pcdiaghg_gpu( nbase, hl, sl, nx, ew, vl, idesc, .false. )
!$acc end host_data
#endif 
!$acc host_data use_device(vl) 
        IF( nbgrp > 1 ) THEN ! results must be brodcast to the other band groups
           CALL mp_bcast( vl, root_bgrp_id, inter_bgrp_comm )
           CALL mp_bcast( ew, root_bgrp_id, inter_bgrp_comm )
        ENDIF
!$acc end host_data 
     ELSE
#if defined(__SCALAPACK) 
!$acc update host(hl, vl, sl) 
        CALL laxlib_pcdiaghg( nbase, hl, sl, nx, ew, vl, idesc)
!$acc update device(vl) 
#else  
!$acc host_data use_device(hl, sl, vl) 
        CALL laxlib_pcdiaghg_gpu( nbase, hl, sl, nx, ew, vl, idesc, .false. )
!$acc end host_data 
#endif
     END IF
     CALL stop_clock( 'cegterg_gpu:diag' )
     !
     ! ... test for convergence
     !
     WHERE( btype(1:nvec) == 1 )
        !
        conv(1:nvec) = ( ( ABS( ew(1:nvec) - e(1:nvec) ) < ethr ) )
        !
     ELSEWHERE
        !
        conv(1:nvec) = ( ( ABS( ew(1:nvec) - e(1:nvec) ) < empty_ethr ) )
        !
     END WHERE
     ! ... next line useful for band parallelization of exact exchange
     IF ( nbgrp > 1 ) CALL mp_bcast(conv,root_bgrp_id,inter_bgrp_comm)
     !
     notcnv = COUNT( .NOT. conv(:) )
     !
     e(1:nvec) = ew(1:nvec)
!$acc update device(ew) 
!$acc kernels present(ew) 
     e_d(1:nvec) = ew(1:nvec)
!$acc end kernels
     !ew_d = ew
     !
     ! ... if overall convergence has been achieved, or the dimension of
     ! ... the reduced basis set is becoming too large, or in any case if
     ! ... we are at the last iteration refresh the basis set. i.e. replace
     ! ... the first nvec elements with the current estimate of the
     ! ... eigenvectors;  set the basis dimension to nvec.
     !
     IF ( notcnv == 0 .OR. nbase+notcnv > nvecx .OR. dav_iter == maxter ) THEN
        !
        CALL start_clock( 'cegterg_gpu:last' )
        !
        CALL refresh_evc_gpu()
        !
        IF ( notcnv == 0 ) THEN
           !
           ! ... all roots converged: return
           !
           CALL stop_clock( 'cegterg_gpu:last' )
           !
           EXIT iterate
           !
        ELSE IF ( dav_iter == maxter ) THEN
           !
           ! ... last iteration, some roots not converged: return
           !
           CALL stop_clock( 'cegterg_gpu:last' )
           !
           EXIT iterate
           !
        END IF
        !
        ! ... refresh psi, H*psi and S*psi
        !
!$acc host_data use_device(psi_w) 
        CALL dev_memcpy(psi_w, evc_d, (/1, npwx*npol /), 1 , (/ 1, nvec /), 1) ! need if refresh_evc_gpu
!$acc end host_data
        !
        IF ( uspp ) THEN
           !
           CALL refresh_spsi_gpu(spsi_w, psi_w) 
           ! 
        END IF
        !
        CALL refresh_hpsi_gpu()
        !
        ! ... refresh the reduced hamiltonian
        !
        nbase = nvec
        !
        CALL desc_init( nvec, nx, la_proc, idesc, rank_ip, irc_ip, nrc_ip )
        !
        IF( la_proc ) THEN
           !
           ! note that nx has been changed by desc_init
           ! we need to re-alloc with the new size.
           !
!$acc exit data finalize delete(vl, hl, sl) 
           DEALLOCATE( vl, hl, sl )
           ALLOCATE( vl( nx, nx ), STAT=ierr )
!$acc enter data create(vl) 
           IF( ierr /= 0 ) &
              CALL errore( ' pcegterg ',' cannot allocate vl ', ABS(ierr) )
           ALLOCATE( hl( nx, nx ), STAT=ierr )
!$acc enter data create(hl) 
           IF( ierr /= 0 ) &
              CALL errore( ' pcegterg ',' cannot allocate hl ', ABS(ierr) )
           ALLOCATE( sl( nx, nx ), STAT=ierr )
!$acc enter data create(sl) 
           IF( ierr /= 0 ) &
              CALL errore( ' pcegterg ',' cannot allocate sl ', ABS(ierr) )
           !
        END IF
        !
        CALL set_h_from_e_gpu( )
        !
        CALL set_to_identity_gpu( vl, idesc )
        CALL set_to_identity_gpu( sl, idesc )
        !
        CALL stop_clock( 'cegterg_gpu:last' )
        !
     END IF
     !
  END DO iterate
  !
!$acc exit data finalize delete(vl,hl, sl,ew, psi_w, spsi_w,hpsi_w) 
  DEALLOCATE( vl, hl, sl )
  !
  DEALLOCATE( rank_ip )
  DEALLOCATE( ic_notcnv )
  DEALLOCATE( irc_ip )
  DEALLOCATE( nrc_ip )
  DEALLOCATE( notcnv_ip )
  DEALLOCATE( conv )
  DEALLOCATE( ew )
  DEALLOCATE( e )
  DEALLOCATE(psi_w) 
  DEALLOCATE(spsi_w)
  DEALLOCATE(hpsi_w)
 
  !
  CALL stop_clock( 'cegterg_gpu' )
  !
  RETURN
  !
  !
CONTAINS
  !
  SUBROUTINE set_to_identity_gpu( distmat, idesc )
     IMPLICIT NONE
     INTEGER, INTENT(IN)  :: idesc(LAX_DESC_SIZE)
     COMPLEX(DP), INTENT(OUT) :: distmat(:,:)
     ! 
     INTEGER :: i,nd
!$acc data present(distmat) 
!$acc kernels
     distmat = ( 0_DP , 0_DP )
!$acc end kernels
     IF( idesc(LAX_DESC_MYC) == idesc(LAX_DESC_MYR) .AND. idesc(LAX_DESC_ACTIVE_NODE) > 0 ) THEN
        nd =  idesc(LAX_DESC_NC)
        !$acc parallel loop 
        DO i = 1, nd
           distmat( i, i ) = ( 1_DP , 0_DP )
        END DO
     END IF 
!$acc end data
     RETURN
  END SUBROUTINE set_to_identity_gpu
   !
  SUBROUTINE reorder_v()
     !
     IMPLICIT NONE
     INTEGER :: ipc
     INTEGER :: nc, ic
     INTEGER :: nl, npl
     !
     np = 0
     !
     notcnv_ip = 0
     !
     n = 0
     !
     DO ipc = 1, idesc(LAX_DESC_NPC)
        !
        nc = nrc_ip( ipc )
        ic = irc_ip( ipc )
        !
        npl = 0
        !
        IF( ic <= nvec ) THEN
           !
           DO nl = 1, min( nvec - ic + 1, nc )
              !
              n  = n  + 1
              !
              IF ( .NOT. conv(n) ) THEN
                 !
                 ! ... this root not yet converged ... 
                 !
                 np  = np  + 1
                 npl = npl + 1
                 IF( npl == 1 ) ic_notcnv( ipc ) = np
                 !
                 ! ... reorder eigenvectors so that coefficients for unconverged
                 ! ... roots come first. This allows to use quick matrix-matrix 
                 ! ... multiplications to set a new basis vector (see below)
                 !
                 notcnv_ip( ipc ) = notcnv_ip( ipc ) + 1
                 !
                 IF ( npl /= nl ) THEN
                    IF( la_proc .AND. idesc(LAX_DESC_MYC) == ipc-1 ) THEN
!$acc kernels 
                       vl( :, npl) = vl( :, nl )
!$acc end kernels 
                    END IF
                 END IF
                 !
                 ! ... for use in g_psi
                 !
                 ew(nbase+np) = e(n)
                 !   
              END IF
              !
           END DO
           !
        END IF
        !
     END DO
     !
  END SUBROUTINE reorder_v
  !
  !
  SUBROUTINE hpsi_dot_v_gpu()
     !
     use cublas
     !
     IMPLICIT NONE
     INTEGER :: ipc, ipr
     INTEGER :: nr, ir, ic, notcl, root, np, ipol, ib
     COMPLEX(DP) :: beta
     !
     COMPLEX(DP),         ALLOCATABLE :: vtmp(:,:), ptmp(:,:) 
     COMPLEX(DP)                      :: ps1, ps2
     INTEGER                          :: idx1, idx2, offsvec, npol_mio, npwx_mio 
     !
     ALLOCATE( vtmp( nx, nx ) )
     ALLOCATE( ptmp( npwx*npol, nx ) )
!$acc enter data create(vtmp, ptmp) 
!$acc kernels present(vtmp, ptmp)  
     vtmp = ZERO 
     ptmp = ZERO
!$acc end kernels


     DO ipc = 1, idesc(LAX_DESC_NPC)
        !
        IF( notcnv_ip( ipc ) > 0 ) THEN

           notcl = notcnv_ip( ipc )
           ic    = ic_notcnv( ipc )

           beta = ZERO

           DO ipr = 1, idesc(LAX_DESC_NPR)
              !
              nr = nrc_ip( ipr )
              ir = irc_ip( ipr )
              !
              root = rank_ip( ipr, ipc )
! FIXME

!$acc host_data use_device(vtmp, ptmp,psi_w, spsi_w, hpsi_w, vl) 
              IF( ipr-1 == idesc(LAX_DESC_MYR) .AND. ipc-1 == idesc(LAX_DESC_MYC) .AND. la_proc ) & 
                 CALL dev_memcpy(vtmp(1:,1:notcl),vl(1:,1:notcl),[1,nx],1,[1,notcl]) 
              !
              CALL mp_bcast( vtmp(:,1:notcl), root, ortho_parent_comm )
              !
              !
              IF ( uspp ) THEN
                 !
                 CALL ZGEMM( 'N', 'N', kdim, notcl, nr, ONE, &
                    spsi_w(1, ir), kdmx, vtmp, nx, beta, psi_w(1,nb1+ic-1), kdmx )
                 !
              ELSE
                 !
                 CALL ZGEMM( 'N', 'N', kdim, notcl, nr, ONE, &
                    psi_w(1, ir), kdmx, vtmp, nx, beta, psi_w(1,nb1+ic-1), kdmx )
                 !
              END IF
              !
              CALL ZGEMM( 'N', 'N', kdim, notcl, nr, ONE, &
                      hpsi_w(1, ir), kdmx, vtmp, nx, beta, ptmp, kdmx )
              !
              beta = ONE
              !
!$acc end host_data
           END DO
           !
           offsvec = nbase + ic -1  
           npwx_mio = npwx*npol 
           !$acc kernels present(ptmp,ew, psi_w) 
           !!$acc parallel present(ptmp, ew,psi_w) 
           !!$acc loop gang private(ps1, idx2)  
           DO np = offsvec+1, offsvec + notcl
              idx2 = np - offsvec 
              ps1 = ew(np) 
              !!$acc loop private(ps2) vector 
              DO k = 1, npwx_mio
                   ps2 = ps1 * psi_w(k,np) 
                   psi_w(k, np) = ptmp(k, idx2) - ps2
              END DO
           END DO
           !$acc end kernels
           !!$acc end parallel 
           !
           ! clean up garbage if there is any
!$acc kernels present(psi_w) 
           IF (npw < npwx) psi_w(npw+1:npwx,nbase+ic:nbase+notcl+ic-1) = ZERO
           IF (npol == 2)  psi_w(npwx+npw+1:2*npwx,nbase+ic:nbase+notcl+ic-1) = ZERO
!$acc end kernels 
           !
        END IF
        !
     END DO


!$acc exit data delete(vtmp, ptmp) 
     DEALLOCATE( vtmp)
     DEALLOCATE( ptmp)

     RETURN
  END SUBROUTINE hpsi_dot_v_gpu
  !
  !
  SUBROUTINE refresh_evc_gpu( )
     !
     use cublas
     !
     IMPLICIT NONE
     INTEGER :: ipc, ipr
     INTEGER :: nr, nc, ir, ic, root
     COMPLEX(DP) :: beta
     !
     COMPLEX(DP), ALLOCATABLE :: work(:,:)
     !
     ALLOCATE( work( nx, nx ) )
!$acc data present (vl, psi_w) create(work) 
!$acc host_data use_device(vl, psi_w, work) 
!$acc kernels 
     work = ZERO
!$acc end kernels
     !
     DO ipc = 1, idesc(LAX_DESC_NPC)
        !
        nc = nrc_ip( ipc )
        ic = irc_ip( ipc )
        !
        IF( ic <= nvec ) THEN
           !
           nc = min( nc, nvec - ic + 1 )
           !
           beta = ZERO

           DO ipr = 1, idesc(LAX_DESC_NPR)
              !
              nr = nrc_ip( ipr )
              ir = irc_ip( ipr )
              !
              root = rank_ip( ipr, ipc )

              IF( ipr-1 == idesc(LAX_DESC_MYR) .AND. ipc-1 == idesc(LAX_DESC_MYC) .AND. la_proc ) THEN
                 !
                 !  this proc sends his block
                 !
                 CALL mp_bcast( vl(:,1:nc), root, ortho_parent_comm )
                 !
!!$acc kernels present(work , vl) 
!                 work (1:nx,1:nc) = vl(1:nx,1:nc) ! FIXME!
!!$acc end kernels 
                 CALL dev_memcpy(work(1:,1:), vl(1:,1:),[1,nx],1,[1,nc])  
                 !
                 CALL ZGEMM( 'N', 'N', kdim, nc, nr, ONE, &
                          psi_w(1,ir), kdmx, work, nx, beta, evc_d(1,ic), kdmx )
                 !
              ELSE
                 !
                 !  all other procs receive
                 !
                 CALL mp_bcast( work(:,1:nc), root, ortho_parent_comm )
                 !
                 CALL ZGEMM( 'N', 'N', kdim, nc, nr, ONE, &
                          psi_w(1,ir), kdmx, work, nx, beta, evc_d(1,ic), kdmx )
                 !
              END IF
              !
              beta = ONE
              !
           END DO
           !
        END IF
        !
     END DO
     !
!$acc end host_data
!$acc end data 
     DEALLOCATE( work)
     !
     RETURN
  END SUBROUTINE refresh_evc_gpu
  !
  !
  SUBROUTINE refresh_spsi_gpu( spsi_out, psi_in)
     !
     use cublas
     !
     IMPLICIT NONE
     COMPLEX(DP) :: spsi_out(:,:), psi_in(:,:) 
     INTEGER :: ipc, ipr
     INTEGER :: nr, nc, ir, ic, root, nvec_, kdmx_, nvecx_, idx1
     COMPLEX(DP) :: beta
     !
     COMPLEX(DP), ALLOCATABLE :: work(:,:)
     !
     ALLOCATE( work( nx, nx ) )
     nvec_ = nvec
     nvecx_ = nvecx
     kdmx_ = npol*npwx 
     !
!$acc data present(spsi_out(1:kdmx_,1:nvecx), psi_in(1:kdmx_,1:nvecx_)) create(work(:nx,:nx))
!$acc kernels present(work(:nx,:nx)) 
     work = ZERO
!$acc end kernels  
     DO ipc = 1, idesc(LAX_DESC_NPC)
        !
        nc = nrc_ip( ipc )
        ic = irc_ip( ipc )
        !
        IF( ic <= nvec ) THEN
           !
           nc = min( nc, nvec - ic + 1 )
           !
           beta = ZERO
           !
           DO ipr = 1, idesc(LAX_DESC_NPR)
              !
              nr = nrc_ip( ipr )
              ir = irc_ip( ipr )
              !
              root = rank_ip( ipr, ipc )

              IF( ipr-1 == idesc(LAX_DESC_MYR) .AND. ipc-1 == idesc(LAX_DESC_MYC) .AND. la_proc ) THEN
                 !
                 !  this proc sends his block
                 !
!$acc host_data use_device(vl) 
                 CALL mp_bcast( vl(:,1:nc), root, ortho_parent_comm )
!$acc end host_data
                 !
!$acc kernels present(vl, work) 
                 work(:,1:nc) = vl(:,1:nc)
!$acc end kernels
                 !
!$acc host_data use_device(psi_in,spsi_out, work) 
                 CALL ZGEMM( 'N', 'N', kdim, nc, nr, ONE, &
                          spsi_out(1,ir), kdmx, work, nx, beta, psi_in(1,nvec+ic), kdmx )
!$acc end host_data 
                 !
              ELSE
                 !
                 !  all other procs receive
                 !
!$acc host_data use_device(work, psi_in, spsi_out) 
                 CALL mp_bcast( work(:,1:nc), root, ortho_parent_comm )
                 !
                 CALL ZGEMM( 'N', 'N', kdim, nc, nr, ONE, &
                          spsi_out(1,ir), kdmx, work, nx, beta, psi_in(1,nvec+ic), kdmx )
!$acc end host_data 
                 !
              END IF
              !
              beta = ONE
              !
           END DO
           !
        END IF
        !
     END DO
     ! 
     !call start_clock("final_copy") 
     !!$cuf kernel do(2) <<<*,*>>>
     !DO j = 1, nvec
     !   DO i = 1, npwx*npol
     !      print *, "CIAO ", i, j
     !      spsi_d(i,j) = psi_w(i,nvec+j)
     !   END DO
     !END DO
     ! 
!!$acc kernels present(spsi_out(1:kdmx_, 1:nvecx_), psi_in(1:kdmx_,1:nvecx_))  
!     DO j =nvec_ + 1, 2 * nvec_ 
!        idx1 = j - nvec_ 
!        spsi_out(1:kdmx_,idx1) = psi_in(1:kdmx_, j)  
!     END DO 
!!$acc end kernels 
!$acc host_data use_device(spsi_out, psi_in) 
     call dev_memcpy(spsi_out(1:,1:), psi_in(1:,nvec_ + 1:),[1,kdmx_],1,[1,nvec_],1)
!$acc end host_data 
!$acc end data 
     DEALLOCATE( work)
     RETURN
  END SUBROUTINE refresh_spsi_gpu
  !
  !
  SUBROUTINE refresh_hpsi_gpu( )
     !
     use cublas
     !
     IMPLICIT NONE
     INTEGER :: ipc, ipr
     INTEGER :: nr, nc, ir, ic, root
     COMPLEX(DP) :: beta
     !
     COMPLEX(DP), ALLOCATABLE :: work(:,:)
!
INTEGER :: i, j
     !
     ALLOCATE( work( nx, nx ) )
!$acc data create(work(:nx,:nx)) present( hpsi_w(:npol*npwx,:nvecx), psi_w(:npol*npwx,:nvecx))
!$acc kernels
     work  = ZERO
!$acc end kernels 
     !
     DO ipc = 1, idesc(LAX_DESC_NPC)
        !
        nc = nrc_ip( ipc )
        ic = irc_ip( ipc )
        !
        IF( ic <= nvec ) THEN
           !
           nc = min( nc, nvec - ic + 1 )
           !
           beta = ZERO
           !
           DO ipr = 1, idesc(LAX_DESC_NPR)
              !
              nr = nrc_ip( ipr )
              ir = irc_ip( ipr )
              !
              root = rank_ip( ipr, ipc )

              IF( ipr-1 == idesc(LAX_DESC_MYR) .AND. ipc-1 == idesc(LAX_DESC_MYC) .AND. la_proc ) THEN
                 !
                 !  this proc sends his block
                 !
!$acc host_data use_device(vl) 
                 CALL mp_bcast( vl(:,1:nc), root, ortho_parent_comm )
!$acc end host_data
                 !
!$acc kernels present(work, vl) 
                 work(:,1:nc) = vl(:,1:nc)
!$acc end kernels 
                 !
!$acc host_data use_device(psi_w,hpsi_w, work) 
                 CALL ZGEMM( 'N', 'N', kdim, nc, nr, ONE, &
                          hpsi_w(1,ir), kdmx, work, nx, beta, psi_w(1,nvec+ic), kdmx )
!$acc end host_data
                 !
              ELSE
                 !
                 !  all other procs receive
                 !
!$acc host_data use_device(psi_w, hpsi_w, work) 
                 CALL mp_bcast( work(:,1:nc), root, ortho_parent_comm )
                 !
                 CALL ZGEMM( 'N', 'N', kdim, nc, nr, ONE, &
                          hpsi_w(1,ir), kdmx, work, nx, beta, psi_w(1,nvec+ic), kdmx )
!$acc end host_data
                 !
              END IF
              !
              beta = ONE
              !
           END DO
           !
        END IF
        !
     END DO
     !
     !!$cuf kernel do(2) <<<*,*>>>
     !$acc kernels present(psi_w,hpsi_w) 
     DO j = 1, nvec
           hpsi_w(1:npwx*npol,j) = psi_w(1:npwx*npol,nvec+j)
     END DO
     !$acc end kernels 
     !
!$acc end data 
     DEALLOCATE( work)
     !
     RETURN
  END SUBROUTINE refresh_hpsi_gpu
  !
  !
  SUBROUTINE compute_distmat_gpu( dm, v, w )
     !
     !  This subroutine compute <vi|wj> and store the
     !  result in distributed matrix dm
     !
     use cublas
     !
     IMPLICIT NONE
     INTEGER :: ipc, ipr
     INTEGER :: nr, nc, ir, ic, root
     COMPLEX(DP), INTENT(OUT) :: dm( :, : )
     COMPLEX(DP), INTENT(IN) :: v(:,:), w(:,:)
     COMPLEX(DP), ALLOCATABLE :: work(:,:)
     !
     ALLOCATE( work( nx, nx ) )
!$acc data present(dm, v, w) create(work) 
!$acc host_data use_device(dm, v, w, work) 
!$acc kernels 
     work = ZERO
!$acc end kernels 
     !
     !
     !  Only upper triangle is computed, then the matrix is hermitianized
     !
     DO ipc = 1, idesc(LAX_DESC_NPC) !  loop on column procs
        !
        nc = nrc_ip( ipc )
        ic = irc_ip( ipc )
        !
        DO ipr = 1, ipc ! idesc(LAX_DESC_NPR) ! ipc ! use symmetry for the loop on row procs
           !
           nr = nrc_ip( ipr )
           ir = irc_ip( ipr )
           !
           !  rank of the processor for which this block (ipr,ipc) is destinated
           !
           root = rank_ip( ipr, ipc )

           ! use blas subs. on the matrix block
           CALL ZGEMM( 'C', 'N', nr, nc, kdim, ONE , &
                       v(1,ir), kdmx, w(1,ic), kdmx, ZERO, work, nx )
           !
           ! accumulate result on dm of root proc.
           !
           CALL mp_root_sum( work, dm, root, ortho_parent_comm )
        END DO
        !
     END DO
     if (ortho_parent_comm.ne.intra_bgrp_comm .and. nbgrp > 1) then
             !$acc  kernels 
             do ir = 1, nx
             do ic = 1, nx
             dm(ir,ic) = dm(ir,ic)/nbgrp
             enddo
             enddo
             !$acc end kernels 
     endif
     !
     !  The matrix is hermitianized using upper triangle
     !
     CALL laxlib_zsqmher( nbase, dm, nx, idesc )
!$acc end host_data
!$acc end data 
     !
     DEALLOCATE( work )
     !
     RETURN
  END SUBROUTINE compute_distmat_gpu
  !
  !
  SUBROUTINE update_distmat_gpu( dm, v, w )
     !
     use cublas
     !
     IMPLICIT NONE
     INTEGER :: ipc, ipr
     INTEGER :: nr, nc, ir, ic, root, icc, ii
     COMPLEX(DP), INTENT(INOUT) :: dm( :, : )
     COMPLEX(DP), INTENT(IN) :: v(:,:), w(:,:)
     !
     COMPLEX(DP), ALLOCATABLE :: work(:,:)
     !
     ALLOCATE( work( nx, nx ) )
!$acc data present(dm, v, w) create(work) 
!$acc host_data use_device(dm,v,w,work) 
!$acc kernels 
     work = ZERO
!$acc end kernels 
     !
     !
     DO ipc = 1, idesc(LAX_DESC_NPC)
        !
        nc = nrc_ip( ipc )
        ic = irc_ip( ipc )
        !
        IF( ic+nc-1 >= nb1 ) THEN
           !
           nc = MIN( nc, ic+nc-1 - nb1 + 1 )
           IF( ic >= nb1 ) THEN
              ii = ic
              icc = 1
           ELSE
              ii = nb1
              icc = nb1-ic+1
           END IF
           !
           ! icc to nc is the local index of the unconverged bands
           ! ii is the global index of the first unconverged bands
           !
           DO ipr = 1, ipc ! idesc(LAX_DESC_NPR) use symmetry
              !
              nr = nrc_ip( ipr )
              ir = irc_ip( ipr )
              !
              root = rank_ip( ipr, ipc )
              !
              CALL ZGEMM( 'C', 'N', nr, nc, kdim, ONE, v(1, ir), &
                          kdmx, w(1,ii), kdmx, ZERO, work, nx )
              !
              !
!$acc kernels 
              IF (ortho_parent_comm.ne.intra_bgrp_comm .and. nbgrp > 1) work  = work / nbgrp
!$acc end kernels 
              !
              IF(  (idesc(LAX_DESC_ACTIVE_NODE) > 0) .AND. &
                   (ipr-1 == idesc(LAX_DESC_MYR)) .AND. (ipc-1 == idesc(LAX_DESC_MYC)) ) THEN
                 CALL mp_root_sum( work(:,1:nc), dm(:,icc:icc+nc-1), root, ortho_parent_comm )
              ELSE
                 CALL mp_root_sum( work(:,1:nc), dm, root, ortho_parent_comm )
              END IF

           END DO
           !
        END IF
        !
     END DO
     !
     CALL laxlib_zsqmher( nbase+notcnv, dm, nx, idesc )
     !
!$acc end host_data 
!$acc end data
     DEALLOCATE( work)
     RETURN
  END SUBROUTINE update_distmat_gpu
  !
  !
  SUBROUTINE set_e_from_h_gpu()
     IMPLICIT NONE
     INTEGER :: nc, ic, i
     e_d(1:nbase) = 0_DP
     IF( idesc(LAX_DESC_MYC) == idesc(LAX_DESC_MYR) .AND. la_proc ) THEN
        nc = idesc(LAX_DESC_NC)
        ic = idesc(LAX_DESC_IC)
        !$acc kernels loop present(e_d, hl) 
        DO i = 1, nc
           e_d( i + ic - 1 ) = REAL( hl( i, i ) )
        END DO
        !$acc end kernels loop 
     END IF
     CALL mp_sum( e_d(1:nbase), ortho_parent_comm )
     RETURN
  END SUBROUTINE set_e_from_h_gpu
  !
  SUBROUTINE set_h_from_e_gpu()
     IMPLICIT NONE
     INTEGER :: nc, ic, i
     IF( la_proc ) THEN
!$acc kernels present(hl, e_d) 
        hl  = ZERO
        IF( idesc(LAX_DESC_MYC) == idesc(LAX_DESC_MYR) ) THEN
           nc = idesc(LAX_DESC_NC)
           ic = idesc(LAX_DESC_IC)
           !$acc loop 
           DO i = 1, nc
              hl(i,i) = CMPLX( e_d( i + ic - 1 ), 0_DP ,kind=DP)
           END DO
        END IF
!$acc end kernels
     END IF
     RETURN
  END SUBROUTINE set_h_from_e_gpu
  !
	END SUBROUTINE pcegterg_gpu
#endif
