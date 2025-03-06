!
! Copyright (C) 2001-2024 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE orthoatwfc (orthogonalize_wfc)
  !-----------------------------------------------------------------------
  !
  !! This routine calculates atomic wavefunctions,
  !! orthogonalizes them if "orthogonalize_wfc" is .true.,
  !! saves them into buffer "iunsat" (that must be opened on input)
  !! Useful for options "wannier" and "one_atom_occupations"
  !
  USE kinds,            ONLY : DP
  USE buffers,          ONLY : save_buffer
  USE io_global,        ONLY : stdout
  USE io_files,         ONLY : iunsat, nwordatwfc
  USE ions_base,        ONLY : nat
  USE basis,            ONLY : natomwfc
  USE klist,            ONLY : nks, xk, ngk, igk_k
  USE wvfct,            ONLY : npwx
  USE uspp,             ONLY : nkb, vkb
  USE becmod,           ONLY : allocate_bec_type_acc, deallocate_bec_type_acc, &
                               bec_type, becp, calbec
  USE control_flags,    ONLY : gamma_only, use_gpu, offload_type, offload_type2
  USE noncollin_module, ONLY : noncolin, npol
  USE uspp_init,        ONLY : init_us_2
  IMPLICIT NONE
  !
  LOGICAL, INTENT(in) :: orthogonalize_wfc
  !
  INTEGER :: ik, ibnd, info, i, j, k, na, nb, nt, isym, n, ntemp, m, &
       l, lm, ltot, ntot, ipol, npw
  ! ik: the k point under consideration
  ! ibnd: counter on bands
  LOGICAL :: normalize_only = .FALSE.
  COMPLEX(DP) , ALLOCATABLE ::  wfcatom (:,:)
  COMPLEX(DP) , ALLOCATABLE :: swfcatom (:,:)
  
  normalize_only=.FALSE.
  ALLOCATE(  wfcatom(npwx*npol,natomwfc) )
  ALLOCATE( swfcatom(npwx*npol,natomwfc) )
  !$acc enter data create(wfcatom, swfcatom)
#if defined(__OPENMP_GPU)
  !$omp target data map(alloc:wfcatom, swfcatom)
#endif
  !
  ! Allocate the array becp = <beta|wfcatom>
  CALL allocate_bec_type_acc (nkb,natomwfc, becp) 
  
  DO ik = 1, nks
     
     IF (noncolin) THEN
       CALL atomic_wfc_nc_updown (ik, wfcatom)
     ELSE
       CALL atomic_wfc (ik, wfcatom)
     ENDIF
#if defined(__OPENMP_GPU)
     !$omp target update to(wfcatom)
#endif
     npw = ngk (ik)
     !
     CALL init_us_2 (npw, igk_k(1,ik), xk (1, ik), vkb, use_gpu)
     !
#if defined(__OPENMP_GPU)
     !$omp target data map(to:vkb)
#endif
     CALL calbec( offload_type2, npw, vkb, wfcatom, becp )
#if defined(__OPENMP_GPU)
     !$omp end target data
     CALL s_psi_omp( npwx, npw, natomwfc, wfcatom, swfcatom )
#else
     CALL s_psi_acc( npwx, npw, natomwfc, wfcatom, swfcatom )
#endif
     !
     IF (orthogonalize_wfc) THEN
             CALL ortho_swfc ( npw, normalize_only, natomwfc, wfcatom, swfcatom, .FALSE. )
     ENDIF
     !
     ! write S * atomic wfc to unit iunsat
     !
     !$acc update host(swfcatom)
#if defined(__OPENMP_GPU)
     !$omp target update from(swfcatom)
#endif
     CALL save_buffer (swfcatom, nwordatwfc, iunsat, ik)
     !
  ENDDO
  !$acc exit data delete(wfcatom, swfcatom)
#if defined(__OPENMP_GPU)  
  !$omp end target data
#endif
  DEALLOCATE (swfcatom)
  DEALLOCATE ( wfcatom)
  CALL deallocate_bec_type_acc ( becp )
  !
  RETURN
     
END SUBROUTINE orthoatwfc
