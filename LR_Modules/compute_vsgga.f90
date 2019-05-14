!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE compute_vsgga( rhoout, grho, vsgga )
  !----------------------------------------------------------------------------
  !
  USE constants,            ONLY : e2
  USE kinds,                ONLY : DP
  USE gvect,                ONLY : ngm, g
  USE noncollin_module,     ONLY : noncolin, nspin_gga
  USE funct,                ONLY : dft_is_gradient, get_igcc, init_gga_xc
  USE xc_gga,               ONLY : gcxc, gcx_spin, gcc_spin, gcc_spin_more
  USE spin_orb,             ONLY : domag
  USE fft_base,             ONLY : dfftp
  !
  IMPLICIT NONE
  !
  REAL(DP),    INTENT(IN)    :: rhoout(dfftp%nnr,nspin_gga)
  REAL(DP),    INTENT(IN)    :: grho(3,dfftp%nnr,nspin_gga)
  REAL(DP),    INTENT(OUT)   :: vsgga(dfftp%nnr)
  !
  INTEGER :: k, ipol, is
  !
  REAL(DP), ALLOCATABLE :: h(:,:,:), dh(:)
  REAL(DP), ALLOCATABLE :: vaux(:,:)
  !
  LOGICAL  :: igcc_is_lyp
  REAL(DP) :: grho2(2), sx(1), sc(1), v2c(1,2), &
              v1x(1,2), v2x(1,2), v1c(1,2), &
              arho, zeta(1), rh(1), grh2(1)
  REAL(DP) :: v2cud(1), r(1,2), &
              grhor(1,2), grhoud(1), gr(2)
  !
  !^^^TEMPORARY
  real(DP), dimension(1,2) :: r_vv, s2_vv
  !^^^
  !
  REAL(DP), PARAMETER :: vanishing_charge = 1.D-6, &
                         vanishing_mag    = 1.D-12
  REAL(DP), PARAMETER :: epsr = 1.D-6, epsg = 1.D-10
  !
  !
  IF ( .NOT. dft_is_gradient() ) RETURN
  
  CALL init_gga_xc()
  
  IF ( .NOT. (noncolin.and.domag) ) &
     call errore('compute_vsgga','routine called in the wrong case',1)

  igcc_is_lyp = (get_igcc() == 3)
  !
  ALLOCATE( h(3,dfftp%nnr,nspin_gga)  )
  ALLOCATE( vaux(dfftp%nnr,nspin_gga) )

  DO k = 1, dfftp%nnr
     !
     rh(1) = rhoout(k,1) + rhoout(k,2)
     !
     arho=abs(rh(1))
     !
     IF ( arho > vanishing_charge ) THEN
        !
        grho2(:) = grho(1,k,:)**2 + grho(2,k,:)**2 + grho(3,k,:)**2
        !
        IF ( grho2(1) > epsg .OR. grho2(2) > epsg ) THEN
           !CALL gcx_spin( rhoout(k,1), rhoout(k,2), grho2(1), &
           !               grho2(2), sx, v1xup, v1xdw, v2xup, v2xdw )
           r_vv(1,1)=rhoout(k,1) ; r_vv(1,2)=rhoout(k,2)
           s2_vv(1,1)=grho2(1)   ; s2_vv(1,2)=grho2(2)
           call gcx_spin( 1, r_vv, s2_vv, sx, v1x, v2x )
           !
           IF ( igcc_is_lyp ) THEN
              !
              r(1,1) = rhoout(k,1)
              r(1,2) = rhoout(k,2)
              !
              grhor(1,1) = grho(1,k,1)**2 + grho(2,k,1)**2 + grho(3,k,1)**2
              grhor(1,2) = grho(1,k,2)**2 + grho(2,k,2)**2 + grho(3,k,2)**2
              !
              grhoud = grho(1,k,1) * grho(1,k,2) + &
                       grho(2,k,1) * grho(2,k,2) + &
                       grho(3,k,1) * grho(3,k,2)
              !
              CALL gcc_spin_more( 1, r, grhor, grhoud, sc, v1c, v2c, v2cud )
              !
           ELSE
              !
              zeta = ( rhoout(k,1) - rhoout(k,2) ) / rh
              !
              grh2 = ( grho(1,k,1) + grho(1,k,2) )**2 + &
                     ( grho(2,k,1) + grho(2,k,2) )**2 + &
                     ( grho(3,k,1) + grho(3,k,2) )**2
              !
              CALL gcc_spin( 1, rh, zeta, grh2, sc, v1c, v2c(:,1) )
              !
              v2c(:,2) = v2c(:,1)
              v2cud = v2c(:,1)
              !
           END IF
        ELSE
           !
           sc  = 0.D0
           sx  = 0.D0
           v1x = 0.D0
           v2x = 0.D0
           v1c = 0.D0
           v2c = 0.D0
           v2cud = 0.D0
        ENDIF
     ELSE
        !
        sc  = 0.D0
        sx  = 0.D0
        v1x = 0.D0
        v2x = 0.D0
        v1c = 0.D0
        v2c = 0.D0
        v2c = 0.D0
        v2cud = 0.D0
        !
     ENDIF
     !
     ! ... first term of the gradient correction : D(rho*Exc)/D(rho)
     !
     vaux(k,1) = e2 * ( v1x(1,1) + v1c(1,1) )
     vaux(k,2) = e2 * ( v1x(1,2) + v1c(1,2) )
     !
     ! ... h contains D(rho*Exc)/D(|grad rho|) * (grad rho) / |grad rho|
     !
     DO ipol = 1, 3
        !
        gr(:) = grho(ipol,k,:)
        h(ipol,k,1) = e2 * ( ( v2x(1,1) + v2c(1,1) ) * gr(1) + v2cud(1) * gr(2) )
        h(ipol,k,2) = e2 * ( ( v2x(1,2) + v2c(1,2) ) * gr(2) + v2cud(1) * gr(1) )
        !
     END DO
     !
  END DO
  !
  ALLOCATE( dh( dfftp%nnr ) )
  !
  ! ... second term of the gradient correction :
  ! ... \sum_alpha (D / D r_alpha) ( D(rho*Exc)/D(grad_alpha rho) )
  !
  DO is = 1, nspin_gga
     !
     CALL fft_graddot( dfftp, h(1,1,is), g, dh )
     !
     vaux(:,is) = vaux(:,is) - dh(:)
     !
  END DO

  vsgga(:)=(vaux(:,1)-vaux(:,2))

  !
  DEALLOCATE( dh )
  DEALLOCATE( h )
  DEALLOCATE( vaux )
  !
  RETURN
  !
END SUBROUTINE compute_vsgga
!
