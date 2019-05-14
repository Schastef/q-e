!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE cg_setupdgc
  !-----------------------------------------------------------------------
  ! Setup all arrays needed in the gradient correction case
  ! This version requires on input allocated array
  !
  USE kinds,     ONLY: dp
  USE constants, ONLY: e2
  USE scf,       ONLY: rho, rho_core, rhog_core, rhoz_or_updw
  USE funct,     ONLY: dft_is_gradient, init_gga_xc
  USE xc_gga,    ONLY: gcxc, gcx_spin, gcc_spin
  USE fft_base,  ONLY: dfftp
  USE gvect,     ONLY: ngm, g
  USE lsda_mod,  ONLY: nspin
  USE uspp,      ONLY: nlcc_any
  USE cgcom
  !
  IMPLICIT NONE
  INTEGER k, is
  real(DP) &
       &       grho2(2), rh(1), zeta(1), grh2(1),  &
       &       sx(1),sc(1),v1x(1,2),v2x(1,2),v1c(1,2),v2c(1,2), &
       &       vrrx(1,2),vsrx(1,2),vssx(1,2),                   &
       &       vrrc(1,2),vsrc(1,2),vssc(1),                     &
       &       vrzc(1,2), rho_v(1), grho2_v(1), s2_vvv(3,2)   !^^^
  !
  !
  !^^^
  real(dp), dimension(1,2) :: r_vv, s2_vv
  !^^^
  !
  REAL(DP), PARAMETER :: epsr=1.0d-6, epsg=1.0d-10
  !
  IF (.not. dft_is_gradient() ) RETURN
  CALL start_clock('setup_dgc')
  !
  dvxc_rr(:,:,:) = 0.d0
  dvxc_sr(:,:,:) = 0.d0
  dvxc_ss(:,:,:) = 0.d0
  dvxc_s (:,:,:) = 0.d0
  grho (:,:,:) = 0.d0
  !
  CALL init_gga_xc()
  !
  !    add rho_core to the charge density rho
  !
  IF (nlcc_any) THEN
     rho%of_r(:,1) = rho_core(:)  + rho%of_r(:,1)
     rho%of_g(:,1) = rhog_core(:) + rho%of_g(:,1)
  ENDIF
  !
  !    for LSDA, convert rho to (up,down) (gradient grho is (up,down))
  !
  IF (nspin == 2) CALL rhoz_or_updw( rho, 'r_and_g', '->updw' )
  DO is=1,nspin
     CALL fft_gradient_g2r (dfftp, rho%of_g(1,is), g, grho(1,1,is))
  ENDDO
  !
  IF (nspin==1) THEN
     DO k = 1,dfftp%nnr
        grho2(1)=grho(1,k,1)**2+grho(2,k,1)**2+grho(3,k,1)**2
        IF (abs(rho%of_r(k,1))>epsr.and.grho2(1)>epsg) THEN
           !
           !
           rho_v(1)=rho%of_r(k,nspin)
           grho2_v(1)=grho2(1)
           CALL gcxc( 1, rho_v,grho2_v(1),sx,sc,v1x(:,1),v2x(:,1),v1c(:,1),v2c(:,1) ) !^^^
           !
           !
           CALL dgcxc( 1, rho%of_r(k,nspin),grho2(1),vrrx(:,1),vsrx(:,1), &
                         vssx(:,1),vrrc(:,1),vsrc(:,1),vssc)
           dvxc_rr(k,1,1) = e2 * ( vrrx(1,1) + vrrc(1,1) )
           dvxc_sr(k,1,1) = e2 * ( vsrx(1,1) + vsrc(1,1) )
           dvxc_ss(k,1,1) = e2 * ( vssx(1,1) + vssc(1) )
           dvxc_s (k,1,1) = e2 * ( v2x(1,1) + v2c(1,1) )
        ENDIF
     ENDDO
  ELSE
     DO k = 1,dfftp%nnr
        grho2(1) = grho(1,k,1)**2 + grho(2,k,1)**2 + grho(3,k,1)**2
        grho2(2) = grho(1,k,2)**2 + grho(2,k,2)**2 + grho(3,k,2)**2
        rh(1)=rho%of_r(k,1) + rho%of_r(k,2)
        grh2(1)= (grho(1,k,1) + grho(1,k,2))**2              &
                           + (grho(2,k,1)+grho(2,k,2))**2 &
                           + (grho(3,k,1)+grho(3,k,2))**2
        !
        !CALL gcx_spin(rho%of_r(k,1),rho%of_r(k,2),grho2(1),grho2(2),sx, &
        !     v1xup,v1xdw,v2xup,v2xdw)
        !
        r_vv(1,1) = rho%of_r(k,1) ; r_vv(1,2) = rho%of_r(k,2)
        s2_vv(1,1) = grho2(1)   ; s2_vv(1,2) = grho2(2)
        call gcx_spin( 1, r_vv, s2_vv, sx, v1x, v2x )
        !
        !
        !
        s2_vvv(:,:) = grho(:,k,:)
        !
        CALL dgcxc_spin( 1, r_vv, s2_vvv, &
                         vrrx, vsrx, vssx, vrrc, vsrc, vssc, vrzc )
        !
        IF (rh(1)>epsr) THEN
           zeta=(rho%of_r(k,1)-rho%of_r(k,2))/rh(1)
           CALL gcc_spin( 1, rh, zeta, grh2, sc, v1c, v2c(:,1) )
           !
           dvxc_rr(k,1,1)=e2*(vrrx(1,1)+vrrc(1,1)+vrzc(1,1)*(1.d0-zeta(1))/rh(1))
           dvxc_rr(k,1,2)=e2*(vrrc(1,1)-vrzc(1,1)*(1.d0+zeta(1))/rh(1))
           dvxc_rr(k,2,1)=e2*(vrrc(1,2)+vrzc(1,2)*(1.d0-zeta(1))/rh(1))
           dvxc_rr(k,2,2)=e2*(vrrx(1,2)+vrrc(1,2)-vrzc(1,2)*(1.d0+zeta(1))/rh(1))
           !
           dvxc_s(k,1,1)=e2*(v2x(1,1)+v2c(1,1))
           dvxc_s(k,1,2)=e2*v2c(1,1)
           dvxc_s(k,2,1)=e2*v2c(1,1)
           dvxc_s(k,2,2)=e2*(v2x(1,2)+v2c(1,1))
        ELSE
           dvxc_rr(k,1,1)=0.d0
           dvxc_rr(k,1,2)=0.d0
           dvxc_rr(k,2,1)=0.d0
           dvxc_rr(k,2,2)=0.d0
           !
           dvxc_s(k,1,1)=0.d0
           dvxc_s(k,1,2)=0.d0
           dvxc_s(k,2,1)=0.d0
           dvxc_s(k,2,2)=0.d0
        ENDIF
        dvxc_sr(k,1,1)=e2*(vsrx(1,1)+vsrc(1,1))
        dvxc_sr(k,1,2)=e2*vsrc(1,1)
        dvxc_sr(k,2,1)=e2*vsrc(1,2)
        dvxc_sr(k,2,2)=e2*(vsrx(1,2)+vsrc(1,2))
        !
        dvxc_ss(k,1,1)=e2*(vssx(1,1)+vssc(1))
        dvxc_ss(k,1,2)=e2*vssc(1)
        dvxc_ss(k,2,1)=e2*vssc(1)
        dvxc_ss(k,2,2)=e2*(vssx(1,2)+vssc(1))
     ENDDO
  ENDIF
  !   restore rho to its input value
  IF (nspin == 2) CALL rhoz_or_updw( rho, 'r_and_g', '->rhoz' )
  IF (nlcc_any) THEN
     rho%of_r(:,1)  = rho%of_r(:,1) - rho_core(:)
     rho%of_g(:,1) = rho%of_g(:,1) - rhog_core(:)
  ENDIF
  CALL stop_clock('setup_dgc')
  !
  RETURN
END SUBROUTINE cg_setupdgc
