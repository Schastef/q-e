!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine setup_dgc
  !-----------------------------------------------------------------------
  !
  ! Allocate and setup all variable needed in the gradient correction case
  !
  ! GGA+LSDA is allowed. ADC (September 1999).
  ! GGA+LSDA+NLCC is allowed. ADC (November 1999).
  ! GGA+noncollinear+NLCC is allowed. ADC (June 2007).
  !
  USE constants,            ONLY : e2
  USE fft_base,             ONLY : dfftp
  USE fft_interfaces,       ONLY : fwfft
  USE gvect,                ONLY : ngm, g
  USE spin_orb,             ONLY : domag
  USE scf,                  ONLY : rho, rho_core, rhog_core, rhoz_or_updw
  USE noncollin_module,     ONLY : noncolin, ux, nspin_gga, nspin_mag
  USE wavefunctions,        ONLY : psic
  USE kinds,                ONLY : DP
  USE funct,                ONLY : dft_is_gradient, init_gga_xc
  USE xc_gga,               ONLY : gcxc, gcx_spin, gcc_spin
  USE uspp,                 ONLY : nlcc_any
  USE gc_lr,                ONLY : grho, gmag, dvxc_rr, dvxc_sr, &
                                   dvxc_ss, dvxc_s, vsgga, segni

  implicit none
  integer :: k, is, ipol, jpol, ir
  real(DP) :: grho2(2), rh(1), zeta(1), grh2(1), fac, rho_v(1), grho2_v(1,2), &  !^^^
              amag, seg, seg0, sgn(2)
  !
  real(dp), dimension(1,nspin_gga) :: v1x, v2x, v1c, v2c
  real(dp), dimension(1,nspin_gga) :: vrrx, vsrx, vssx, vrrc, vsrc, vrzc
  real(dp), dimension(1) :: vssc, sc, sx
  !
  !^^^TEMPORARY
  real(dp), dimension(1,2) :: r_vv, s2_vv
  real(dp), dimension(1,3,2) :: s2_vvv
  !^^^
  !
  COMPLEX(DP), ALLOCATABLE :: rhogout(:,:)
  real(DP), allocatable :: rhoout(:,:)
  real (DP), parameter :: epsr = 1.0d-6, epsg = 1.0d-10

  IF ( .NOT. dft_is_gradient() ) RETURN
 
  CALL start_clock ('setup_dgc')
  
  CALL init_gga_xc()
  
  IF (noncolin.AND.domag) THEN
     allocate (segni (dfftp%nnr))
     allocate (vsgga (dfftp%nnr))
     allocate (gmag (3, dfftp%nnr, nspin_mag))
     gmag=0.0_dp
  ENDIF

  IF(.NOT.ALLOCATED(dvxc_rr)) ALLOCATE (dvxc_rr(dfftp%nnr, nspin_gga , nspin_gga))
  IF(.NOT.ALLOCATED(dvxc_sr)) ALLOCATE (dvxc_sr(dfftp%nnr, nspin_gga , nspin_gga))
  IF(.NOT.ALLOCATED(dvxc_ss)) ALLOCATE (dvxc_ss(dfftp%nnr, nspin_gga , nspin_gga))
  IF(.NOT.ALLOCATED(dvxc_s))  ALLOCATE (dvxc_s (dfftp%nnr, nspin_gga , nspin_gga))
  IF(.NOT.ALLOCATED(grho))    ALLOCATE (grho   (  3    , dfftp%nnr, nspin_gga))
  IF(.NOT.ALLOCATED(rhoout))  ALLOCATE (rhoout ( dfftp%nnr, nspin_gga))

  dvxc_rr(:,:,:) = 0.d0
  dvxc_sr(:,:,:) = 0.d0
  dvxc_ss(:,:,:) = 0.d0
  dvxc_s (:,:,:) = 0.d0
  grho   (:,:,:) = 0.d0
  sgn(1)=1.d0  ;   sgn(2)=-1.d0
  !
  !    add rho_core
  !
  fac = 1.d0 / DBLE (nspin_gga)
  IF (noncolin.and.domag) THEN
     allocate(rhogout(ngm,nspin_mag))
     call compute_rho(rho%of_r,rhoout,segni,dfftp%nnr)
     DO is = 1, nspin_gga
        !
        if (nlcc_any) rhoout(:,is)  = fac * rho_core(:)  + rhoout(:,is)

        psic(:) = rhoout(:,is)
        !
        CALL fwfft ('Rho', psic, dfftp)
        !
        rhogout(:,is) = psic(dfftp%nl(:))
        !
        CALL fft_gradient_g2r(dfftp, rhogout(1,is), g, grho(1,1,is) )
        !
     END DO
     DEALLOCATE(rhogout)
  ELSE
     !
     ! for convenience, if LSDA, rhoout is kept in (up,down) format
     !
     do is = 1, nspin_gga
        rhoout(:,is) = ( rho%of_r(:,1) + sgn(is)*rho%of_r(:,nspin_gga) )*0.5d0
     enddo
     !
     ! if LSDA rho%of_g is temporarily converted in (up,down) format
     !
     call rhoz_or_updw(rho, 'only_g', '->updw')
     !
     if (nlcc_any) then
        do is = 1, nspin_gga
           rhoout(:,is)   = fac * rho_core(:)  + rhoout(:,is)
           rho%of_g(:,is) = fac * rhog_core(:) + rho%of_g(:,is)
        enddo
     endif
     !
     do is = 1, nspin_gga
        call fft_gradient_g2r (dfftp, rho%of_g (1, is), g, grho (1, 1, is) )
     enddo
  END IF

  do k = 1, dfftp%nnr
     grho2 (1) = grho (1, k, 1) **2 + grho (2, k, 1) **2 + grho (3, k, 1) **2
     if (nspin_gga == 1) then
        if (abs (rhoout(k, 1) ) > epsr .and. grho2(1) > epsg) then
           !
           rho_v(1)=rhoout(k,1)
           grho2_v(:,1)=grho2(1)
           !
           call gcxc( 1, rho_v, grho2_v(:,1), sx, sc, v1x(:,1), &
                                             v2x(:,1), v1c(:,1), v2c(:,1) )
           !
           !WRITE(*,*) 'ttest gxcx: setup_dgc'
           !
           call dgcxc( 1, rho_v, grho2_v, vrrx(:,1), vsrx(:,1), vssx(:,1), &
                                        vrrc(:,1), vsrc(:,1), vssc ) !^^^
           dvxc_rr(k,1,1) = e2 * (vrrx(1,1) + vrrc(1,1))
           dvxc_sr(k,1,1) = e2 * (vsrx(1,1) + vsrc(1,1))
           dvxc_ss(k,1,1) = e2 * (vssx(1,1) + vssc(1))
           dvxc_s(k,1,1)  = e2 * (v2x(1,1)  + v2c(1,1))
        endif
     else
        grho2 (2) = grho(1, k, 2) **2 + grho(2, k, 2) **2 + grho(3, k, 2) **2
        rh = rhoout (k, 1) + rhoout (k, 2)

        grh2 = (grho (1, k, 1) + grho (1, k, 2) ) **2 + (grho (2, k, 1) &
             + grho (2, k, 2) ) **2 + (grho (3, k, 1) + grho (3, k, 2) ) ** 2

        !call gcx_spin (rhoout (k, 1), rhoout (k, 2), grho2 (1), grho2 (2), & !^^^
        !     sx, v1xup, v1xdw, v2xup, v2xdw)
        r_vv(1,1)=rhoout(k,1) ; r_vv(1,2)=rhoout(k,2)
        s2_vv(1,1)=grho2(1)   ; s2_vv(1,2)=grho2(2)
        !
        call gcx_spin( 1, r_vv, s2_vv, sx, v1x, v2x )
        !
        !
        s2_vvv(1,1:3,1)=grho(1:3,k,1)
        s2_vvv(1,1:3,2)=grho(1:3,k,2)
        call dgcxc_spin ( 1, r_vv, s2_vvv, vrrx, vsrx, vssx, vrrc, vsrc, vssc, vrzc )
        
        if (rh(k) > epsr) then
           zeta = (rhoout(k,1) - rhoout(k,2) ) / rh
           !
           call gcc_spin( 1, rh, zeta, grh2, sc, v1c, v2c(:,1) )  !^^^
           !
           dvxc_rr(k,1,1) = e2 * (vrrx(1,1) + vrrc(1,1) + vrzc(1,1) * (1.d0 - zeta(1)) / rh(1))
           dvxc_rr(k,1,2) = e2 * (vrrc(1,1) - vrzc(1,1) * (1.d0 + zeta(1)) / rh(1))
           dvxc_rr(k,2,1) = e2 * (vrrc(1,2) + vrzc(1,2) * (1.d0 - zeta(1)) / rh(1))
           dvxc_rr(k,2,2) = e2 * (vrrx(1,2) + vrrc(1,2) - vrzc(1,2) * (1.d0 + zeta(1)) / rh(1))
           dvxc_s(k,1,1) = e2 * (v2x(1,1) + v2c(1,1))
           dvxc_s(k,1,2) = e2 * v2c(1,1)
           dvxc_s(k,2,1) = e2 * v2c(1,1)
           dvxc_s(k,2,2) = e2 * (v2x(1,2) + v2c(1,1))
        else
           dvxc_rr(k,1,1) = 0.d0
           dvxc_rr(k,1,2) = 0.d0
           dvxc_rr(k,2,1) = 0.d0
           dvxc_rr(k,2,2) = 0.d0
           dvxc_s(k,1,1) = 0.d0
           dvxc_s(k,1,2) = 0.d0
           dvxc_s(k,2,1) = 0.d0
           dvxc_s(k,2,2) = 0.d0
        endif
        dvxc_sr(k,1,1) = e2 * (vsrx(1,1) + vsrc(1,1))
        dvxc_sr(k,1,2) = e2 * vsrc(1,1)
        dvxc_sr(k,2,1) = e2 * vsrc(1,2)
        dvxc_sr(k,2,2) = e2 * (vsrx(1,2) + vsrc(1,2))
        dvxc_ss(k,1,1) = e2 * (vssx(1,1) + vssc(1))
        dvxc_ss(k,1,2) = e2 * vssc(1)
        dvxc_ss(k,2,1) = e2 * vssc(1)
        dvxc_ss(k,2,2) = e2 * (vssx(1,2) + vssc(1))
        
     endif
  enddo
  !
  if (noncolin.and.domag) then
     call compute_vsgga(rhoout, grho, vsgga)
  else
     if (nlcc_any) then
        do is = 1, nspin_gga
           rho%of_g(:,is) = rho%of_g(:,is) - fac * rhog_core(:)
        enddo
     endif
     !
     CALL rhoz_or_updw(rho, 'only_g', '->rhoz')
     !
  endif

  DEALLOCATE(rhoout)

  CALL stop_clock ('setup_dgc')

  RETURN

end subroutine setup_dgc
