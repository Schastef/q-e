!-----------------------------------------------------------------------
!------- DRIVERS FOR DERIVATIVES OF XC POTENTIAL (GGA CASE) ------------
!-----------------------------------------------------------------------
!
!---------------------------------------------------------------------------
SUBROUTINE dgcxc( length, r_in, s2_in, vrrx, vsrx, vssx, vrrc, vsrc, vssc )
  !-------------------------------------------------------------------------
  !
  USE xc_gga,       ONLY: gcxc, igcx_l, igcc_l
  USE funct,        ONLY: init_gga_xc  
  USE kinds,        ONLY: DP
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: length
  REAL(DP), INTENT(IN), DIMENSION(length) :: r_in
  REAL(DP), INTENT(IN), DIMENSION(length) :: s2_in
  REAL(DP), INTENT(OUT), DIMENSION(length) :: vrrx, vsrx, vssx
  REAL(DP), INTENT(OUT), DIMENSION(length) :: vrrc, vsrc, vssc
  !
  ! ... local variables
  !
  INTEGER :: ir
  REAL(DP), DIMENSION(1) :: r, s2, dr, s, ds
  REAL(DP), DIMENSION(1) :: v1xp, v2xp, v1cp, v2cp, &
                            v1xm, v2xm, v1cm, v2cm
  REAL(DP), DIMENSION(1) :: sx, sc
  !
  !
  CALL init_gga_xc()
  !
  DO ir = 1, length
     !
     r  = r_in(ir)
     s2 = s2_in(ir)
     s  = SQRT(s2)
     dr = MIN( 1.d-4, 1.d-2*r )
     ds = MIN( 1.d-4, 1.d-2*s )
     !
     CALL gcxc( 1, r+dr, s2, sx, sc, v1xp, v2xp, v1cp, v2cp )
     !
     CALL gcxc( 1, r-dr, s2, sx, sc, v1xm, v2xm, v1cm, v2cm )
     !
     vrrx(ir) = 0.5_DP  * (v1xp(1) - v1xm(1)) / dr(1)
     vrrc(ir) = 0.5_DP  * (v1cp(1) - v1cm(1)) / dr(1)
     vsrx(ir) = 0.25_DP * (v2xp(1) - v2xm(1)) / dr(1)
     vsrc(ir) = 0.25_DP * (v2cp(1) - v2cm(1)) / dr(1)
     !
     CALL gcxc( 1, r, (s+ds)**2, sx, sc, v1xp, v2xp, v1cp, v2cp )
     !
     CALL gcxc( 1, r, (s-ds)**2, sx, sc, v1xm, v2xm, v1cm, v2cm )
     !
     vsrx(ir) = vsrx(ir) + 0.25_DP * (v1xp(1) - v1xm(1)) / ds(1) / s(1)
     vsrc(ir) = vsrc(ir) + 0.25_DP * (v1cp(1) - v1cm(1)) / ds(1) / s(1)
     vssx(ir) = 0.5_DP * (v2xp(1) - v2xm(1)) / ds(1) / s(1)
     vssc(ir) = 0.5_DP * (v2cp(1) - v2cm(1)) / ds(1) / s(1)
     !
  ENDDO
  !
  RETURN
  !
END SUBROUTINE dgcxc
!
!--------------------------------------------------------------------------
SUBROUTINE dgcxc_spin( length, r_in, g_in, vrrx, vrsx, vssx, vrrc, vrsc, &
                       vssc, vrzc )
  !------------------------------------------------------------------------
  !! This routine computes the derivative of the exchange and correlation
  !! potentials with respect to the density, the gradient and zeta
  !
  USE xc_gga,       ONLY: gcx_spin, gcc_spin, igcx_l, igcc_l
  USE funct,        ONLY: init_gga_xc
  USE kinds,        ONLY: DP
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: length
  REAL(DP), INTENT(IN), DIMENSION(length,2) :: r_in
  REAL(DP), INTENT(IN), DIMENSION(length,3,2) :: g_in
  ! input: the charges and the gradient
  REAL(DP), INTENT(OUT), DIMENSION(length,2) :: vrrx, vrsx, vssx
  REAL(DP), INTENT(OUT), DIMENSION(length,2) :: vrrc, vrsc, vrzc
  REAL(DP), INTENT(OUT), DIMENSION(length) :: vssc
  ! output: derivatives of the exchange and of the correlation
  !
  ! ... local variables
  !
  INTEGER :: ir
  REAL(DP), DIMENSION(1,2) :: r, s, s2
  REAL(DP), DIMENSION(1,2) :: drup, drdw, dsup, dsdw
  REAL(DP), DIMENSION(1,2) :: v1xp, v2xp, v1cp, v1xm, v2xm, v1cm
  !
  REAL(DP), DIMENSION(1) :: rt, zeta, zetapm, st, s2t
  REAL(DP), DIMENSION(1) :: dr, ds, dzeta
  REAL(DP), DIMENSION(1) :: sx, sc, v2cp, v2cm
  !
  ! charge densities and square gradients
  ! delta charge densities and gra
  ! delta gradients
  ! energies
  ! exchange potentials
  ! exchange potentials
  ! coorelation potentials
  ! coorelation potentials
  !
  REAL(DP), PARAMETER :: eps = 1.D-6
  !
  !write(*,*) g_in(1,1:3,1)
  !
  CALL init_gga_xc()
  !
  DO ir = 1, length
     !
     drup = 0.0_DP
     drdw = 0.0_DP
     dsup = 0.0_DP
     dsdw = 0.0_DP
     !
     r(1,:) = r_in(ir,:)
     rt(1) = r(1,1) + r(1,2)
     !
     IF ( rt(1) > eps ) THEN
        zeta(1) = (r(1,1) - r(1,2)) / rt(1)
     ELSE
        zeta(1) = 2.0_DP
     ENDIF
     !
     s2(1,1) = g_in(ir,1,1)**2 + g_in(ir,2,1)**2 + g_in(ir,3,1)**2 !up
     s2(1,2) = g_in(ir,1,2)**2 + g_in(ir,2,2)**2 + g_in(ir,3,2)**2 !down
     s2t(1)  = (g_in(ir,1,1) + g_in(ir,1,2))**2 + &
               (g_in(ir,2,1) + g_in(ir,2,2))**2 + &
               (g_in(ir,3,1) + g_in(ir,3,2))**2
     s(1,1) = SQRT(s2(1,1))
     s(1,2) = SQRT(s2(1,2))
     st(1)  = SQRT(s2t(1))
     !
     ! up part of exchange
     !
     IF ( r(1,1)>eps .AND. s(1,1)>eps ) THEN
        !
        drup(1,1) = MIN(1.D-4, 1.D-2 * r(1,1))
        dsup(1,1) = MIN(1.D-4, 1.D-2 * s(1,1))
        !
        ! ... derivatives of exchange: up part
        !
        CALL gcx_spin( 1, r+drup, s2, sx, v1xp, v2xp )
        !
        CALL gcx_spin( 1, r-drup, s2, sx, v1xm, v2xm )
        !
        vrrx(ir,1) = 0.5_DP  * (v1xp(1,1) - v1xm(1,1)) / drup(1,1)
        vrsx(ir,1) = 0.25_DP * (v2xp(1,1) - v2xm(1,1)) / drup(1,1)
        !
        CALL gcx_spin( 1, r, (s+dsup)**2, sx, v1xp, v2xp )
        !
        CALL gcx_spin( 1, r, (s-dsup)**2, sx, v1xm, v2xm )
        !
        vrsx(ir,1) = vrsx(ir,1) + 0.25_DP * (v1xp(1,1) - v1xm(1,1)) / dsup(1,1) / s(1,1)
        vssx(ir,1) = 0.5_DP * (v2xp(1,1) - v2xm(1,1)) / dsup(1,1) / s(1,1)
        !
     ELSE
        !
        vrrx(ir,1) = 0.0_DP
        vrsx(ir,1) = 0.0_DP
        vssx(ir,1) = 0.0_DP
        !
     ENDIF
     !
     !
     IF ( r(1,2)>eps .AND. s(1,2)>eps ) THEN
        !
        drdw(1,2) = MIN(1.D-4, 1.D-2 * r(1,2))
        dsdw(1,2) = MIN(1.D-4, 1.D-2 * s(1,2))
        !
        ! ... derivatives of exchange: up part
        !
        CALL gcx_spin( 1, r+drdw, s2, sx, v1xp, v2xp )
        !
        CALL gcx_spin( 1, r-drdw, s2, sx, v1xm, v2xm )
        !
        vrrx(ir,2) = 0.5_DP  * (v1xp(1,2) - v1xm(1,2)) / drdw(1,2)
        vrsx(ir,2) = 0.25_DP * (v2xp(1,2) - v2xm(1,2)) / drdw(1,2)
        !
        CALL gcx_spin( 1, r, (s+dsdw)**2, sx, v1xp, v2xp )
        !
        CALL gcx_spin( 1, r, (s-dsdw)**2, sx, v1xm, v2xm )
        !
        vrsx(ir,2) = vrsx(ir,2) + 0.25_DP*(v1xp(1,2) - v1xm(1,2)) / &
                                                  dsdw(1,2) / s(1,2)
        vssx(ir,2) = 0.5_DP*(v2xp(1,2) - v2xm(1,2)) / dsdw(1,2) / s(1,2)
        !
     ELSE
        !
        vrrx(ir,2) = 0.0_DP
        vrsx(ir,2) = 0.0_DP
        vssx(ir,2) = 0.0_DP
        !
     ENDIF
     !
     ! ... derivatives of correlation
     !
     IF ( rt(1)>eps .AND. ABS(zeta(1))<=1._DP .AND. st(1)>eps ) THEN
        !
        dr = MIN(1.D-4, 1.D-2 * rt)
        !
        CALL gcc_spin( 1, rt+dr, zeta, s2t, sc, v1cp, v2cp )
        !
        CALL gcc_spin( 1, rt-dr, zeta, s2t, sc, v1cm, v2cm )
        !
        vrrc(ir,1) = 0.5_DP * (v1cp(1,1) - v1cm(1,1)) / dr(1)
        vrrc(ir,2) = 0.5_DP * (v1cp(1,2) - v1cm(1,2)) / dr(1)
        !
        ds = MIN(1.D-4, 1.D-2 * st)
        !
        CALL gcc_spin( 1, rt, zeta, (st+ds)**2, sc, v1cp, v2cp )
        !
        CALL gcc_spin( 1, rt, zeta, (st-ds)**2, sc, v1cm, v2cm )
        !
        vrsc(ir,1) = 0.5_DP * (v1cp(1,1) - v1cm(1,1)) / ds(1) / st(1)
        vrsc(ir,2) = 0.5_DP * (v1cp(1,2) - v1cm(1,2)) / ds(1) / st(1)
        vssc(ir)   = 0.5_DP * (v2cp(1)   - v2cm(1)  ) / ds(1) / st(1)
        !
        !dzeta(1) = MIN(1.D-4, 1.D-2 * ABS(zeta(1)) )
        dzeta(1) = 1.D-6
        !
        ! ... If zeta is too close to +-1 the derivative is evaluated at a
        ! slightly smaller value.
        !
        zeta = SIGN( MIN( ABS( zeta ), ( 1.0_DP - 2.0_DP*dzeta ) ) , zeta )
        !
        zetapm = zeta+dzeta
        CALL gcc_spin( 1, rt, zetapm, s2t, sc, v1cp, v2cp )
        !
        zetapm = zeta-dzeta
        CALL gcc_spin( 1, rt, zetapm, s2t, sc, v1cm, v2cm )
        !
        vrzc(ir,1) = 0.5_DP * (v1cp(1,1) - v1cm(1,1)) / dzeta(1)
        vrzc(ir,2) = 0.5_DP * (v1cp(1,2) - v1cm(1,2)) / dzeta(1)
        !
     ELSE
        !
        vrrc = 0.0_DP
        vrsc = 0.0_DP
        vssc = 0.0_DP
        vrzc = 0.0_DP
        !
     ENDIF
     !
  ENDDO
  !
  !
  RETURN
  !
END SUBROUTINE dgcxc_spin
!
!
!-----------------------------------------------------------------------
SUBROUTINE d3gcxc( r, s2, vrrrx, vsrrx, vssrx, vsssx, &
                   vrrrc, vsrrc, vssrc, vsssc )
  !-----------------------------------------------------------------------
  !
  !    wat20101006: Calculates all derivatives of the exchange (x) and
  !                 correlation (c) potential in third order.
  !                 of the Exc.
  !
  !    input:       r = rho, s2=|\nabla rho|^2
  !    definition:  E_xc = \int ( f_x(r,s2) + f_c(r,s2) ) dr
  !    output:      vrrrx = d^3(f_x)/d(r)^3
  !                 vsrrx = d^3(f_x)/d(|\nabla r|)d(r)^2 / |\nabla r|
  !                 vssrx = d/d(|\nabla r|) [ &
  !                           d^2(f_x)/d(|\nabla r|)d(r) / |\nabla r| ] &
  !                                                           / |\nabla r|
  !                 vsssx = d/d(|\nabla r|) [ &
  !                           d/d(|\nabla r|) [ &
  !                           d(f_x)/d(|\nabla r|) / |\nabla r| ] &
  !                                                   / |\nabla r| ] &
  !                                                   / |\nabla r|
  !                 same for (c)
  !
  USE kinds, ONLY : DP
  IMPLICIT NONE
  REAL(DP) :: r, s2, vrrrx, vsrrx, vssrx, vsssx, &
              vrrrc, vsrrc, vssrc, vsssc
  REAL(DP) :: dr, s, ds
  !
  REAL(DP) :: vrrx_rp, vsrx_rp, vssx_rp, vrrc_rp, vsrc_rp, vssc_rp, &
              vrrx_rm, vsrx_rm, vssx_rm, vrrc_rm, vsrc_rm, vssc_rm, &
              vrrx_sp, vsrx_sp, vssx_sp, vrrc_sp, vsrc_sp, vssc_sp, &
              vrrx_sm, vsrx_sm, vssx_sm, vrrc_sm, vsrc_sm, vssc_sm
  !
  s = sqrt (s2)
  dr = min (1.d-4, 1.d-2 * r)
  ds = min (1.d-4, 1.d-2 * s)
  !
  !call dgcxc (r+dr, s2, vrrx_rp, vsrx_rp, vssx_rp, vrrc_rp, vsrc_rp, vssc_rp)  !^^^ RIMETTILE TUTTE E QUATTRO
  !call dgcxc (r-dr, s2, vrrx_rm, vsrx_rm, vssx_rm, vrrc_rm, vsrc_rm, vssc_rm)
  !
  !call dgcxc (r, (s+ds)**2, vrrx_sp, vsrx_sp, vssx_sp, vrrc_sp, vsrc_sp, vssc_sp)
  !call dgcxc (r, (s-ds)**2, vrrx_sm, vsrx_sm, vssx_sm, vrrc_sm, vsrc_sm, vssc_sm)
  !
  vrrrx = 0.5d0 * (vrrx_rp - vrrx_rm) / dr
  vsrrx = 0.25d0 * (vsrx_rp - vsrx_rm) / dr &
                + 0.25d0 * (vrrx_sp - vrrx_sm) / ds / s
  vssrx = 0.25d0 * (vssx_rp - vssx_rm) / dr &
                + 0.25d0 * (vsrx_sp - vsrx_sm) / ds / s
  vsssx = 0.5d0 * (vssx_sp - vssx_sm) / ds / s
  !
  vrrrc = 0.5d0 * (vrrc_rp - vrrc_rm) / dr
  vsrrc = 0.25d0 * (vsrc_rp - vsrc_rm) / dr &
                + 0.25d0 * (vrrc_sp - vrrc_sm) / ds / s
  vssrc = 0.25d0 * (vssc_rp - vssc_rm) / dr &
                + 0.25d0 * (vsrc_sp - vsrc_sm) / ds / s
  vsssc = 0.5d0 * (vssc_sp - vssc_sm) / ds / s
  !
  return
  !
END SUBROUTINE d3gcxc
