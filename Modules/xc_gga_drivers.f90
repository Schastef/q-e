MODULE xc_gga
!
USE kinds,     ONLY: DP 
!
IMPLICIT NONE
!
PRIVATE
SAVE
!
!  GGA exchange-correlation drivers
PUBLIC :: gcxc, gcx_spin, gcc_spin, gcc_spin_more, &
          select_gga_functionals
!
PUBLIC :: igcx_l, igcc_l
PUBLIC :: exx_started_g, exx_fraction_g
PUBLIC :: screening_parameter_l, gau_parameter_l
!
!  indexes defining xc functionals
INTEGER  :: igcx_l, igcc_l
!
!  variables for hybrid exchange
LOGICAL  :: exx_started_g
REAL(DP) :: exx_fraction_g
!
!  screening_ and gau_parameters
REAL(DP) :: screening_parameter_l, gau_parameter_l
!
!
 CONTAINS
!
!
!----------------------------------------------------------------------------
!----- Select functionals by the corresponding indexes ----------------------
!----------------------------------------------------------------------------
SUBROUTINE select_gga_functionals( igcx, igcc, exx_fraction, screening_parameter, &
                                    gau_parameter )
   !
   IMPLICIT NONE
   !
   INTEGER,  INTENT(IN) :: igcx, igcc
   REAL(DP), INTENT(IN), OPTIONAL :: exx_fraction
   REAL(DP), INTENT(IN), OPTIONAL :: screening_parameter
   REAL(DP), INTENT(IN), OPTIONAL :: gau_parameter
   !
   ! exchange-correlation indexes
   igcx_l = igcx
   igcc_l = igcc
   !
   ! hybrid exchange vars
   exx_started_g  = .FALSE.
   exx_fraction_g = 0._DP
   IF ( PRESENT(exx_fraction) ) THEN
      exx_started_g  = .TRUE.
      exx_fraction_g = exx_fraction
   ENDIF
   !
   ! screening_ and gau_parameter
   screening_parameter_l = 0.0_DP
   gau_parameter_l = 0.0_DP
   !
   IF ( PRESENT(screening_parameter) ) THEN          !^^^ SISTEMA, metti compatiblita' con indici e 
      screening_parameter_l = screening_parameter         ! vedi ordine variabili optional
   ENDIF
   !
   IF ( PRESENT(gau_parameter) ) THEN
      gau_parameter_l = gau_parameter
   ENDIF
   !
   RETURN
   !
END SUBROUTINE select_gga_functionals
!
!
!-----------------------------------------------------------------------
!------- GRADIENT CORRECTIONS DRIVERS ----------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE gcxc( length, rho_vi, grho_vi, sx, sc, v1x, v2x, v1c, v2c )
  !---------------------------------------------------------------------
  !! Gradient corrections for exchange and correlation - Hartree a.u. 
  !! See comments at the beginning of module for implemented cases
  !
  ! Input:  rho, grho=|\nabla rho|^2
  ! Definition:  E_x = \int E_x(rho,grho) dr
  ! Output: sx = E_x(rho,grho)
  !         v1x= D(E_x)/D(rho)
  !         v2x= D(E_x)/D( D rho/D r_alpha ) / |\nabla rho|
  !         sc, v1c, v2c as above for correlation
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: length
  REAL(DP), INTENT(IN),  DIMENSION(length) :: rho_vi, grho_vi
  REAL(DP), INTENT(OUT), DIMENSION(length) :: sx, sc, v1x, v2x, v1c, v2c
  !
  ! ... local variables
  !
  REAL(DP), DIMENSION(length) :: rho, grho
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: sx_, v1x_, v2x_
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: sxsr, v1xsr, v2xsr
  REAL(DP), PARAMETER :: small = 1.E-10_DP
  !
  !
  rho  = rho_vi
  grho = grho_vi
  WHERE ( rho <= small )
     rho  = 0.5_DP
     grho = 0.2_DP
  END WHERE
  !
  IF ( igcx_l == 28 ) THEN
     !
     ALLOCATE( sx_(length)                )
     ALLOCATE( v1x_(length), v2x_(length) )
     !
  ELSEIF ( (igcx_l == 12 .OR. igcx_l == 20) .AND. exx_started_g ) THEN
     !
     ALLOCATE( sxsr(length)                 )
     ALLOCATE( v1xsr(length), v2xsr(length) )
     !
  ENDIF
  !
  !
  ! ... EXCHANGE
  !  
  SELECT CASE( igcx_l )
  CASE( 1 )
     !
     CALL becke88( length, rho, grho, sx, v1x, v2x )
     !
  CASE( 2 )
     !
     CALL ggax( length, rho, grho, sx, v1x, v2x )
     !
  CASE( 3 )
     !
     CALL pbex( length, rho, grho, 1, sx, v1x, v2x )
     !
  CASE( 4 )
     !
     CALL pbex( length, rho, grho, 2, sx, v1x, v2x )
     !
  CASE( 5 )
     !
     IF (igcc_l == 5) CALL hcth( length, rho, grho, sx, v1x, v2x )
     !
  CASE( 6 )
     !
     CALL optx( length, rho, grho, sx, v1x, v2x )
     !
  ! case igcx_l == 7 (meta-GGA) must be treated in a separate call to another
  ! routine: needs kinetic energy density in addition to rho and grad rho
  CASE( 8 ) ! 'PBE0'
     !
     CALL pbex( length, rho, grho, 1, sx, v1x, v2x )
     IF (exx_started_g) THEN
        sx  = (1.0_DP - exx_fraction_g) * sx
        v1x = (1.0_DP - exx_fraction_g) * v1x
        v2x = (1.0_DP - exx_fraction_g) * v2x
     ENDIF
     !
  CASE( 9 ) ! 'B3LYP'
     !
     CALL becke88( length, rho, grho, sx, v1x, v2x )
     IF (exx_started_g) THEN
        sx  = 0.72_DP * sx
        v1x = 0.72_DP * v1x
        v2x = 0.72_DP * v2x
     ENDIF
     !
  CASE( 10 ) ! 'pbesol'
     !
     CALL pbex( length, rho, grho, 3, sx, v1x, v2x )
     !
  CASE( 11 ) ! 'wc'
     !
     CALL wcx( length, rho, grho, sx, v1x, v2x )
     !
  CASE( 12 ) ! 'pbexsr'
     !
     CALL pbex( length, rho, grho, 1, sx, v1x, v2x )
     !
     IF (exx_started_g) THEN
       CALL pbexsr( length, rho, grho, sxsr, v1xsr, v2xsr, screening_parameter_l )
       sx  = sx  - exx_fraction_g * sxsr
       v1x = v1x - exx_fraction_g * v1xsr
       v2x = v2x - exx_fraction_g * v2xsr
     ENDIF
     !
  CASE( 13 ) ! 'rPW86'
     !
     CALL rPW86( length, rho, grho, sx, v1x, v2x )
     !
  CASE( 16 ) ! 'C09x'
     !
     CALL c09x( length, rho, grho, sx, v1x, v2x )
     !
  CASE( 17 ) ! 'sogga'
     !
     CALL sogga( length, rho, grho, sx, v1x, v2x )
     !
  CASE( 19 ) ! 'pbeq2d'
     !
     CALL pbex( length, rho, grho, 4, sx, v1x, v2x )
     !
  CASE( 20 ) ! 'gau-pbe'
     !
     CALL pbex( length, rho, grho, 1, sx, v1x, v2x )
     IF (exx_started_g) THEN
       CALL pbexgau( length, rho, grho, sxsr, v1xsr, v2xsr, gau_parameter_l )
       sx  = sx  - exx_fraction_g * sxsr
       v1x = v1x - exx_fraction_g * v1xsr
       v2x = v2x - exx_fraction_g * v2xsr
     ENDIF
     !
  CASE( 21 ) ! 'pw86'
     !
     CALL pw86( length, rho, grho, sx, v1x, v2x )
     !
  CASE( 22 ) ! 'b86b'
     !
     CALL becke86b( length, rho, grho, sx, v1x, v2x )
     ! CALL b86b( rho, grho, 1, sx, v1x, v2x )
     !
  CASE( 23 ) ! 'optB88'
     !
     CALL pbex( length, rho, grho, 5, sx, v1x, v2x )
     !
  CASE( 24 ) ! 'optB86b'
     !
     CALL pbex( length, rho, grho, 6, sx, v1x, v2x )
     ! CALL b86b (rho, grho, 2, sx, v1x, v2x)
     !
  CASE( 25 ) ! 'ev93'
     !
     CALL pbex( length, rho, grho, 7, sx, v1x, v2x )
     !
  CASE( 26 ) ! 'b86r'
     !
     CALL b86b( length, rho, grho, 3, sx, v1x, v2x )
     !
  CASE( 27 ) ! 'cx13'
     !
     CALL cx13( length, rho, grho, sx, v1x, v2x )
     !
  CASE( 28 ) ! 'X3LYP'
     !
     CALL becke88( length, rho, grho, sx, v1x, v2x )
     CALL pbex( length, rho, grho, 1, sx_, v1x_, v2x_ )
     IF (exx_started_g) THEN
        sx  = REAL(0.765*0.709,DP) * sx
        v1x = REAL(0.765*0.709,DP) * v1x
        v2x = REAL(0.765*0.709,DP) * v2x
        sx  = sx  + REAL(0.235*0.709,DP) * sx_
        v1x = v1x + REAL(0.235*0.709,DP) * v1x_
        v2x = v2x + REAL(0.235*0.709,DP) * v2x_
     ENDIF
     !
  CASE( 29, 31 ) ! 'cx0'or `cx0p'
     !
     CALL cx13( length, rho, grho, sx, v1x, v2x )
     IF (exx_started_g) THEN
        sx  = (1.0_DP - exx_fraction_g) * sx
        v1x = (1.0_DP - exx_fraction_g) * v1x
        v2x = (1.0_DP - exx_fraction_g) * v2x
     ENDIF
     !
  CASE( 30 ) ! 'r860'
     !
     CALL rPW86( length, rho, grho, sx, v1x, v2x )
     !
     IF (exx_started_g) then
        sx  = (1.0_DP - exx_fraction_g) * sx
        v1x = (1.0_DP - exx_fraction_g) * v1x
        v2x = (1.0_DP - exx_fraction_g) * v2x
     ENDIF
     !
  CASE( 38 ) ! 'BR0'
     !
     CALL b86b( length, rho, grho, 3, sx, v1x, v2x )
     IF (exx_started_g) THEN
        sx  = (1.0_DP - exx_fraction_g) * sx
        v1x = (1.0_DP - exx_fraction_g) * v1x
        v2x = (1.0_DP - exx_fraction_g) * v2x
     ENDIF
     !
  CASE( 40 ) ! 'c090'
     !
     CALL c09x( length, rho, grho, sx, v1x, v2x )
     IF (exx_started_g) THEN
        sx  = (1.0_DP - exx_fraction_g) * sx
        v1x = (1.0_DP - exx_fraction_g) * v1x
        v2x = (1.0_DP - exx_fraction_g) * v2x
     ENDIF
     !
  CASE( 41 ) ! 'B86BPBEX'
     !
     CALL becke86b( length, rho, grho, sx, v1x, v2x )
     IF (exx_started_g) THEN
        sx  = (1.0_DP - exx_fraction_g) * sx
        v1x = (1.0_DP - exx_fraction_g) * v1x
        v2x = (1.0_DP - exx_fraction_g) * v2x
     ENDIF
     !
  CASE( 42 ) ! 'BHANDHLYP'
     !
     CALL becke88( length, rho, grho, sx, v1x, v2x )
     IF (exx_started_g) THEN
        sx  = (1.0_DP - exx_fraction_g) * sx
        v1x = (1.0_DP - exx_fraction_g) * v1x
        v2x = (1.0_DP - exx_fraction_g) * v2x
     ENDIF
     !
  CASE DEFAULT
     !
     sx = 0.0_DP
     v1x = 0.0_DP
     v2x = 0.0_DP
     !
  END SELECT
  !
  !
  ! ... CORRELATION
  !
  SELECT CASE( igcc_l )
  CASE( 1 )
     !
     CALL perdew86( length, rho, grho, sc, v1c, v2c )
     !
  CASE( 2 )
     !
     CALL ggac( length, rho, grho, sc, v1c, v2c )
     !
  CASE( 3 )
     !
     CALL glyp( length, rho, grho, sc, v1c, v2c )
     !
  CASE( 4 )
     !
     CALL pbec( length, rho, grho, 1, sc, v1c, v2c )
     !
  ! igcc_l == 5 (HCTH) is calculated together with case igcx_l=5
  ! igcc_l == 6 (meta-GGA) is treated in a different routine
  CASE( 7 ) !'B3LYP'
     !
     CALL glyp( length, rho, grho, sc, v1c, v2c )
     IF (exx_started_g) THEN
        sc  = 0.81_DP * sc
        v1c = 0.81_DP * v1c
        v2c = 0.81_DP * v2c
     ENDIF
     !
  CASE( 8 ) ! 'PBEsol'
     !
     CALL pbec( length, rho, grho, 2, sc, v1c, v2c )
     !
  ! igcc_l ==  9 set to 5, back-compatibility
  ! igcc_l == 10 set to 6, back-compatibility
  ! igcc_l == 11 M06L calculated in another routine
  CASE( 12 ) ! 'Q2D'
     !
     CALL pbec( length, rho, grho, 3, sc, v1c, v2c )
     !
  CASE( 13 ) !'X3LYP'
     !
     CALL glyp( length, rho, grho, sc, v1c, v2c )
     IF (exx_started_g) THEN
        sc  = 0.871_DP * sc
        v1c = 0.871_DP * v1c
        v2c = 0.871_DP * v2c
     ENDIF
     !
  CASE DEFAULT
     !
     sc = 0.0_DP
     v1c = 0.0_DP
     v2c = 0.0_DP
     !
  END SELECT
  !
  !
  WHERE ( rho <= small )
     sx  = 0.0_DP
     sc  = 0.0_DP
     v1x = 0.0_DP
     v2x = 0.0_DP
     v1c = 0.0_DP
     v2c = 0.0_DP
  END WHERE
  !
  IF ( igcx_l == 28 ) THEN
     DEALLOCATE( sx_        )
     DEALLOCATE( v1x_, v2x_ )
  ELSEIF ( (igcx_l == 12 .OR. igcx_l == 20) .AND. exx_started_g ) THEN
     DEALLOCATE( sxsr         )
     DEALLOCATE( v1xsr, v2xsr )
  ENDIF
  !
  RETURN
  !
END SUBROUTINE gcxc
!
!
!===============> SPIN <===============!
!
!-------------------------------------------------------------------------
SUBROUTINE gcx_spin( length, rho_in, grho2_in, sx_tot, v1x, v2x )
  !-----------------------------------------------------------------------
  !! Gradient corrections for exchange - Hartree a.u.
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: length
  !! Length of the input/output arrays
  REAL(DP), INTENT(IN),  DIMENSION(length,2) :: rho_in
  !! Up and down charge density
  REAL(DP), INTENT(IN),  DIMENSION(length,2) :: grho2_in
  !! Up and down gradient of the charge
  REAL(DP), INTENT(OUT), DIMENSION(length) :: sx_tot
  !! Energy
  REAL(DP), INTENT(OUT), DIMENSION(length,2) :: v1x
  !! Derivatives of exchange wr. rho
  REAL(DP), INTENT(OUT), DIMENSION(length,2) :: v2x
  !! Derivatives of exchange wr. grho
  !
  ! ... local variables
  !
  INTEGER :: is, iflag
  REAL(DP), DIMENSION(length,2) :: rho, grho2
  REAL(DP), DIMENSION(length,2) :: sx, null_v
  REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: sxsr
  REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: v1xsr, v2xsr
  !
  REAL(DP), PARAMETER :: small=1.E-10_DP
  REAL(DP), PARAMETER :: rho_trash=0.5_DP, grho2_trash=0.2_DP
  !
  !
  rho  = rho_in
  grho2 = grho2_in
  sx_tot = 0.0_DP
  null_v = 1.0_DP
  !
  DO is = 1, 2
     WHERE ( rho_in(:,is) <= small .OR. SQRT(ABS(grho2_in(:,is)) ) <= small )
        rho(:,is) = rho_trash
        grho2(:,is) = grho2_trash
        null_v(:,is) = 0.0_DP
     END WHERE
  ENDDO
  !
  IF ( igcx_l==12 .OR. igcx_l==20 .OR. igcx_l==28 ) THEN
     ALLOCATE( sxsr(length,2) )
     ALLOCATE( v1xsr(length,2), v2xsr(length,2) )
  ENDIF
  !
  ! ... exchange
  !
  SELECT CASE( igcx_l )
  CASE( 0 )
     !
     sx_tot = 0.0_DP
     v1x = 0.0_DP
     v2x = 0.0_DP
     !
  CASE( 1 )
     !
     CALL becke88_spin( length, rho, grho2, sx, v1x, v2x )
     !
     sx_tot = sx(:,1)*null_v(:,1) + sx(:,2)*null_v(:,2)
     !
  CASE( 2 )
     !
     rho = 2.0_DP * rho
     grho2 = 4.0_DP * grho2
     !
     CALL ggax( length, rho(:,1), grho2(:,1), sx(:,1), v1x(:,1), v2x(:,1) )
     CALL ggax( length, rho(:,2), grho2(:,2), sx(:,2), v1x(:,2), v2x(:,2) )
     !
     sx_tot = 0.5_DP * ( sx(:,1)*null_v(:,1) + sx(:,2)*null_v(:,2) )
     v2x = 2.0_DP * v2x
     !
  CASE( 3, 4, 8, 10, 12, 20, 23, 24, 25 )
     ! igcx_l=3: PBE, igcx_l=4: revised PBE, igcx_l=8: PBE0, igcx_l=10: PBEsol
     ! igcx_l=12: HSE,  igcx_l=20: gau-pbe, igcx_l=23: obk8, igcx_l=24: ob86, igcx_l=25: ev93
     !
     iflag = 1
     IF ( igcx_l== 4 ) iflag = 2
     IF ( igcx_l==10 ) iflag = 3
     IF ( igcx_l==23 ) iflag = 5
     IF ( igcx_l==24 ) iflag = 6
     IF ( igcx_l==25 ) iflag = 7
     !
     rho = 2.0_DP * rho
     grho2 = 4.0_DP * grho2
     !
     CALL pbex( length, rho(:,1), grho2(:,1), iflag, sx(:,1), v1x(:,1), v2x(:,1) )
     CALL pbex( length, rho(:,2), grho2(:,2), iflag, sx(:,2), v1x(:,2), v2x(:,2) )
     !
     sx_tot = 0.5_DP * ( sx(:,1)*null_v(:,1) + sx(:,2)*null_v(:,2) )
     v2x = 2.0_DP * v2x
     !
     IF ( igcx_l == 8 .AND. exx_started_g ) THEN
        !
        sx_tot = (1.0_DP - exx_fraction_g) * sx_tot
        v1x = (1.0_DP - exx_fraction_g) * v1x
        v2x = (1.0_DP - exx_fraction_g) * v2x
        !
     ELSEIF ( igcx_l == 12 .AND. exx_started_g ) THEN
        !
        CALL pbexsr( length, rho(:,1), grho2(:,1), sxsr(:,1), v1xsr(:,1), &
                                                        v2xsr(:,1), screening_parameter_l )
        CALL pbexsr( length, rho(:,2), grho2(:,2), sxsr(:,2), v1xsr(:,2), &
                                                        v2xsr(:,2), screening_parameter_l )
        !
        sx_tot = sx_tot - exx_fraction_g * 0.5_DP * ( sxsr(:,1)*null_v(:,1) + sxsr(:,2)*null_v(:,2) )
        v1x = v1x - exx_fraction_g * v1xsr
        v2x = v2x - exx_fraction_g * v2xsr * 2.0_DP
        !
     ELSEIF ( igcx_l == 20 .AND. exx_started_g ) THEN
        ! gau-pbe
        !CALL pbexgau_lsd( length, rho, grho2, sxsr, v1xsr, v2xsr, gau_parameter_l )
        CALL pbexgau( length, rho(:,1), grho2(:,1), sxsr(:,1), v1xsr(:,1), &
                                                              v2xsr(:,1), gau_parameter_l )
        CALL pbexgau( length, rho(:,2), grho2(:,2), sxsr(:,2), v1xsr(:,2), &
                                                              v2xsr(:,2), gau_parameter_l )
        !
        sx_tot = sx_tot - exx_fraction_g * 0.5_DP * ( sxsr(:,1)*null_v(:,1) + sxsr(:,2)*null_v(:,2) )
        v1x = v1x - exx_fraction_g * v1xsr
        v2x = v2x - exx_fraction_g * v2xsr * 2.0_DP
        !
     ENDIF
     !
  CASE( 9 )                    ! B3LYP
     !
     CALL becke88_spin( length, rho, grho2, sx, v1x, v2x )
     !
     sx_tot = sx(:,1)*null_v(:,1) + sx(:,2)*null_v(:,2)
     !
     IF ( exx_started_g ) THEN
        sx_tot = 0.72_DP * sx_tot
        v1x = 0.72_DP * v1x
        v2x = 0.72_DP * v2x
     ENDIF
     !
  CASE( 11 )                   ! 'Wu-Cohen'
     !
     rho = 2.0_DP * rho
     grho2 = 4.0_DP * grho2
     !
     CALL wcx( length, rho(:,1), grho2(:,1), sx(:,1), v1x(:,1), v2x(:,1) )
     CALL wcx( length, rho(:,2), grho2(:,2), sx(:,2), v1x(:,2), v2x(:,2) )
     !
     sx_tot = 0.5_DP * ( sx(:,1)*null_v(:,1) + sx(:,2)*null_v(:,2) )
     v2x = 2.0_DP * v2x
     !
  CASE( 13 )                   ! 'revised PW86 for vdw-df2'
     !
     rho = 2.0_DP * rho
     grho2 = 4.0_DP * grho2
     !
     CALL rPW86( length, rho(:,1), grho2(:,1), sx(:,1), v1x(:,1), v2x(:,1) )
     CALL rPW86( length, rho(:,2), grho2(:,2), sx(:,2), v1x(:,2), v2x(:,2) )
     !
     sx_tot = 0.5_DP * ( sx(:,1)*null_v(:,1) + sx(:,2)*null_v(:,2) )
     v2x = 2.0_DP * v2x
     !
  CASE( 16 )                   ! 'c09x for vdw-df-c09.'
     !
     rho = 2.0_DP * rho
     grho2 = 4.0_DP * grho2
     !
     CALL c09x( length, rho(:,1), grho2(:,1), sx(:,1), v1x(:,1), v2x(:,1) )
     CALL c09x( length, rho(:,2), grho2(:,2), sx(:,2), v1x(:,2), v2x(:,2) )
     !
     sx_tot = 0.5_DP * ( sx(:,1)*null_v(:,1) + sx(:,2)*null_v(:,2) )
     v2x = 2.0_DP * v2x
     !
  CASE( 21 )                   ! 'PW86'
     !
     rho = 2.0_DP * rho
     grho2 = 4.0_DP * grho2
     !
     CALL pw86( length, rho(:,1), grho2(:,1), sx(:,1), v1x(:,1), v2x(:,1) )
     CALL pw86( length, rho(:,2), grho2(:,2), sx(:,2), v1x(:,2), v2x(:,2) )
     !
     sx_tot = 0.5_DP * ( sx(:,1)*null_v(:,1) + sx(:,2)*null_v(:,2) )
     v2x = 2.0_DP * v2x
     !
  CASE( 22 )                   ! 'B86B'
     !
     rho = 2.0_DP * rho
     grho2 = 4.0_DP * grho2
     !
     CALL becke86b( length, rho(:,1), grho2(:,1), sx(:,1), v1x(:,1), v2x(:,1) )
     CALL becke86b( length, rho(:,2), grho2(:,2), sx(:,2), v1x(:,2), v2x(:,2) )
     !
     sx_tot = 0.5_DP * ( sx(:,1)*null_v(:,1) + sx(:,2)*null_v(:,2) )
     v2x = 2.0_DP * v2x
     !
   CASE( 26 )                  ! 'B86R for rev-vdW-DF2'
     !
     rho = 2.0_DP * rho
     grho2 = 4.0_DP * grho2
     !
     CALL b86b( length, rho(:,1), grho2(:,1), 3, sx(:,1), v1x(:,1), v2x(:,1) )
     CALL b86b( length, rho(:,2), grho2(:,2), 3, sx(:,2), v1x(:,2), v2x(:,2) )
     !
     sx_tot = 0.5_DP * ( sx(:,1)*null_v(:,1) + sx(:,2)*null_v(:,2) )
     v2x = 2.0_DP * v2x
     !
  CASE( 27 )                   ! 'cx13 for vdw-df-cx'
     !
     rho = 2.0_DP * rho
     grho2 = 4.0_DP * grho2
     !
     CALL cx13( length, rho(:,1), grho2(:,1), sx(:,1), v1x(:,1), v2x(:,1) )
     CALL cx13( length, rho(:,2), grho2(:,2), sx(:,2), v1x(:,2), v2x(:,2) )
     !
     sx_tot = 0.5_DP * ( sx(:,1)*null_v(:,1) + sx(:,2)*null_v(:,2) )
     v2x = 2.0_DP * v2x
     !
  CASE( 28 )                   ! X3LYP
     !
     CALL becke88_spin( length, rho, grho2, sx, v1x, v2x )
     !
     rho = 2.0_DP * rho
     grho2 = 4.0_DP * grho2
     !
     CALL pbex( length, rho(:,1), grho2(:,1), 1, sxsr(:,1), v1xsr(:,1), v2xsr(:,1) )
     CALL pbex( length, rho(:,2), grho2(:,2), 1, sxsr(:,2), v1xsr(:,2), v2xsr(:,2) )
     !
     sx_tot = 0.5_DP * ( sxsr(:,1)*null_v(:,1) + sxsr(:,2)*null_v(:,2) ) * 0.235_DP + &
                       (   sx(:,1)*null_v(:,1) +   sx(:,2)*null_v(:,2) ) * 0.765_DP
     v1x = v1xsr * 0.235_DP + v1x * 0.765_DP
     v2x = v2xsr * 0.235_DP * 2.0_DP + v2x * 0.765_DP
     !
     IF ( exx_started_g ) THEN
        sx_tot = 0.709_DP * sx_tot
        v1x = 0.709_DP * v1x
        v2x = 0.709_DP * v2x
     ENDIF
     !
  CASE( 29, 31 )               ! 'cx0 for vdw-df-cx0' or `cx0p for vdW-DF-cx0p'
     !
     rho = 2.0_DP * rho
     grho2 = 4.0_DP * grho2
     !
     CALL cx13( length, rho(:,1), grho2(:,1), sx(:,1), v1x(:,1), v2x(:,1) )
     CALL cx13( length, rho(:,2), grho2(:,2), sx(:,2), v1x(:,2), v2x(:,2) )
     !
     sx_tot = 0.5_DP * ( sx(:,1)*null_v(:,1) + sx(:,2)*null_v(:,2) )
     v2x = 2.0_DP * v2x
     !
     IF ( exx_started_g ) THEN
        sx_tot = (1.0_DP - exx_fraction_g) * sx_tot
        v1x = (1.0_DP - exx_fraction_g) * v1x
        v2x = (1.0_DP - exx_fraction_g) * v2x
     ENDIF
     !
  CASE( 30 )                   ! 'R860' = 'rPW86-0' for vdw-df2-0'
     !
     rho = 2.0_DP * rho
     grho2 = 4.0_DP * grho2
     !
     CALL rPW86( length, rho(:,1), grho2(:,1), sx(:,1), v1x(:,1), v2x(:,1) )
     CALL rPW86( length, rho(:,2), grho2(:,2), sx(:,2), v1x(:,2), v2x(:,2) )
     !
     sx_tot = 0.5_DP * ( sx(:,1)*null_v(:,1) + sx(:,2)*null_v(:,2) )
     v2x = 2.0_DP * v2x
     !
     IF ( exx_started_g ) THEN
        sx_tot = (1.0_DP - exx_fraction_g) * sx_tot
        v1x = (1.0_DP - exx_fraction_g) * v1x
        v2x = (1.0_DP - exx_fraction_g) * v2x
     ENDIF
     !
  CASE( 38 )                  ! 'br0 for vdw-df2-BR0' etc
     !
     rho = 2.0_DP * rho
     grho2 = 4.0_DP * grho2
     !
     CALL b86b( length, rho(:,1), grho2(:,1), 3, sx(:,1), v1x(:,1), v2x(:,1) )
     CALL b86b( length, rho(:,2), grho2(:,2), 3, sx(:,2), v1x(:,2), v2x(:,2) )     
     !
     sx_tot = 0.5_DP * ( sx(:,1)*null_v(:,1) + sx(:,2)*null_v(:,2) )
     v2x = 2.0_DP * v2x
     !
     IF ( exx_started_g ) THEN
        sx_tot = (1.0_DP - exx_fraction_g) * sx_tot
        v1x = (1.0_DP - exx_fraction_g) * v1x
        v2x = (1.0_DP - exx_fraction_g) * v2x
     ENDIF  
     !
  CASE( 40 )                  ! 'c090 for vdw-df-c090' etc
     !
     rho = 2.0_DP * rho
     grho2 = 4.0_DP * grho2
     !
     CALL c09x( length, rho(:,1), grho2(:,1), sx(:,1), v1x(:,1), v2x(:,1) )
     CALL c09x( length, rho(:,2), grho2(:,2), sx(:,2), v1x(:,2), v2x(:,2) )  
     !
     sx_tot = 0.5_DP * ( sx(:,1)*null_v(:,1) + sx(:,2)*null_v(:,2) )
     v2x = 2.0_DP * v2x
     !
     IF ( exx_started_g ) THEN
        sx_tot = (1.0_DP - exx_fraction_g) * sx_tot
        v1x = (1.0_DP - exx_fraction_g) * v1x
        v2x = (1.0_DP - exx_fraction_g) * v2x
     ENDIF
     !
  CASE( 41 )                 ! B86X for B86BPBEX hybrid
     !
     rho = 2.0_DP * rho
     grho2 = 4.0_DP * grho2
     !
     CALL becke86b( length, rho(:,1), grho2(:,1), sx(:,1), v1x(:,1), v2x(:,1) )
     CALL becke86b( length, rho(:,2), grho2(:,2), sx(:,2), v1x(:,2), v2x(:,2) )
     !
     sx_tot = 0.5_DP * ( sx(:,1)*null_v(:,1) + sx(:,2)*null_v(:,2) )
     v2x = 2.0_DP * v2x
     !
     IF ( exx_started_g ) THEN
        sx_tot = (1.0_DP - exx_fraction_g) * sx_tot
        v1x = (1.0_DP - exx_fraction_g) * v1x
        v2x = (1.0_DP - exx_fraction_g) * v2x
     ENDIF
     !
  CASE( 42 )                ! B88X for BHANDHLYP
     !
     rho = 2.0_DP * rho
     grho2 = 4.0_DP * grho2
     !
     CALL becke88( length, rho(:,1), grho2(:,1), sx(:,1), v1x(:,1), v2x(:,1) )
     CALL becke88( length, rho(:,2), grho2(:,2), sx(:,2), v1x(:,2), v2x(:,2) )
     !
     sx_tot = 0.5_DP * ( sx(:,1)*null_v(:,1) + sx(:,2)*null_v(:,2) )
     v2x = 2.0_DP * v2x
     !
     IF ( exx_started_g ) THEN
        sx_tot = (1.0_DP - exx_fraction_g) * sx_tot
        v1x = (1.0_DP - exx_fraction_g) * v1x
        v2x = (1.0_DP - exx_fraction_g) * v2x
     ENDIF
     !
  ! case igcx_l == 5 (HCTH) and 6 (OPTX) not implemented
  ! case igcx_l == 7 (meta-GGA) must be treated in a separate call to another
  ! routine: needs kinetic energy density in addition to rho and grad rho
  CASE DEFAULT
     !
     CALL errore( 'gcx_spin', 'not implemented', igcx_l )
     !
  END SELECT
  !
  DO is = 1, 2
    v1x(:,is) = v1x(:,is) * null_v(:,is)
    v2x(:,is) = v2x(:,is) * null_v(:,is)
  ENDDO
  !
  IF ( igcx_l==12 .OR. igcx_l==20 .OR. igcx_l==28 ) THEN
     DEALLOCATE( sxsr )
     DEALLOCATE( v1xsr, v2xsr )
  ENDIF
  !
  RETURN
  !
END SUBROUTINE gcx_spin
!
!
!---------------------------------------------------------------------
SUBROUTINE gcc_spin( length, rho_in, zeta_io, grho_in, sc, v1c, v2c )
  !-------------------------------------------------------------------
  !! Gradient corrections for correlations - Hartree a.u. 
  !! Implemented:  Perdew86, GGA (PW91), PBE
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: length
  !! the length of the I/O arrays
  REAL(DP), INTENT(IN), DIMENSION(length) :: rho_in
  !! the total charge
  REAL(DP), INTENT(INOUT), DIMENSION(length) :: zeta_io
  !! the magnetization
  REAL(DP), INTENT(IN), DIMENSION(length) :: grho_in
  !! the gradient of the charge squared
  REAL(DP), INTENT(OUT), DIMENSION(length) :: sc
  !! correlation energies
  REAL(DP), INTENT(OUT), DIMENSION(length,2) :: v1c
  !! derivative of correlation wr. rho
  REAL(DP), INTENT(OUT), DIMENSION(length) :: v2c
  !! derivatives of correlation wr. grho
  !
  ! ... local variables
  !
  REAL(DP), DIMENSION(length) :: rho, zeta, grho
  REAL(DP), DIMENSION(length) :: null_v
  REAL(DP), PARAMETER :: small=1.E-10_DP, epsr=1.E-6_DP
  REAL(DP), PARAMETER :: rho_trash=0.5_DP, zeta_trash=0.2_DP, grho_trash=0.05_DP
  !
  !
  rho = rho_in
  grho = grho_in
  null_v = 1.0_DP
  !
  WHERE ( ABS(zeta_io) <= 1.0_DP ) &
     zeta_io = SIGN( MIN( ABS(zeta_io), (1.0_DP-epsr) ), zeta_io )
  !
  zeta = zeta_io
  !
  WHERE ( ABS(zeta_io)>1.0_DP .OR. rho<=small .OR. SQRT(ABS(grho))<=small )
     rho = rho_trash
     grho = grho_trash
     zeta = zeta_trash
     null_v = 0.0_DP
  END WHERE
  !
  SELECT CASE( igcc_l )
  CASE( 0 )
     !
     sc  = 0.0_DP
     v1c = 0.0_DP
     v2c = 0.0_DP
     !
  CASE( 1 )
     !
     CALL perdew86_spin( length, rho, zeta, grho, sc, v1c, v2c )
     !
  CASE( 2 )
     !
     CALL ggac_spin( length, rho, zeta, grho, sc, v1c, v2c )
     !
  CASE( 4 )
     !
     CALL pbec_spin( length, rho, zeta, grho, 1, sc, v1c, v2c )
     !
  CASE( 8 )
     !
     CALL pbec_spin( length, rho, zeta, grho, 2, sc, v1c, v2c )
     !
  CASE DEFAULT
     !
     CALL errore( 'lsda_functionals (gcc_spin)', 'not implemented', igcc_l )
     !
  END SELECT
  !
  sc = sc * null_v
  v1c(:,1) = v1c(:,1) * null_v
  v1c(:,2) = v1c(:,2) * null_v
  v2c = v2c * null_v
  !
  !
  RETURN
  !
END SUBROUTINE gcc_spin
!
!---------------------------------------------------------------------------
SUBROUTINE gcc_spin_more( length, rho_in, grho_in, grho_ud_in, &
                                               sc, v1c, v2c, v2c_ud )
  !-------------------------------------------------------------------------
  !! Gradient corrections for exchange and correlation.
  !
  !! * Exchange:
  !!    * Becke88;
  !!    * GGAX.
  !! * Correlation:
  !!    * Perdew86;
  !!    * Lee, Yang & Parr;
  !!    * GGAC.
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: length
  !! length of the I/O arrays
  REAL(DP), INTENT(IN), DIMENSION(length,2) :: rho_in
  !! the total charge
  REAL(DP), INTENT(IN), DIMENSION(length,2) :: grho_in
  !! the gradient of the charge squared
  REAL(DP), INTENT(IN), DIMENSION(length) :: grho_ud_in
  !! gradient off-diagonal term up-down
  REAL(DP), INTENT(OUT), DIMENSION(length) :: sc
  !! correlation energies
  REAL(DP), INTENT(OUT), DIMENSION(length,2) :: v1c
  !! derivative of correlation wr. rho
  REAL(DP), INTENT(OUT), DIMENSION(length,2) :: v2c
  !! derivative of correlation wr. grho
  REAL(DP), INTENT(OUT), DIMENSION(length) :: v2c_ud
  !! derivative of correlation wr. grho, off-diag. term
  !
  ! ... local variables
  !
  REAL(DP), DIMENSION(length,2) :: rho, grho
  REAL(DP), DIMENSION(length) :: grho_ud, null_v
  REAL(DP), PARAMETER :: small=1.E-20_DP
  REAL(DP), PARAMETER :: rho_trash=0.5_DP, grho_trash=0.2_DP, &
                         grhoud_trash=0.05_DP
  !
  rho = rho_in
  grho = grho_in
  grho_ud = grho_ud_in
  sc  = 0.0_DP
  v1c = 0.0_DP
  v2c = 0.0_DP
  v2c_ud = 0.0_DP
  null_v = 1.0_DP
  !
  WHERE ( rho_in(:,1)+rho_in(:,2) < small )
     rho(:,1) = rho_trash
     rho(:,2) = rho_trash
     grho(:,1) = grho_trash
     grho(:,2) = grho_trash
     grho_ud(:) = grhoud_trash
     null_v(:) = 0.0_DP
  END WHERE
  !
  CALL lsd_glyp( length, rho, grho, grho_ud, sc, v1c, v2c, v2c_ud )
  !
  SELECT CASE( igcc_l )
  CASE( 3 )
     !
     ! ... void
     !
  CASE( 7 )
     !
     IF ( exx_started_g ) THEN
        sc = 0.81_DP * sc
        v1c = 0.81_DP * v1c
        v2c = 0.81_DP * v2c
        v2c_ud = 0.81_DP * v2c_ud
     ENDIF
     !
  CASE( 13 )
     !
     IF ( exx_started_g ) THEN
        sc = 0.871_DP * sc
        v1c = 0.871_DP * v1c
        v2c = 0.871_DP * v2c
        v2c_ud = 0.871_DP * v2c_ud
     ENDIF
     !
  CASE DEFAULT
     !
     CALL errore( " gcc_spin_more ", " gradient correction not implemented ", 1 )
     !
  END SELECT
  !
  sc = sc * null_v
  v1c(:,1) = v1c(:,1) * null_v(:)
  v2c(:,1) = v2c(:,1) * null_v(:)
  v1c(:,2) = v1c(:,2) * null_v(:)
  v2c(:,2) = v2c(:,2) * null_v(:)
  v2c_ud = v2c_ud * null_v
  !
  !
  RETURN
  !
END SUBROUTINE gcc_spin_more
!
!
END MODULE xc_gga
