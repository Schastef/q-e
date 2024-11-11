!----------------------------------------------------------
SUBROUTINE write_epsilon(npe,drhoscfh)
  !----------------------------------------------
  ! F. Macheda (2024)
  ! ---------------------------------------------
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : ionode
  USE fft_base,             ONLY : dfftp 
  USE units_ph,             ONLY : iurhoun 
  USE gvect,                ONLY : gg
  USE noncollin_module,     ONLY : nspin_mag
  USE fft_interfaces,       ONLY : fwfft
  USE dv_of_drho_lr,        ONLY : dv_of_drho
  USE constants,            ONLY : fpi,e2
  USE cell_base,            ONLY : tpiba
  USE qpoint,               ONLY : xq
  USE control_ph,           ONLY : extpot
  use mp,                   ONLY : mp_sum
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE control_lr,           ONLY : lmacro
  !
  IMPLICIT NONE
  INTEGER, INTENT(in) :: npe
  complex(DP), intent(in) :: drhoscfh (dfftp%nnr,nspin_mag,npe)
  complex(DP), allocatable :: dvscf_toprint(:,:,:), drho_toprint(:,:,:)
  INTEGER, EXTERNAL :: find_free_unit
  INTEGER :: is, ii
  COMPLEX(DP) :: epsm1, eps, chi, chi0, chiRPA, barchi
  REAL(DP) :: vq
  COMPLEX(DP) :: drhoaux, dvaux
  INTEGER :: my_rnk
  !
  allocate (dvscf_toprint ( dfftp%nnr, nspin_mag , npe))
  allocate (drho_toprint ( dfftp%nnr, nspin_mag , npe))
  !
  call zcopy (dfftp%nnr*nspin_mag,drhoscfh,1,dvscf_toprint,1)
  call zcopy (dfftp%nnr*nspin_mag,drhoscfh,1,drho_toprint,1)
  !
  call dv_of_drho (dvscf_toprint(1,1,1)) !nlcc cannot be applied since they depend on the perturbation, that we are taking as a scalar here
  !
  CALL fwfft ('Rho', dvscf_toprint(:,1,1), dfftp)
  CALL fwfft ('Rho', drho_toprint(:,1,1), dfftp)
  !
  drho_toprint = drho_toprint / extpot
  dvscf_toprint = dvscf_toprint / extpot
  !
  drhoaux=0d0
  dvaux=0d0
  IF (gg(1) < 1d-8) THEN
    drhoaux = drho_toprint(dfftp%nl(1),1,1)
    dvaux   = dvscf_toprint(dfftp%nl(1),1,1)
  ENDIF
  CALL mp_sum(drhoaux,intra_bgrp_comm)
  CALL mp_sum(dvaux,intra_bgrp_comm)
  !
  IF (ionode) THEN
    !      
    iurhoun=find_free_unit()
    OPEN(unit=iurhoun, file="drhodv.dat")
    !
    WRITE(iurhoun, '(a)') "#  Re(drhoind),Im(drhoind),Re(dVscf),Im(dVscf)"
    !
    WRITE(iurhoun,'(4f18.12)') real(drhoaux), aimag(drhoaux),&
                               real(dvaux)  , aimag(dvaux)
    !
    vq = fpi*e2/tpiba**2/dot_product(xq,xq)
    !
    IF (.NOT. lmacro) THEN
      !
      WRITE(iurhoun, '(a)') "#  epsm1"
      !
      epsm1 = 1d0+dvaux
      WRITE(iurhoun,'(2f18.12)') REAL(epsm1), AIMAG(epsm1)
      !
      WRITE(iurhoun, '(a)') "#  1/epsm1"
      WRITE(iurhoun,'(2f18.12)') REAL(1d0/epsm1), AIMAG(1d0/epsm1)
      !
      WRITE(iurhoun, '(a)') "#  chi"
      !
      chi = drhoaux
      WRITE(iurhoun,'(2f18.12)') REAL(chi), AIMAG(chi)
      !
      WRITE(iurhoun, '(a)') "#  (1+vq*chi)"
      WRITE(iurhoun,'(2f18.12)') REAL(1d0+vq*chi), AIMAG(1d0+vq*chi)
      !
    ELSE
      !
      chi0 = drhoaux/(1d0+dvaux)
      WRITE(iurhoun, '(a)') "#  bar chi0"
      WRITE(iurhoun,'(2f18.12)') REAL(chi0), AIMAG(chi0)
      !
      WRITE(iurhoun, '(a)') "#  bar chi"
      !
      barchi = drhoaux
      WRITE(iurhoun,'(2f18.12)') REAL(barchi), AIMAG(barchi)
      !
      WRITE(iurhoun, '(a)') "# 1/epsm1"
      WRITE(iurhoun,'(2f18.12)') REAL(1d0-vq*drhoaux), AIMAG(1d0-vq*drhoaux)
      !
    ENDIF
    !
  ENDIF
  !
  CLOSE(iurhoun)
  DEALLOCATE (dvscf_toprint)
  DEALLOCATE (drho_toprint)
  !
END SUBROUTINE write_epsilon
!----------------------------------------------------------
SUBROUTINE write_drhoun
  !----------------------------------------------
  ! F. Macheda (2024)
  ! ---------------------------------------------
  USE kinds,          ONLY : DP
  USE io_global,      ONLY : stdout, ionode
  USE units_ph,       ONLY : iudumpdrho, iudumpeff 
  USE modes,          ONLY : u
  USE ions_base,      ONLY : tau, nat, ityp
  USE gvect,          ONLY : mill, ig_l2g
  USE qpoint,         ONLY : xq
  USE cell_base,      ONLY : tpiba, omega
  USE uspp_param,     ONLY : upf
  USE constants,      ONLY : tpi, e2 
  USE mp_bands,       ONLY : intra_bgrp_comm
  USE control_ph,     ONLY : current_iq, extpot
  USE eqv,            ONLY : vlocq
  USE gvecs,          ONLY : ngms, ngms_g
  USE dynmat,         ONLY : dyn
  use mp,             ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  INTEGER  :: irr, ipert, imode0, is
  INTEGER  :: ig, i, j
  INTEGER  :: ipol, jpol, na, nb, mode, nt
  INTEGER, ALLOCATABLE :: itmp_mill(:,:)
  INTEGER  :: mill_igg(3)
  COMPLEX(DP), ALLOCATABLE :: drho0(:,:)
  COMPLEX(DP), ALLOCATABLE :: aux2_g(:,:,:)
  COMPLEX(DP) :: Im_i=(0._dp,1._dp)
  COMPLEX(DP), ALLOCATABLE :: phase(:)
  REAL(DP)  :: zval
  REAL(DP)  :: arg 
  CHARACTER(LEN=6), EXTERNAL :: int_to_char
  complex(DP), ALLOCATABLE :: aux1 (:,:,:)  
  COMPLEX(DP) :: gtau
  INTEGER, EXTERNAL :: find_free_unit
  complex(DP) :: phi (3, nat), dynsym(3*nat,3*nat), tmpphi(3,nat) 
  REAL(DP) :: absq
  REAL(KIND = DP) :: sa(3, 3)
  REAL(KIND = DP) :: sb(3, 3)
  REAL(KIND = DP) :: srr(3, 3)
  !
  absq = SQRT(DOT_PRODUCT(xq,xq))
  !
  ALLOCATE(phase(nat))
  !
  ALLOCATE(aux2_g(3,nat,ngms_g))
  ALLOCATE(itmp_mill(3,ngms_g))
  !
  itmp_mill =0
  !
  DO ig=1,ngms
    !
    itmp_mill( 1, ig_l2g( ig ) ) = mill(1, ig )
    itmp_mill( 2, ig_l2g( ig ) ) = mill(2, ig )
    itmp_mill( 3, ig_l2g( ig ) ) = mill(3, ig )
    !
  ENDDO
  !
  CALL mp_sum( itmp_mill , intra_bgrp_comm )
  !
  !
  phi = (0.0d0, 0.0d0)
  !
  DO na=1,nat
    !
    DO ipol=1,3
      !
      i=(na-1)*3+ipol
      DO nb=1,nat
        !
        DO jpol=1,3
          !
          j=(nb-1)*3+jpol
          phi(ipol,na) = phi(ipol,na) + dyn (i, j) * CONJG(u (i, j)) 
          !
        ENDDO
        !
      ENDDO
      !
    ENDDO
    !
  ENDDO
  !
  CALL symvectorq(nat,phi)
  !
  IF ( ionode ) then
    !
    iudumpdrho = find_free_unit()
    iudumpeff = find_free_unit()
    !
    iudumpdrho = find_free_unit()
    OPEN(unit=iudumpdrho,file='drho'//TRIM(int_to_char(current_iq))//'.dat')
    WRITE(iudumpdrho,*) '# Re(drhox),Im(drhox),Re(drhoy),Im(drhoy),Re(drhoz),Im(drhoz) '
    iudumpeff = find_free_unit()
    OPEN(unit=iudumpeff,file='effcharge'//TRIM(int_to_char(current_iq))//'.dat')
    WRITE(iudumpeff,*) '# Re(Zx),Im(Zx),Re(Zy), Im(Zy),Re(Zz),Im(Zz) '
    !
    DO ig=1,ngms_g 
      !
      IF(ABS(itmp_mill(1,ig)).le.1d-8.and.&
         ABS(itmp_mill(2,ig)).le.1d-8.and.&
         ABS(itmp_mill(3,ig)).le.1d-8)THEN 
        !
        DO na = 1, nat
          !
          arg = (xq(1) * tau(1,na) + &
                 xq(2) * tau(2,na) + &
                 xq(3) * tau(3,na))*tpi
          !
          phase(na)= CMPLX( COS( arg ),  SIN( arg ) ,kind=DP)
          !
        ENDDO
        !
        phi = phi / extpot
        !
        DO na=1,nat
          !
          zval=upf(ityp(na))%zp
          !
          write(iudumpdrho,'(6f18.12)'),&
              (-sqrt(e2)/omega*(phi(ipol,na)*phase(na)+&
              Im_i*xq(ipol)*tpiba*zval), ipol=1,3)
          write(iudumpeff,'(6f18.12)')&
              (1d0/Im_i/(absq*tpiba)*(phi(ipol,na)*phase(na)+&
              Im_i*xq(ipol)*tpiba*zval), ipol=1,3)
          !
        ENDDO 
        !
      ENDIF
      !
    ENDDO
    !
    CLOSE(iudumpdrho)
    !
  ENDIF
  !
  DEALLOCATE(phase)
  DEALLOCATE(itmp_mill)
  !
  WRITE(stdout, *) "------------------------------------------------------------------------"
  WRITE(stdout, *) " The code is printing the induced charges"
  WRITE(stdout, *) " Please refer to:"
  WRITE(stdout, *) " Macheda F., Barone P. & Mauri, F. (2024), "
  WRITE(stdout, *) " First-principles calculations of dynamical Born effective charges, quadrupoles,"
  WRITE(stdout, *) " and higher order terms from the charge response in large semiconducting and metallic systems"
  WRITE(stdout, *) " Physical Review B, 110, 094306. https://doi.org/10.1103/PhysRevB.110.094306"
  WRITE(stdout, *) "------------------------------------------------------------------------"
  CONTAINS
   !--------------------------------------------------------------------------
   SUBROUTINE symvectorq( nat, vect )
     !-----------------------------------------------------------------------
     !! Symmetrize a function \(f(i,na)\) (e.g. the forces in cartesian axis),
     !! where \(i\) is the cartesian component, \(na\) the atom index.
     USE kinds,          ONLY : DP
     USE cell_base,      ONLY : at, bg
     USE symm_base,      ONLY : s, nsym, t_rev, irt, invs 
     USE lr_symm_base,   ONLY : minus_q, irotmq, nsymq, rtau
     USE qpoint,         ONLY : xq
     !
     IMPLICIT NONE
     !
     INTEGER, INTENT(IN) :: nat
     !! number of atoms
     COMPLEX(DP), INTENT(INOUT) :: vect(3,nat)
     !! vector function to symmetrize
     !
     ! ... local variables
     !
     INTEGER :: na, isym, nar
     COMPLEX(DP), ALLOCATABLE :: phi(:,:), work(:), phip(:,:)
     REAL(DP) :: arg, fase
     COMPLEX(DP) :: faseq(48)
     INTEGER :: iflb (nat), isymq, kpol, sna, irot
     !
     IF (nsym == 1) RETURN
     !
     ALLOCATE (phip(3,nat))
     ALLOCATE (phi(3,nat))
     ALLOCATE (work(3))
     !
     ! bring vector to crystal axis
     !
     DO na = 1, nat
        phi(:,na) =  vect(1,na)*at(1,:) + &
                     vect(2,na)*at(2,:) + &
                     vect(3,na)*at(3,:)
     END DO
     !
     !    If no other symmetry is present we quit here
     !
     if ( (nsymq == 1) .and. (.not.minus_q) ) return
     !
     !    Then we impose the symmetry q -> -q+G if present
     !
     if (minus_q) then
        do na = 1, nat
              do ipol = 1, 3
                    work(:) = (0.d0, 0.d0)
                    sna = irt (irotmq, na)
                    arg = 0.d0
                    do kpol = 1, 3
                       arg = arg + (xq (kpol) * (-rtau (kpol, irotmq, na) ) ) 
                    enddo
                    arg = arg * tpi
                    fase = CMPLX(cos (arg), sin (arg) ,kind=DP)
                    do kpol = 1, 3
                          work (ipol) = work (ipol) + &
                               s (ipol, kpol, irotmq) &
                               * phi (kpol, sna) * fase
                    enddo
                    phip (ipol, na) = (phi (ipol, na) + &
                         CONJG( work (ipol) ) ) * 0.5d0
              enddo
           enddo
        phi = phip
     endif
     !
     !    Here we symmetrize with respect to the small group of q
     !
     if (nsymq == 1) return
     !
     iflb (:) = 0
     do na = 1, nat
           if (iflb (na) == 0) then
              work(:) = (0.d0, 0.d0)
              do isymq = 1, nsymq
                 irot = isymq
                 sna = irt (irot, na)
                 arg = 0.d0
                 do ipol = 1, 3
                    arg = arg + (xq (ipol) * (-rtau (ipol, irot, na) ) )
                 enddo
                 arg = arg * tpi
                 faseq (isymq) = CMPLX(cos (arg), sin (arg) ,kind=DP)
                 do ipol = 1, 3
                       do kpol = 1, 3
                             IF (t_rev(isymq)==1) THEN
                                work (ipol) = work (ipol) + &
                                     s (ipol, kpol, irot) &
                              * CONJG(phi (kpol, sna) * faseq (isymq))
                             ELSE
                                work (ipol) = work (ipol) + &
                                     s (ipol, kpol, irot) &
                                    * phi (kpol, sna) * faseq (isymq)
                             ENDIF
                       enddo
                 enddo
              enddo
              do isymq = 1, nsymq
                 irot = isymq
                 sna = irt (irot, na)
                 do ipol = 1, 3
                       phi (ipol, sna) = (0.d0, 0.d0)
                       do kpol = 1, 3
                             IF (t_rev(isymq)==1) THEN
                                phi(ipol,sna)=phi(ipol,sna) &
                                + s(ipol,kpol,invs(irot))&
                                  * CONJG(work (kpol)*faseq (isymq))
                             ELSE
                                phi(ipol,sna)=phi(ipol,sna) &
                                + s(ipol,kpol,invs(irot))&
                                  * work (kpol) * CONJG(faseq (isymq) )
                             ENDIF
                       enddo
                 enddo
                 iflb (sna) = 1
              enddo
           endif
     enddo
     phi (:, :) = phi (:, :) / DBLE(nsymq)
     !
     ! bring vector back to cartesian axis
     !
     DO na = 1, nat
        vect(:,na) = phi(1,na)*bg(:,1) + &
                     phi(2,na)*bg(:,2) + &
                     phi(3,na)*bg(:,3)
     END DO
     !
     DEALLOCATE (phi)
     DEALLOCATE (work)
     !
   END SUBROUTINE symvectorq
   !
END SUBROUTINE write_drhoun
!
SUBROUTINE Vaeps_dvloc(pot, mode, ind_ig)
  !
  USE kinds,          ONLY : DP
  USE fft_base,       ONLY : dffts
  USE modes,          ONLY : u
  USE ions_base,      ONLY : nat, ityp, ntyp => nsp
  USE gvect,          ONLY : g, mill, eigts1, eigts2, eigts3, ngm
  USE qpoint,         ONLY : xq, eigqts
  USE cell_base,      ONLY : tpiba, tpiba2, omega
  USE uspp_param,     ONLY : upf
  USE gvecs,          ONLY : ngms
  !
  IMPLICIT NONE
  COMPLEX(DP), INTENT(INOUT) :: pot
  INTEGER, INTENT(IN) :: mode
  INTEGER, INTENT(IN) :: ind_ig
  !! Index to be passed
  REAL (DP), ALLOCATABLE :: vlocq(:,:)  ! ngm, ntyp)
  INTEGER :: na, mu, ig, itmp, nt
  INTEGER, ALLOCATABLE :: nl_d(:)
  complex(DP) :: gtau, gu, fact, u1, u2, u3, gu0
  complex(DP) , allocatable :: aux1 (:)
  REAL(DP) :: zval

  ALLOCATE( nl_d(dffts%ngm) )
  nl_d  = dffts%nl
  allocate (aux1(dffts%nnr)) !aux1 is dvlocin
  ALLOCATE(vlocq(ngm,ntyp))

  do nt = 1, ntyp
    zval=upf(nt)%zp
    CALL setlocq_coul (xq, zval, tpiba2, ngm, g, omega, vlocq(:,nt))
  enddo

  aux1 = 0.0d0

  do na = 1, nat
     fact = tpiba * (0.d0, -1.d0) * eigqts (na)
     mu = 3 * (na - 1)
     if ( abs (u (mu + 1, mode) ) + abs (u (mu + 2, mode) ) + &
          abs (u (mu + 3, mode) ) > 1.0d-12) then
        nt = ityp (na)
        u1 = u (mu + 1, mode)
        u2 = u (mu + 2, mode)
        u3 = u (mu + 3, mode)
        gu0 = xq (1) * u1 + xq (2) * u2 + xq (3) * u3
        do ig = 1, ngms
           gtau = eigts1 (mill(1,ig), na) * eigts2 (mill(2,ig), na) * &
                  eigts3 (mill(3,ig), na)
           gu = gu0 + g (1, ig) * u1 + g (2, ig) * u2 + g (3, ig) * u3
           aux1 (dffts%nl (ig) ) = aux1 (dffts%nl (ig) ) + vlocq (ig, nt) &
                * gu * fact * gtau
        enddo
     endif
  enddo
  !
  pot = pot - aux1(ind_ig)
  DEALLOCATE(aux1)

END SUBROUTINE Vaeps_dvloc
!----------------------------------------------------------------------
subroutine setlocq_coul (xq, zp, tpiba2, ngm, g, omega, vloc)
 !----------------------------------------------------------------------
 !! Fourier transform of the Coulomb potential - For all-electron
 !! calculations, in specific cases only, for testing purposes.
 !
 USE kinds, ONLY: DP
 USE constants, ONLY : fpi, e2, eps8
 implicit none
 !
 integer, intent(in) :: ngm
 real(DP) :: xq (3), zp, tpiba2, omega, g(3,ngm)
 real(DP), intent (out) :: vloc(ngm)
 !
 real(DP) :: g2a
 integer :: ig

 do ig = 1, ngm
  g2a = (xq (1) + g (1, ig) ) **2 + (xq (2) + g (2, ig) ) **2 + &
        (xq (3) + g (3, ig) ) **2
  if (g2a < eps8) then
       vloc (ig) = 0.d0
  else
       vloc (ig) = - fpi * zp *e2 / omega / tpiba2 / g2a
  endif
 enddo

end subroutine setlocq_coul
!----------------------------------------------------------------------
!----------------------------------------------------------
SUBROUTINE init_rho(npe,drhoscf,drhoscfh,iq_dummy)
  !----------------------------------------------
  ! F. Macheda (2024)
  ! ---------------------------------------------
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : ionode
  USE fft_base,             ONLY : dfftp, dffts 
  USE noncollin_module,     ONLY : nspin_mag
  USE fft_interfaces,       ONLY : fft_interpolate
  USE save_ph,              ONLY : tmp_dir_save
  USE output,               ONLY : fildrho
  USE units_ph,             ONLY : iudrho, lrdrho
  USE io_files,             ONLY : prefix, diropn
  USE qpoint,               ONLY : xq
  USE cell_base,            ONLY : at
  USE gvecs,                ONLY : doublegrid
  USE dfile_autoname,       ONLY : dfile_name
  !
  IMPLICIT NONE
  INTEGER, INTENT(in) :: npe
  complex(DP), intent(inout) :: drhoscfh (dfftp%nnr,nspin_mag,npe)
  complex(DP), intent(inout) :: drhoscf (dffts%nnr,nspin_mag,npe)
  INTEGER, INTENT(in) :: iq_dummy
  !
  INTEGER :: ipert, is
  character(len=256) :: filename
  logical :: exst
  !
  do ipert = 1, npe
    if (fildrho.ne.' ') then
      IF (ionode) THEN
        INQUIRE(UNIT = iudrho, OPENED = exst)
        IF (exst) CLOSE (UNIT = iudrho, STATUS='keep')
        filename = dfile_name(xq, at, fildrho, TRIM(tmp_dir_save)//prefix, generate=.true., index_q=iq_dummy)
        CALL diropn (iudrho, filename, lrdrho, exst)
      ENDIF ! ionode
      !     
      call davcio_drho (drhoscfh(1,1,ipert), lrdrho, iudrho, 1, -1)
      !
    endif
  enddo
  !
  if (doublegrid) then
     do is = 1, nspin_mag
        do ipert = 1, npe
           call fft_interpolate (dfftp, drhoscfh(:,is,ipert), dffts, drhoscf(:,is,ipert))
        enddo
     enddo
  else
     call zcopy (npe*nspin_mag*dfftp%nnr, drhoscfh, 1, drhoscf, 1)
  endif
  !
END SUBROUTINE init_rho
!----------------------------------------------------------

