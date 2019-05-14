!
! Copyright (C) 2004-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-------------------------------------------------------------------
module funct
!-------------------------------------------------------------------
! This module contains data defining the DFT functional in use
! and a number of functions and subroutines to manage them.
! Data are PRIVATE and are accessed and set only by function calls.
! Basic drivers to compute XC quantities are also included.
!
!  setting routines:   set_dft_from_name (previously which_dft)
!                      set_dft_from_indices
!                      enforce_input_dft
!                      start_exx
!                      stop_exx
!                      set_finite_size_volume
!  retrieve functions: get_dft_name, get_dft_short, get_dft_long
!                      get_iexch
!                      get_icorr
!                      get_igcx
!                      get_igcc
!                      get_exx_fraction
!                      write_dft_name
!  logical functions:  dft_is_gradient
!                      dft_is_meta
!                      dft_is_hybrid
!                      dft_is_nonlocc
!                      exx_is_active
!                      dft_has_finite_size_correction
!
!  XC computation drivers: xc, xc_spin, gcxc, gcx_spin, gcc_spin, gcc_spin_more
!  derivatives of XC-gga computation drivers: dgcxc, dgcxc_spin
!
  USE io_global, ONLY: stdout
  USE kinds,     ONLY: DP
  IMPLICIT NONE
  PRIVATE
  SAVE
  ! subroutines/functions managing dft name and indices
  PUBLIC  :: set_dft_from_indices, set_dft_from_name
  PUBLIC  :: enforce_input_dft, write_dft_name
  PUBLIC  :: get_dft_name, get_dft_short, get_dft_long,&
             get_nonlocc_name
  PUBLIC  :: get_iexch, get_icorr, get_igcx, get_igcc, get_meta, get_inlc
  PUBLIC  :: dft_is_gradient, dft_is_meta, dft_is_hybrid, dft_is_nonlocc, igcc_is_lyp
  PUBLIC  :: set_auxiliary_flags

  ! additional subroutines/functions for hybrid functionals
  PUBLIC  :: start_exx, stop_exx, get_exx_fraction, exx_is_active
  PUBLIC  :: set_exx_fraction, dft_force_hybrid
  PUBLIC  :: set_screening_parameter, get_screening_parameter
  PUBLIC  :: set_gau_parameter, get_gau_parameter

  ! additional subroutines/functions for finite size corrections
  PUBLIC  :: dft_has_finite_size_correction, set_finite_size_volume
  ! rpa specific
  PUBLIC  :: init_dft_exxrpa, enforce_dft_exxrpa

  ! driver subroutines computing XC
  PUBLIC  :: init_lda_xc, init_gga_xc
  PUBLIC  :: tau_xc , tau_xc_spin
  PUBLIC  :: tau_xc_array, tau_xc_array_spin
  PUBLIC  :: nlc
  !
  ! PRIVATE variables defining the DFT functional
  !
  PRIVATE :: dft, iexch, icorr, igcx, igcc, imeta, inlc
  PRIVATE :: qe_or_libxc
  PRIVATE :: discard_input_dft
  PRIVATE :: isgradient, ismeta, ishybrid
  PRIVATE :: exx_fraction, exx_started
  PRIVATE :: has_finite_size_correction, &
             finite_size_cell_volume,  finite_size_cell_volume_set
  !
  character (len=25) :: dft = 'not set'
  !
  ! "dft" is the exchange-correlation functional label, as set by the user,
  ! using either set_dft_from_name or set_dft_from_indices. It can contain
  ! either the short names or a series of the keywords listed below.
  ! All operations on names are case-insensitive.
  !
  !           short name       complete name       Short description
  !              "pz"    = "sla+pz"            = Perdew-Zunger LDA
  !              "bp"    = "b88+p86"           = Becke-Perdew grad.corr.
  !              "pw91"  = "sla+pw+ggx+ggc"    = PW91 (aka GGA)
  !              "blyp"  = "sla+b88+lyp+blyp"  = BLYP
  !              "pbe"   = "sla+pw+pbx+pbc"    = PBE
  !              "revpbe"= "sla+pw+rpb+pbc"    = revPBE (Zhang-Yang)
  !              "pw86pbe" = "sla+pw+pw86+pbc" = PW86 exchange + PBE correlation
  !              "b86bpbe" = "sla+pw+b86b+pbc" = B86b exchange + PBE correlation
  !              "pbesol"= "sla+pw+psx+psc"    = PBEsol
  !              "q2d"   = "sla+pw+q2dx+q2dc"  = PBEQ2D
  !              "hcth"  = "nox+noc+hcth+hcth" = HCTH/120
  !              "olyp"  = "nox+lyp+optx+blyp" = OLYP
  !              "wc"    = "sla+pw+wcx+pbc"    = Wu-Cohen
  !              "sogga" = "sla+pw+sox+pbec"   = SOGGA
  !              "optbk88"="sla+pw+obk8+p86"   = optB88
  !              "optb86b"="sla+pw+ob86+p86"   = optB86
  !              "ev93"  = "sla+pw+evx+nogc"   = Engel-Vosko
  !              "tpss"  = "sla+pw+tpss+tpss"  = TPSS Meta-GGA
  !              "m06l"  = "nox+noc+m6lx+m6lc" = M06L Meta-GGA
  !              "tb09"  = "sla+pw+tb09+tb09"  = TB09 Meta-GGA
  !              "pbe0"  = "pb0x+pw+pb0x+pbc"  = PBE0
  !              "b86bx" = "pb0x+pw+b86x+pbc"  = B86bPBE hybrid
  !              "bhahlyp" = "pb0x+pw+b88x+blyp"  = Becke half-and-half LYP
  !              "hse"   = "sla+pw+hse+pbc"    = Heyd-Scuseria-Ernzerhof (HSE 06, see note below)
  !              "b3lyp" = "b3lp+b3lp+b3lp+b3lp"= B3LYP
  !              "b3lypv1r"    = "b3lp+b3lpv1r+b3lp+b3lp"= B3LYP-VWN1-RPA
  !              "x3lyp" = "x3lp+x3lp+x3lp+x3lp"= X3LYP
  !              "vwn-rpa"     = "sla+vwn-rpa" = VWN LDA using vwn1-rpa parametriz
  !              "gaupbe"= "sla+pw+gaup+pbc"   = Gau-PBE (also "gaup")
  !              "vdw-df"       ="sla+pw+rpb +vdw1"   = vdW-DF1
  !              "vdw-df2"      ="sla+pw+rw86+vdw2"   = vdW-DF2
  !              "vdw-df-c09"   ="sla+pw+c09x+vdw1"   = vdW-DF-C09
  !              "vdw-df2-c09"  ="sla+pw+c09x+vdw2"   = vdW-DF2-C09
  !              "vdw-df-obk8"  ="sla+pw+obk8+vdw1"   = vdW-DF-obk8 (optB88-vdW)
  !              "vdw-df-ob86"  ="sla+pw+ob86+vdw1"   = vdW-DF-ob86 (optB86b-vdW)
  !              "vdw-df2-b86r" ="sla+pw+b86r+vdw2"   = vdW-DF2-B86R (rev-vdw-df2)
  !              "vdw-df-cx"    ="sla+pw+cx13+vdW1"   = vdW-DF-cx
  !              "vdw-df-cx0"    ="sla+pw+cx13+vdW1+HF/4"   = vdW-DF-cx-0
  !              "vdw-df2-0"     ="sla+pw+rw86+vdw2+HF/4"   = vdW-DF2-0
  !              "vdw-df2-br0"  ="sla+pw+b86r+vdW2+HF/4"   = vdW-DF2-b86r-0
  !              "vdw-df-c090"   ="sla+pw+c09x+vdw1+HF/4"   = vdW-DF-C09-0
  !              "vdw-df-x"     ="sla+pw+????+vdwx"   = vdW-DF-x, reserved Thonhauser, not implemented
  !              "vdw-df-y"     ="sla+pw+????+vdwy"   = vdW-DF-y, reserved Thonhauser, not implemented
  !              "vdw-df-z"     ="sla+pw+????+vdwz"   = vdW-DF-z, reserved Thonhauser, not implemented
  !              "rvv10" = "sla+pw+rw86+pbc+vv10"     = rVV10
  !
  ! Any nonconflicting combination of the following keywords is acceptable:
  !
  ! Exchange:    "nox"    none                           iexch=0
  !              "sla"    Slater (alpha=2/3)             iexch=1 (default)
  !              "sl1"    Slater (alpha=1.0)             iexch=2
  !              "rxc"    Relativistic Slater            iexch=3
  !              "oep"    Optimized Effective Potential  iexch=4
  !              "hf"     Hartree-Fock                   iexch=5
  !              "pb0x"   (Slater*0.75+HF*0.25)          iexch=6 for PBE0 and vdW-DF-cx0 and vdW-DF2-0 etc
  !              "b3lp"   B3LYP(Slater*0.80+HF*0.20)     iexch=7
  !              "kzk"    Finite-size corrections        iexch=8
  !              "x3lp"   X3LYP(Slater*0.782+HF*0.218)   iexch=9
  !              "kli"    KLI aproximation for exx       iexch=10
  !
  ! Correlation: "noc"    none                           icorr=0
  !              "pz"     Perdew-Zunger                  icorr=1 (default)
  !              "vwn"    Vosko-Wilk-Nusair              icorr=2
  !              "lyp"    Lee-Yang-Parr                  icorr=3
  !              "pw"     Perdew-Wang                    icorr=4
  !              "wig"    Wigner                         icorr=5
  !              "hl"     Hedin-Lunqvist                 icorr=6
  !              "obz"    Ortiz-Ballone form for PZ      icorr=7
  !              "obw"    Ortiz-Ballone form for PW      icorr=8
  !              "gl"     Gunnarson-Lunqvist             icorr=9
  !              "kzk"    Finite-size corrections        icorr=10
  !              "vwn-rpa" Vosko-Wilk-Nusair, alt param  icorr=11
  !              "b3lp"   B3LYP (0.19*vwn+0.81*lyp)      icorr=12
  !              "b3lpv1r"  B3LYP-VWN-1-RPA
  !                         (0.19*vwn_rpa+0.81*lyp)      icorr=13
  !              "x3lp"   X3LYP (0.129*vwn_rpa+0.871*lyp)icorr=14
  !
  ! Gradient Correction on Exchange:
  !              "nogx"   none                           igcx =0 (default)
  !              "b88"    Becke88 (beta=0.0042)          igcx =1
  !              "ggx"    Perdew-Wang 91                 igcx =2
  !              "pbx"    Perdew-Burke-Ernzenhof exch    igcx =3
  !              "rpb"    revised PBE by Zhang-Yang      igcx =4
  !              "hcth"   Cambridge exch, Handy et al    igcx =5
  !              "optx"   Handy's exchange functional    igcx =6
  !              "pb0x"   PBE0 (PBE exchange*0.75)       igcx =8
  !              "b3lp"   B3LYP (Becke88*0.72)           igcx =9
  !              "psx"    PBEsol exchange                igcx =10
  !              "wcx"    Wu-Cohen                       igcx =11
  !              "hse"    HSE screened exchange          igcx =12
  !              "rw86"   revised PW86                   igcx =13
  !              "pbe"    same as PBX, back-comp.        igcx =14
  !              "c09x"   Cooper 09                      igcx =16
  !              "sox"    sogga                          igcx =17
  !              "q2dx"   Q2D exchange grad corr         igcx =19
  !              "gaup"   Gau-PBE hybrid exchange        igcx =20
  !              "pw86"   Perdew-Wang (1986) exchange    igcx =21
  !              "b86b"   Becke (1986) exchange          igcx =22
  !              "obk8"   optB88  exchange               igcx =23
  !              "ob86"   optB86b exchange               igcx =24
  !              "evx"    Engel-Vosko exchange           igcx =25
  !              "b86r"   revised Becke (b86b)           igcx =26
  !              "cx13"   consistent exchange            igcx =27
  !              "x3lp"   X3LYP (Becke88*0.542 +
  !                              Perdew-Wang91*0.167)    igcx =28
  !              "cx0"    vdW-DF-cx+HF/4 (cx13-0)        igcx =29 
  !              "r860"   rPW86+HF/4 (rw86-0)            igcx =30 (for DF0)
  !              "cx0p"   vdW-DF-cx+HF/5 (cx13-0p)       igcx =31 
  !              "ahcx"   vdW-DF-cx based not yet in use igcx =32 reserved PH
  !              "ahf2"   vdW-DF2 based not yet in use   igcx =33 reserved PH
  !              "ahpb"   PBE based not yet in use       igcx =34 reserved PH
  !              "ahps"   PBE-sol based not in use       igcx =35 reserved PH
  !              "cx14"   Exporations                    igcx =36 reserved PH
  !              "cx15"   Exporations                    igcx =37 reserved PH
  !              "br0"    vdW-DF2-b86r+HF/4 (b86r-0)     igcx =38 
  !              "cx16"   Exporations                    igcx =39 reserved PH
  !              "c090"   vdW-DF-c09+HF/4 (c09-0)        igcx =40 
  !              "b86x"   B86b exchange * 0.75           igcx =41
  !              "b88x"   B88 exchange * 0.50            igcx =42
  !
  ! Gradient Correction on Correlation:
  !              "nogc"   none                           igcc =0 (default)
  !              "p86"    Perdew86                       igcc =1
  !              "ggc"    Perdew-Wang 91 corr.           igcc =2
  !              "blyp"   Lee-Yang-Parr                  igcc =3
  !              "pbc"    Perdew-Burke-Ernzenhof corr    igcc =4
  !              "hcth"   Cambridge corr, Handy et al    igcc =5
  !              "b3lp"   B3LYP (Lee-Yang-Parr*0.81)     igcc =7
  !              "psc"    PBEsol corr                    igcc =8
  !              "pbe"    same as PBX, back-comp.        igcc =9
  !              "q2dc"   Q2D correlation grad corr      igcc =12
  !              "x3lp"   X3LYP (Lee-Yang-Parr*0.871)    igcc =13
  !
  ! Meta-GGA functionals
  !              "tpss"   TPSS Meta-GGA                  imeta=1
  !              "m6lx"   M06L Meta-GGA                  imeta=2
  !              "tb09"   TB09 Meta-GGA                  imeta=3
  !              "+meta"  activate MGGA even without MGGA-XC   imeta=4
  !              "scan"   SCAN Meta-GGA                  imeta=5
  !              "sca0"   SCAN0  Meta-GGA                imeta=6
  !
  ! Van der Waals functionals (nonlocal term only)
  !              "nonlc"  none                           inlc =0 (default)
  !              "vdw1"   vdW-DF1                        inlc =1
  !              "vdw2"   vdW-DF2                        inlc =2
  !              "vv10"   rVV10                          inlc =3
  !              "vdwx"   vdW-DF-x                       inlc =4, reserved Thonhauser, not implemented
  !              "vdwy"   vdW-DF-y                       inlc =5, reserved Thonhauser, not implemented
  !              "vdwz"   vdW-DF-z                       inlc =6, reserved Thonhauser, not implemented
  !
  ! Meta-GGA with Van der Waals
  !              "rvv10-scan" rVV10 (with b=15.7) and scan inlc=3 (PRX 6, 041005 (2016))
  !
  ! Note: as a rule, all keywords should be unique, and should be different
  ! from the short name, but there are a few exceptions.
  !
  ! References:
  !              pz      J.P.Perdew and A.Zunger, PRB 23, 5048 (1981)
  !              vwn     S.H.Vosko, L.Wilk, M.Nusair, Can.J.Phys. 58,1200(1980)
  !              vwn1-rpa S.H.Vosko, L.Wilk, M.Nusair, Can.J.Phys. 58,1200(1980)
  !              wig     E.P.Wigner, Trans. Faraday Soc. 34, 67 (1938)
  !              hl      L.Hedin and B.I.Lundqvist, J. Phys. C4, 2064 (1971)
  !              gl      O.Gunnarsson and B.I.Lundqvist, PRB 13, 4274 (1976)
  !              pw      J.P.Perdew and Y.Wang, PRB 45, 13244 (1992)
  !              obpz    G.Ortiz and P.Ballone, PRB 50, 1391 (1994)
  !              obpw    as above
  !              b88     A.D.Becke, PRA 38, 3098 (1988)
  !              p86     J.P.Perdew, PRB 33, 8822 (1986)
  !              pw86    J.P.Perdew, PRB 33, 8800 (1986)
  !              b86b    A.D.Becke, J.Chem.Phys. 85, 7184 (1986)
  !              ob86    Klimes, Bowler, Michaelides, PRB 83, 195131 (2011)
  !              b86r    I. Hamada, Phys. Rev. B 89, 121103(R) (2014)
  !              pbe     J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996)
  !              pw91    J.P.Perdew and Y. Wang, PRB 46, 6671 (1992)
  !              blyp    C.Lee, W.Yang, R.G.Parr, PRB 37, 785 (1988)
  !              hcth    Handy et al, JCP 109, 6264 (1998)
  !              olyp    Handy et al, JCP 116, 5411 (2002)
  !              revPBE  Zhang and Yang, PRL 80, 890 (1998)
  !              pbesol  J.P. Perdew et al., PRL 100, 136406 (2008)
  !              q2d     L. Chiodo et al., PRL 108, 126402 (2012)
  !              rw86    Eamonn D. Murray et al, J. Chem. Theory Comput. 5, 2754 (2009)
  !              wc      Z. Wu and R. E. Cohen, PRB 73, 235116 (2006)
  !              kzk     H.Kwee, S. Zhang, H. Krakauer, PRL 100, 126404 (2008)
  !              pbe0    J.P.Perdew, M. Ernzerhof, K.Burke, JCP 105, 9982 (1996)
  !              hse     Heyd, Scuseria, Ernzerhof, J. Chem. Phys. 118, 8207 (2003)
  !                      Heyd, Scuseria, Ernzerhof, J. Chem. Phys. 124, 219906 (2006).
  !              b3lyp   P.J. Stephens,F.J. Devlin,C.F. Chabalowski,M.J. Frisch
  !                      J.Phys.Chem 98, 11623 (1994)
  !              x3lyp   X. Xu, W.A Goddard III, PNAS 101, 2673 (2004)
  !              vdW-DF       M. Dion et al., PRL 92, 246401 (2004)
  !                           T. Thonhauser et al., PRL 115, 136402 (2015)
  !              vdW-DF2      Lee et al., Phys. Rev. B 82, 081101 (2010)
  !              rev-vdW-DF2  I. Hamada, Phys. Rev. B 89, 121103(R) (2014)
  !              vdW-DF-cx    K. Berland and P. Hyldgaard, PRB 89, 035412 (2014)
  !              vdW-DF-cx0   K. Berland, Y. Jiao, J.-H. Lee, T. Rangel, J. B. Neaton and P. Hyldgaard,
  !                           J. Chem. Phys. 146, 234106 (2017)
  !              vdW-DF-cx0p  Y. Jiao, E. Schröder and P. Hyldgaard, 
  !                           J. Chem. Phys. 148, 194115 (2018)
  !              vdW-DF-obk8  Klimes et al, J. Phys. Cond. Matter, 22, 022201 (2010)
  !              vdW-DF-ob86  Klimes et al, Phys. Rev. B, 83, 195131 (2011)
  !              c09x    V. R. Cooper, Phys. Rev. B 81, 161104(R) (2010)
  !              tpss    J.Tao, J.P.Perdew, V.N.Staroverov, G.E. Scuseria,
  !                      PRL 91, 146401 (2003)
  !              tb09    F Tran and P Blaha, Phys.Rev.Lett. 102, 226401 (2009)
  !              scan    J Sun, A Ruzsinszky and J Perdew, PRL 115, 36402 (2015)
  !              scan0   K Hui and J-D. Chai, JCP 144, 44114 (2016)
  !              sogga   Y. Zhao and D. G. Truhlar, JCP 128, 184109 (2008)
  !              m06l    Y. Zhao and D. G. Truhlar, JCP 125, 194101 (2006)
  !              gau-pbe J.-W. Song, K. Yamashita, K. Hirao JCP 135, 071103 (2011)
  !              rVV10   R. Sabatini et al. Phys. Rev. B 87, 041108(R) (2013)
  !              ev93     Engel-Vosko, Phys. Rev. B 47, 13164 (1993)
  !
  ! NOTE ABOUT HSE: there are two slight deviations with respect to the HSE06
  ! functional as it is in Gaussian code (that is considered as the reference
  ! in the chemistry community):
  ! - The range separation in Gaussian is precisely 0.11 bohr^-1,
  !   instead of 0.106 bohr^-1 in this implementation
  ! - The gradient scaling relation is a bit more complicated
  !   [ see: TM Henderson, AF Izmaylov, G Scalmani, and GE Scuseria,
  !          J. Chem. Phys. 131, 044108 (2009) ]
  ! These two modifications accounts only for a 1e-5 Ha difference for a
  ! single He atom. Info by Fabien Bruneval
  !
  integer, parameter:: notset = -1
  !
  ! internal indices for exchange-correlation
  !    iexch: type of exchange
  !    icorr: type of correlation
  !    igcx:  type of gradient correction on exchange
  !    igcc:  type of gradient correction on correlation
  !    inlc:  type of non local correction on correlation
  !    imeta: type of meta-GGA
  integer :: iexch = notset
  integer :: icorr = notset
  integer :: igcx  = notset
  integer :: igcc  = notset
  integer :: imeta = notset
  integer :: inlc  = notset
  !
  ! Switches to decide between qe (0) and libxc (1) routines for each single xc term
  ! (only for LDA part at the moment)
  ! PROVISIONAL: probably it should be put as an environment variable or something
#if defined(__LIBXC)
  integer :: qe_or_libxc(1:2) = 1
#else
  integer :: qe_or_libxc(1:2) = 0
#endif
  !
  real(DP):: exx_fraction = 0.0_DP
  real(DP):: screening_parameter = 0.0_DP
  real(DP):: gau_parameter = 0.0_DP
  logical :: islda       = .false.
  logical :: isgradient  = .false.
  logical :: ismeta      = .false.
  logical :: ishybrid    = .false.
  logical :: isnonlocc   = .false.
  logical :: exx_started = .false.
  logical :: has_finite_size_correction = .false.
  logical :: finite_size_cell_volume_set = .false.
  real(DP):: finite_size_cell_volume = notset
  logical :: discard_input_dft = .false.
  !
  integer, parameter:: nxc=8, ncc=10, ngcx=42, ngcc=12, nmeta=6, ncnl=6
  character (len=4) :: exc, corr, gradx, gradc, meta, nonlocc
  dimension :: exc (0:nxc), corr (0:ncc), gradx (0:ngcx), gradc (0:ngcc), &
               meta(0:nmeta), nonlocc (0:ncnl)

  data exc  / 'NOX', 'SLA', 'SL1', 'RXC', 'OEP', 'HF', 'PB0X', 'B3LP', 'KZK' /
  data corr / 'NOC', 'PZ', 'VWN', 'LYP', 'PW', 'WIG', 'HL', 'OBZ', &
              'OBW', 'GL' , 'KZK' /

  data gradx / 'NOGX', 'B88', 'GGX', 'PBX',  'RPB', 'HCTH', 'OPTX',&
               'xxxx', 'PB0X', 'B3LP','PSX', 'WCX', 'HSE', 'RW86', 'PBE', &
               'xxxx', 'C09X', 'SOX', 'xxxx', 'Q2DX', 'GAUP', 'PW86', 'B86B', &
               'OBK8', 'OB86', 'EVX', 'B86R', 'CX13', 'X3LP', &
               'CX0', 'R860', 'CX0P', 'AHCX', 'AHF2', &
               'AHPB', 'AHPS', 'CX14', 'CX15', 'BR0', 'CX16', 'C090', &
               'B86X', 'B88X'/

  data gradc / 'NOGC', 'P86', 'GGC', 'BLYP', 'PBC', 'HCTH', 'NONE',&
               'B3LP', 'PSC', 'PBE', 'xxxx', 'xxxx', 'Q2DC' /

  data meta  / 'NONE', 'TPSS', 'M06L', 'TB09', 'META', 'SCAN', 'SCA0' /

  data nonlocc/'NONE', 'VDW1', 'VDW2', 'VV10', 'VDWX', 'VDWY', 'VDWZ' /

#if defined(__LIBXC)
  integer :: libxc_major=0, libxc_minor=0, libxc_micro=0
  public :: libxc_major, libxc_minor, libxc_micro, get_libxc_version
#endif

CONTAINS
  !-----------------------------------------------------------------------
  subroutine set_dft_from_name( dft_ )
    !-----------------------------------------------------------------------
    !
    ! translates a string containing the exchange-correlation name
    ! into internal indices iexch, icorr, igcx, igcc, inlc, imeta
    !
    implicit none
    character(len=*), intent(in) :: dft_
    integer :: len, l, i
    character (len=50):: dftout
    logical :: dft_defined = .false.
    character (len=1), external :: capital
    integer ::  save_iexch, save_icorr, save_igcx, save_igcc, save_meta, save_inlc
    !
    ! Exit if set to discard further input dft
    !
    if ( discard_input_dft ) return
    !
    ! save current status of XC indices
    !
    save_iexch = iexch
    save_icorr = icorr
    save_igcx  = igcx
    save_igcc  = igcc
    save_meta  = imeta
    save_inlc  = inlc
    !
    ! convert to uppercase
    !
    len = len_trim(dft_)
    dftout = ' '
    do l = 1, len
       dftout (l:l) = capital (dft_(l:l) )
    enddo
    !
    ! ----------------------------------------------
    ! FIRST WE CHECK ALL THE SHORT NAMES
    ! Note: comparison is done via exact matching
    ! ----------------------------------------------
    !
    ! special cases : PZ  (LDA is equivalent to PZ)
    IF (('PZ' .EQ. TRIM(dftout) ).OR.('LDA' .EQ. TRIM(dftout) )) THEN
       dft_defined = set_dft_values(1,1,0,0,0,0)

    ! speciale cases : PW ( LDA with PW correlation )
    ELSE IF ( 'PW' .EQ. TRIM(dftout)) THEN
    dft_defined = set_dft_values(1,4,0,0,0,0)
    ! special cases : VWN-RPA
    else IF ('VWN-RPA' .EQ. TRIM(dftout) ) THEN
       dft_defined = set_dft_values(1,11,0,0,0,0)

    ! special cases : OEP no GC part (nor LDA...) and no correlation by default
    else IF ('OEP' .EQ. TRIM(dftout) ) THEN
       dft_defined = set_dft_values(4,0,0,0,0,0)

    else if ("KLI" .EQ. TRIM(dftout)) then
      dft_defined = set_dft_values(10,0,0,0,0,0)

    ! special cases : HF no GC part (nor LDA...) and no correlation by default
    else IF ('HF' .EQ. TRIM(dftout) ) THEN
       dft_defined = set_dft_values(5,0,0,0,0,0)

    else if ('PBE' .EQ. TRIM(dftout) ) then
    ! special case : PBE
       dft_defined = set_dft_values(1,4,3,4,0,0)
    !special case : B88
    else if ('B88' .EQ. TRIM(dftout) ) then
      dft_defined = set_dft_values(1,1,1,0,0,0)
    ! special case : BP = B88 + P86
    else if ('BP'.EQ. TRIM(dftout) ) then
       dft_defined = set_dft_values(1,1,1,1,0,0)

    ! special case : PW91 = GGX + GGC
    else if ('PW91'.EQ. TRIM(dftout) ) then
       dft_defined = set_dft_values(1,4,2,2,0,0)

    elseif ( 'REVPBE' .EQ. TRIM(dftout) ) then
    ! special case : revPBE
       dft_defined = set_dft_values(1,4,4,4,0,0)

    else if ('PBESOL'.EQ. TRIM(dftout) ) then
    ! special case : PBEsol
       dft_defined = set_dft_values(1,4,10,8,0,0)

    ! special cases : BLYP (note, BLYP=>B88)
    else IF ('BLYP' .EQ. TRIM(dftout) ) THEN
       dft_defined = set_dft_values(1,3,1,3,0,0)

    else if ('OPTBK88' .EQ. TRIM(dftout)) then
    ! Special case optB88
       dft_defined = set_dft_values(1,4,23,1,0,0)

    else if ('OPTB86B' .EQ. TRIM(dftout)) then
    ! Special case optB86b
       dft_defined = set_dft_values(1,4,24,1,0,0)

    else if ('PBC'.EQ. TRIM(dftout) ) then
    ! special case : PBC  = PW + PBC
       dft_defined = set_dft_values(1,4,0,4,0,0)

    ! special case : HCTH
    else if ('HCTH'.EQ. TRIM(dftout)) then
       dft_defined = set_dft_values(0,0,5,5,0,0)
       call errore('set_dft_from_name','HCTH yields suspicious results',1)

    ! special case : OLYP = OPTX + LYP
    else if ('OLYP'.EQ. TRIM(dftout)) then
       dft_defined = set_dft_values(0,3,6,3,0,0)
       call errore('set_dft_from_name','OLYP yields suspicious results',1)

    else if ('WC' .EQ. TRIM(dftout) ) then
    ! special case : Wu-Cohen
       dft_defined = set_dft_values(1,4,11,4,0,0)

    elseif ('PW86PBE' .EQ. TRIM(dftout) ) then
       ! special case : PW86PBE
       dft_defined = set_dft_values(1,4,21,4,0,0)

    elseif ('B86BPBE' .EQ. TRIM(dftout) ) then
       ! special case : B86BPBE
       dft_defined = set_dft_values(1,4,22,4,0,0)

    else if ('PBEQ2D' .EQ. TRIM(dftout) .OR. 'Q2D'.EQ. TRIM(dftout) ) then
    ! special case : PBEQ2D
       dft_defined = set_dft_values(1,4,19,12,0,0)

    ! special case : SOGGA = SOX + PBEc
    else if ('SOGGA' .EQ. TRIM(dftout) ) then
       dft_defined = set_dft_values(1,4,17,4,0,0)

    ! special case : Engel-Vosko
    else if ( 'EV93' .EQ. TRIM(dftout) ) THEN
       dft_defined = set_dft_values(1,4,25,0,0,0)

    else if ('RPBE' .EQ. TRIM(dftout) ) then
    ! special case : RPBE
         call errore('set_dft_from_name', &
     &   'RPBE (Hammer-Hansen-Norskov) not implemented (revPBE is)',1)

    else if ('PBE0'.EQ. TRIM(dftout) ) then
    ! special case : PBE0
       dft_defined = set_dft_values(6,4,8,4,0,0)

    else if ('B86BPBEX'.EQ. TRIM(dftout) ) then
    ! special case : B86BPBEX
       dft_defined = set_dft_values(6,4,41,4,0,0)

    else if ('BHAHLYP'.EQ. TRIM(dftout).OR.'BHANDHLYP'.EQ.TRIM(dftout)) then
    ! special case : B86BPBEX
       dft_defined = set_dft_values(6,4,42,3,0,0)

   else if ('HSE' .EQ. TRIM( dftout) ) then
    ! special case : HSE
       dft_defined = set_dft_values(1,4,12,4,0,0)

   else if ( 'GAUP' .EQ. TRIM(dftout) .OR. 'GAUPBE' .EQ. TRIM(dftout) ) then
    ! special case : GAUPBE
       dft_defined = set_dft_values(1,4,20,4,0,0)

    else if ('VDW-DF' .EQ. TRIM(dftout)) then
    ! Special case vdW-DF
       dft_defined = set_dft_values(1,4,4,0,1,0)

    else if ('VDW-DF2' .EQ. TRIM(dftout) ) then
     ! Special case vdW-DF2
        dft_defined = set_dft_values(1,4,13,0,2,0)

    else if ('VDW-DF-C09'  .EQ. TRIM(dftout) ) then
    ! Special case vdW-DF with C09 exchange
       dft_defined = set_dft_values(1,4,16,0,1,0)

    else if ('VDW-DF2-C09' .EQ. TRIM(dftout) ) then
      ! Special case vdW-DF2 with C09 exchange
         dft_defined = set_dft_values(1,4,16,0,2,0)

    else if ('VDW-DF-OBK8' .EQ. TRIM(dftout)) then
    ! Special case vdW-DF-obk8, or vdW-DF + optB88
       dft_defined = set_dft_values(1,4,23,0,1,0)

    else if ('VDW-DF-OB86' .EQ. TRIM(dftout) ) then
    ! Special case vdW-DF-ob86, or vdW-DF + optB86
       dft_defined = set_dft_values(1,4,24,0,1,0)

    else if ('VDW-DF2-B86R' .EQ. TRIM(dftout) ) then
    ! Special case vdW-DF2 with B86R
       dft_defined = set_dft_values(1,4,26,0,2,0)

      else if ('VDW-DF-CX' .EQ. TRIM(dftout)) then
      ! Special case vdW-DF-CX
         dft_defined = set_dft_values(1,4,27,0,1,0)

      else if ('VDW-DF-CX0' .EQ. TRIM(dftout) ) then
      ! Special case vdW-DF-CX0
         dft_defined = set_dft_values(6,4,29,0,1,0)

    else if ('VDW-DF-CX0P' .EQ. TRIM(dftout) ) then
    ! Special case vdW-DF-CX0P
       dft_defined = set_dft_values(6,4,31,0,1,0)

      else if ('VDW-DF2-0' .EQ. TRIM(dftout) ) then
      ! Special case vdW-DF2-0
         dft_defined = set_dft_values(6,4,30,0,2,0)

      else if ('VDW-DF2-BR0' .EQ. TRIM(dftout) ) then
      ! Special case vdW-DF2-BR0
         dft_defined = set_dft_values(6,4,38,0,2,0)

      else if ('VDW-DF-C090' .EQ. TRIM(dftout) ) then
      ! Special case vdW-DF-C090
         dft_defined = set_dft_values(6,4,40,0,1,0)

    else if ('RVV10' .EQ. TRIM(dftout) ) then
    ! Special case rVV10
       dft_defined = set_dft_values(1,4,13,4,3,0)
    else if ('RVV10-SCAN' .EQ. TRIM(dftout) ) then
      ! Special case rVV10+scan
         dft_defined = set_dft_values(0,0,0,0,3,5)
    else if ('B3LYP'.EQ. TRIM(dftout) ) then
    ! special case : B3LYP hybrid
       dft_defined = set_dft_values(7,12,9,7,0,0)

    else if ('B3LYP-V1R'.EQ. TRIM(dftout) ) then
    ! special case : B3LYP-VWN-1-RPA hybrid
       dft_defined = set_dft_values(7,13,9,7,0,0)

    else if ('X3LYP'.EQ. TRIM(dftout) ) then
    ! special case : X3LYP hybrid
       dft_defined = set_dft_values(9,14,28,13,0,0)

    ! special case : TPSS meta-GGA Exc
    else IF ('TPSS'.EQ. TRIM(dftout ) ) THEN
       dft_defined = set_dft_values(1,4,7,6,0,1)

    ! special case : M06L Meta GGA
    else if ( 'M06L' .EQ. TRIM(dftout) ) THEN
       dft_defined = set_dft_values(0,0,0,0,0,2)

    ! special case : TB09 meta-GGA Exc
    else IF ('TB09'.EQ. TRIM(dftout) ) THEN
       dft_defined = set_dft_values(0,0,0,0,0,3)

   ! special case : SCAN Meta GGA
    else if ( 'SCAN' .EQ. TRIM(dftout) ) THEN
       dft_defined = set_dft_values(0,0,0,0,0,5)

     ! special case : SCAN0
      else IF ('SCAN0'.EQ. TRIM(dftout ) ) THEN
        dft_defined = set_dft_values(0,0,0,0,0,6)

    ! special case : PZ/LDA + null meta-GGA
    else IF (('PZ+META'.EQ. TRIM(dftout)) .or. ('LDA+META'.EQ. TRIM(dftout)) ) THEN
       dft_defined = set_dft_values(1,1,0,0,0,4)

    ! special case : PBE + null meta-GGA
    else IF ('PBE+META'.EQ. TRIM(dftout) ) THEN
       dft_defined = set_dft_values(1,4,3,4,0,4)

      else if ('REV-VDW-DF2' .EQ. TRIM(dftout) ) then
        call errore('set_dft_from_name','obsolete XC label, use VDW-DF2-B86R',1)

     else if ('VDW-DF3' .EQ. TRIM(dftout) ) then
        call errore('set_dft_from_name','obsolete XC label, use VDW-DF-OBK8',1)

     else if ('VDW-DF4'.EQ.TRIM(dftout) .OR. 'OPTB86B-VDW'.EQ.TRIM(dftout) ) then
        call errore('set_dft_from_name','obsolete XC label, use VDW-DF-OB86',1)

     else if ('VDW-DF-X' .EQ. TRIM(dftout) ) then
     ! Special case vdW-DF-X
        call errore('set_dft_from_name','functional not yet implemented',1)

     else if ('VDW-DF-Y' .EQ. TRIM(dftout) ) then
     ! Special case vdW-DF-Y
        call errore('set_dft_from_name','functional not yet implemented',1)

     else if ('VDW-DF-Z' .EQ. TRIM(dftout) ) then
     ! Special case vdW-DF-Z
        call errore('set_dft_from_name','functional not yet implemented',1)

     ELSE IF ( 'INDEX:' ==  dftout(1:6)) THEN
     ! Special case for old RRKJ format, containing indices instead of label
        READ( dftout(7:18), '(6i2)') iexch, icorr, igcx, igcc, inlc, imeta
        dft_defined = set_dft_values(iexch, icorr, igcx, igcc, inlc, imeta)
        dftout = exc (iexch) //'-'//corr (icorr) //'-'//gradx (igcx) //'-' &
             &//gradc (igcc) //'-'// nonlocc(inlc)
    END IF

    !
    ! ----------------------------------------------------------------
    ! If the DFT was not yet defined, check every part of the string
    ! ----------------------------------------------------------------
    !
    if (.not. dft_defined) then

      iexch = matching (dftout,  nxc, exc)
      icorr = matching (dftout,  ncc, corr)
      igcx  = matching (dftout, ngcx,gradx)
      igcc  = matching (dftout, ngcc,gradc)
      imeta = matching (dftout,nmeta, meta)
      inlc  = matching (dftout, ncnl, nonlocc)
      
    endif
     
    ! ----------------------------------------------------------------
    ! Last check
    ! No more defaults, the code exits if the dft is not defined
    ! ----------------------------------------------------------------

    ! Back compatibility - TO BE REMOVED

    if (igcx == 14) igcx = 3 ! PBE -> PBX
    if (igcc ==  9) igcc = 4 ! PBE -> PBC

    if (igcx == 6) &
         call infomsg('set_dft_from_name','OPTX untested! please test')

    ! check for unrecognized labels

    if ( iexch<=0.and.icorr<=0.and.igcx<=0.and.igcc<= 0.and.imeta<=0 ) then
        if ( inlc <= 0 .and. trim(dftout) /= 'NOX-NOC') then
           call errore('set_dft_from_name',trim(dftout)//': unrecognized dft',1)
        else
           ! if inlc is the only nonzero index the label is likely wrong
           call errore('set_dft_from_name',trim(dftout)//': strange dft, please check',inlc)
        endif
    endif
    !
    ! Fill variables and exit
    !
    dft = dftout

    !dft_longname = exc (iexch) //'-'//corr (icorr) //'-'//gradx (igcx) //'-' &
    !     &//gradc (igcc) //'-'// nonlocc(inlc)

    call set_auxiliary_flags
    !
    ! check dft has not been previously set differently
    !
    if (save_iexch .ne. notset .and. save_iexch .ne. iexch) then
       write (stdout,*) iexch, save_iexch
       call errore('set_dft_from_name',' conflicting values for iexch',1)
    end if
    if (save_icorr .ne. notset .and. save_icorr .ne. icorr) then
       write (stdout,*) icorr, save_icorr
       call errore('set_dft_from_name',' conflicting values for icorr',1)
    end if
    if (save_igcx .ne. notset .and. save_igcx .ne. igcx) then
       write (stdout,*) igcx, save_igcx
       call errore('set_dft_from_name',' conflicting values for igcx',1)
    end if
    if (save_igcc .ne. notset .and. save_igcc .ne. igcc) then
       write (stdout,*) igcc, save_igcc
       call errore('set_dft_from_name',' conflicting values for igcc',1)
    end if
    if (save_meta .ne. notset .and. save_meta .ne. imeta) then
       write (stdout,*) inlc, save_meta
       call errore('set_dft_from_name',' conflicting values for imeta',1)
    end if
    if (save_inlc .ne. notset .and. save_inlc .ne. inlc) then
       write (stdout,*) inlc, save_inlc
       call errore('set_dft_from_name',' conflicting values for inlc',1)
    end if

    return
  end subroutine set_dft_from_name
  !
  integer function matching(dft, n, name)
  implicit none
  integer, intent(in):: n
  character(len=*), intent(in):: name(0:n)
  character(len=*), intent(in):: dft
  logical, external :: matches
  integer :: i

  matching = notset
  do i = n, 0, -1
     if (matches (name (i), trim(dft)) ) then
        if ( matching == notset ) then
           ! write(*, '("matches",i2,2X,A,2X,A)') i, name(i), trim(dft)
           matching = i
        else
           write(*, '(2(2X,i2,2X,A))') i, trim(name(i)), &
                                   matching, trim(name(matching))
           call errore ('set_dft', 'two conflicting matching values', 1)
        end if
     endif
  end do
  if (matching == notset) matching = 0
  !
  end function matching
  !
  !-----------------------------------------------------------------------
  subroutine set_auxiliary_flags
    !-----------------------------------------------------------------------
    ! set logical flags describing the complexity of the xc functional
    ! define the fraction of exact exchange used by hybrid fuctionals
    !
    isnonlocc = (inlc > 0)
    ismeta    = (imeta > 0)
    isgradient= (igcx > 0) .or. (igcc > 0) .or. ismeta .or. isnonlocc
    islda     = (iexch> 0) .and. (icorr > 0) .and. .not. isgradient
    ! PBE0/DF0
    IF ( iexch==6 .or. igcx ==8 ) exx_fraction = 0.25_DP
    ! CX0P
    IF ( iexch==6 .AND. igcx ==31 ) exx_fraction = 0.20_DP
    ! B86BPBEX
    IF ( iexch==6 .and. igcx ==41 ) exx_fraction = 0.25_DP
    ! BHANDHLYP
    IF ( iexch==6 .and. igcx ==42 ) exx_fraction = 0.5_DP
    ! HSE
    IF ( igcx ==12 ) THEN
       exx_fraction = 0.25_DP
       screening_parameter = 0.106_DP
    END IF
    ! gau-pbe
    IF ( igcx ==20 ) THEN
       exx_fraction = 0.24_DP
       gau_parameter = 0.150_DP
    END IF
    ! HF or OEP
    IF ( iexch==4 .or. iexch==5 ) exx_fraction = 1.0_DP
    ! B3LYP or B3LYP-VWN-1-RPA
    IF ( iexch == 7 ) exx_fraction = 0.2_DP
    ! X3LYP
    IF ( iexch == 9 ) exx_fraction = 0.218_DP
    !
    ishybrid = ( exx_fraction /= 0.0_DP )

    has_finite_size_correction = ( iexch==8 .or. icorr==10)

    return
  end subroutine set_auxiliary_flags
  !
  !-----------------------------------------------------------------------
  logical function set_dft_values(i1,i2,i3,i4,i5,i6)
    !-----------------------------------------------------------------------
    !
    implicit none
    integer :: i1,i2,i3,i4,i5,i6
    !
    iexch=i1
    icorr=i2
    igcx =i3
    igcc =i4
    inlc =i5
    imeta=i6
    !
    set_dft_values = .true.
    !
    return
    !
  end function set_dft_values
  !
  !
  !--------------------------------------------------------------------------
  FUNCTION qe_to_libxc_index( index_qe, term )
    !-----------------------------------------------------------------------
    !! PROVISIONAL: converts q-e indexes of functionals into libxc ones,
    !! when possible.
    !! In the next commit this will be done directly in 'funct.f90' and
    !! including all the cases available.
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: index_qe
    !! index of q-e functional term
    CHARACTER(LEN=*), INTENT(IN) :: term
    !! which term to consider
    INTEGER :: qe_to_libxc_index
    !! index converted in libxc notation
    !
    ! ... exchange term
    IF (term .EQ. 'exch_LDA') THEN
      !
      qe_or_libxc(1) = 1
      !
      SELECT CASE( index_qe )
      CASE( 1 )
        qe_to_libxc_index = 1
      CASE DEFAULT
        qe_to_libxc_index = index_qe
        qe_or_libxc(1) = 0
      END SELECT
      !
      ! ... correlation term
    ELSEIF (term .EQ. 'corr_LDA') THEN
      !
      qe_or_libxc(2) = 1
      !
      SELECT CASE( index_qe )
      CASE( 1 ) !pz
         qe_to_libxc_index = 9
      CASE( 5 ) !wigner
         qe_to_libxc_index = 2
      CASE( 2 ) !vwn
         qe_to_libxc_index = 7
      CASE( 4 ) !pw
         qe_to_libxc_index = 12
      CASE DEFAULT
         qe_to_libxc_index = index_qe
         qe_or_libxc(2) = 0
      END SELECT      
      !
    ELSE
      !
      CALL errore( 'qe_to_libxc_index', 'ERROR: option not available', 3 )
      !
    ENDIF
    !
    RETURN
    !
  END FUNCTION qe_to_libxc_index
  !
  !
  !-----------------------------------------------------------------------
  subroutine enforce_input_dft (dft_, nomsg)
    !
    ! translates a string containing the exchange-correlation name
    ! into internal indices and force any subsequent call to set_dft_from_name
    ! to return without changing them
    !
    implicit none
    character(len=*), intent(in) :: dft_
    logical, intent(in), optional :: nomsg

     call set_dft_from_name (dft_)
     if (dft == 'not set') call errore('enforce_input_dft','cannot fix unset dft',1)
     discard_input_dft = .true.

     if ( present (nomsg) ) return

     write (stdout,'(/,5x,a)') "IMPORTANT: XC functional enforced from input :"
     call write_dft_name
     write (stdout,'(5x,a)') "Any further DFT definition will be discarded"
     write (stdout,'(5x,a/)') "Please, verify this is what you really want"

     return
  end subroutine enforce_input_dft

  !-----------------------------------------------------------------------
  subroutine enforce_dft_exxrpa ( )
    !
    implicit none
    !
    !character(len=*), intent(in) :: dft_
    !logical, intent(in), optional :: nomsg

    iexch = 0; icorr = 0; igcx = 0; igcc = 0
    exx_fraction = 1.0_DP
    ishybrid = ( exx_fraction /= 0.0_DP )

    write (stdout,'(/,5x,a)') "XC functional enforced to be EXXRPA"
    call write_dft_name
    write (stdout,'(5x,a)') "!!! Any further DFT definition will be discarded"
    write (stdout,'(5x,a/)') "!!! Please, verify this is what you really want !"

    return
  end subroutine enforce_dft_exxrpa

  !-----------------------------------------------------------------------
  subroutine init_dft_exxrpa ( )
    !
    implicit none
    !
    exx_fraction = 1.0_DP
    ishybrid = ( exx_fraction /= 0.0_DP )

    write (stdout,'(/,5x,a)') "Only exx_fraction is set to 1.d0"
    write (stdout,'(5x,a)') "XC functional still not changed"
    call write_dft_name

    return
  end subroutine init_dft_exxrpa

  !-----------------------------------------------------------------------
  subroutine start_exx
     if (.not. ishybrid) &
        call errore('start_exx','dft is not hybrid, wrong call',1)
     exx_started = .true.
  end subroutine start_exx
  !-----------------------------------------------------------------------
  subroutine stop_exx
     if (.not. ishybrid) &
        call errore('stop_exx','dft is not hybrid, wrong call',1)
     exx_started = .false.
  end subroutine stop_exx
  !-----------------------------------------------------------------------
  subroutine dft_force_hybrid(request)
    LOGICAL,OPTIONAL,INTENT(inout) :: request
    LOGICAL :: aux
    IF(present(request)) THEN
      aux = ishybrid
      ishybrid = request
      request = aux
    ELSE
      ishybrid= .true.
    ENDIF
 end subroutine dft_force_hybrid
 !-----------------------------------------------------------------------
  function exx_is_active ()
     logical exx_is_active
     exx_is_active = exx_started
  end function exx_is_active
  !-----------------------------------------------------------------------
  subroutine set_exx_fraction (exxf_)
     implicit none
     real(DP):: exxf_
     exx_fraction = exxf_
     write (stdout,'(5x,a,f6.2)') 'EXX fraction changed: ',exx_fraction
  end subroutine set_exx_fraction
  !---------------------------------------------------------------------
  subroutine set_screening_parameter (scrparm_)
     implicit none
     real(DP):: scrparm_
     screening_parameter = scrparm_
     write (stdout,'(5x,a,f12.7)') 'EXX Screening parameter changed: ', &
          & screening_parameter
  end subroutine set_screening_parameter
  !----------------------------------------------------------------------
  function get_screening_parameter ()
     real(DP):: get_screening_parameter
     get_screening_parameter = screening_parameter
     return
  end function get_screening_parameter
  !---------------------------------------------------------------------
  subroutine set_gau_parameter (gauparm_)
     implicit none
     real(DP):: gauparm_
     gau_parameter = gauparm_
     write (stdout,'(5x,a,f12.7)') 'EXX Gau parameter changed: ', &
          & gau_parameter
  end subroutine set_gau_parameter
  !----------------------------------------------------------------------
  function get_gau_parameter ()
     real(DP):: get_gau_parameter
     get_gau_parameter = gau_parameter
     return
  end function get_gau_parameter
  !-----------------------------------------------------------------------
  function get_iexch ()
     integer get_iexch
     get_iexch = iexch
     return
  end function get_iexch
  !-----------------------------------------------------------------------
  function get_icorr ()
     integer get_icorr
     get_icorr = icorr
     return
  end function get_icorr
  !-----------------------------------------------------------------------
  function get_igcx ()
     integer get_igcx
     get_igcx = igcx
     return
  end function get_igcx
  !-----------------------------------------------------------------------
  function get_igcc ()
     integer get_igcc
     get_igcc = igcc
     return
  end function get_igcc
  !-----------------------------------------------------------------------
  function get_meta ()
     integer get_meta
     get_meta = imeta
     return
  end function get_meta
  !-----------------------------------------------------------------------
  function get_inlc ()
     integer get_inlc
     get_inlc = inlc
     return
  end function get_inlc
  !-----------------------------------------------------------------------
  function get_nonlocc_name ()
     character(10) get_nonlocc_name
     get_nonlocc_name = TRIM(nonlocc(inlc))
     return
  end function get_nonlocc_name
 !-----------------------------------------------------------------------
  function dft_is_nonlocc ()
    logical :: dft_is_nonlocc
    dft_is_nonlocc = isnonlocc
    return
  end function dft_is_nonlocc
  !-----------------------------------------------------------------------
  function get_exx_fraction ()
     real(DP):: get_exx_fraction
     get_exx_fraction = exx_fraction
     return
  end function get_exx_fraction
  !-----------------------------------------------------------------------
  function get_dft_name ()
     character (len=25) :: get_dft_name
     get_dft_name = dft
     return
  end function get_dft_name
  !-----------------------------------------------------------------------
  function dft_is_gradient ()
     logical :: dft_is_gradient
     dft_is_gradient = isgradient
     return
  end function dft_is_gradient
  !-----------------------------------------------------------------------
  function dft_is_meta ()
     logical :: dft_is_meta
     dft_is_meta = ismeta
     return
  end function dft_is_meta
  !-----------------------------------------------------------------------
  function dft_is_hybrid ()
     logical :: dft_is_hybrid
     dft_is_hybrid = ishybrid
     return
  end function dft_is_hybrid
  !-----------------------------------------------------------------------
  function igcc_is_lyp ()
     logical :: igcc_is_lyp
     igcc_is_lyp = (get_igcc()==3 .or. get_igcc()==7 .OR. get_igcc()==13)
     return
  end function igcc_is_lyp
  !-----------------------------------------------------------------------
  function dft_has_finite_size_correction ()
     logical :: dft_has_finite_size_correction
     dft_has_finite_size_correction = has_finite_size_correction
     return
  end function  dft_has_finite_size_correction
  !-----------------------------------------------------------------------
  subroutine set_finite_size_volume(volume)
     real, intent (IN) :: volume
     if (.not. has_finite_size_correction) &
         call errore('set_finite_size_volume', &
                      'dft w/o finite_size_correction, wrong call',1)
     if (volume <= 0.d0) &
         call errore('set_finite_size_volume', &
                     'volume is not positive, check omega and/or nk1,nk2,nk3',1)
     finite_size_cell_volume = volume
     finite_size_cell_volume_set = .TRUE.
  end subroutine set_finite_size_volume
  !-----------------------------------------------------------------------
  !
  subroutine get_finite_size_cell_volume(is_present, volume)
     logical, intent(out) :: is_present
     real(dp), intent(out) :: volume
     is_present = finite_size_cell_volume_set
     volume = -1.d0
     if (is_present) volume = finite_size_cell_volume
  end subroutine get_finite_size_cell_volume
  !
  !-----------------------------------------------------------------------
  subroutine set_dft_from_indices(iexch_,icorr_,igcx_,igcc_, inlc_)
     integer :: iexch_, icorr_, igcx_, igcc_, inlc_
     if ( discard_input_dft ) return
     if (iexch == notset) iexch = iexch_
     if (iexch /= iexch_) then
        write (stdout,*) iexch, iexch_
        call errore('set_dft',' conflicting values for iexch',1)
     end if
     if (icorr == notset) icorr = icorr_
     if (icorr /= icorr_) then
        write (stdout,*) icorr, icorr_
        call errore('set_dft',' conflicting values for icorr',1)
     end if
     if (igcx  == notset) igcx = igcx_
     if (igcx /= igcx_) then
        write (stdout,*) igcx, igcx_
        call errore('set_dft',' conflicting values for igcx',1)
     end if
     if (igcc  == notset) igcc = igcc_
     if (igcc /= igcc_) then
        write (stdout,*) igcc, igcc_
        call errore('set_dft',' conflicting values for igcc',1)
     end if
     if (inlc  == notset) inlc = inlc_
     if (inlc /= inlc_) then
        write (stdout,*) inlc, inlc_
        call errore('set_dft',' conflicting values for inlc',1)
     end if
     dft = exc (iexch) //'-'//corr (icorr) //'-'//gradx (igcx) //'-' &
           &//gradc (igcc)//'-'//nonlocc (inlc)
     ! WRITE( stdout,'(a)') dft
     call set_auxiliary_flags
     return
  end subroutine set_dft_from_indices
  !-------------------------------------------------------------------------------
  function get_dft_short ( )
  !---------------------------------------------------------------------
  !
  character (len=12) :: get_dft_short
  character (len=12) :: shortname
  !
  shortname = 'no shortname'
  if (iexch==1.and.igcx==0.and.igcc==0) then
      shortname = TRIM(corr(icorr))
   else if ( iexch==4.and.icorr==0.and.igcx==0.and.igcc==0) then
      shortname = 'OEP'
   else if (iexch==1.and.icorr==11.and.igcx==0.and.igcc==0) then
      shortname = 'VWN-RPA'
  else if (iexch==1.and.icorr==3.and.igcx==1.and.igcc==3) then
     shortname = 'BLYP'
  else if (iexch==1.and.icorr==1.and.igcx==1.and.igcc==0) then
     shortname = 'B88'
  else if (iexch==1.and.icorr==1.and.igcx==1.and.igcc==1) then
     shortname = 'BP'
  else if (iexch==1.and.icorr==4.and.igcx==2.and.igcc==2) then
     shortname = 'PW91'
  else if (iexch==1.and.icorr==4.and.igcx==3.and.igcc==4) then
     shortname = 'PBE'
  else if (iexch==6.and.icorr==4.and.igcx==8.and.igcc==4) then
     shortname = 'PBE0'
  else if (iexch==6.and.icorr==4.and.igcx==41.and.igcc==4) then
     shortname = 'B86BPBEX'
  else if (iexch==6.and.icorr==4.and.igcx==42.and.igcc==3) then
     shortname = 'BHANDHLYP'
  else if (iexch==1.and.icorr==4.and.igcx==4.and.igcc==4) then
     shortname = 'revPBE'
  else if (iexch==1.and.icorr==4.and.igcx==10.and.igcc==8) then
     shortname = 'PBESOL'
  else if (iexch==1.and.icorr==4.and.igcx==19.and.igcc==12) then
     shortname = 'Q2D'
  else if (iexch==1.and.icorr==4.and.igcx==12.and.igcc==4) then
     shortname = 'HSE'
  else if (iexch==1.and.icorr==4.and.igcx==20.and.igcc==4) then
     shortname = 'GAUPBE'
  else if (iexch==1.and.icorr==4.and.igcx==21.and.igcc==4) then
      shortname = 'PW86PBE'
  else if (iexch==1.and.icorr==4.and.igcx==22.and.igcc==4) then
      shortname = 'B86BPBE'
  else if (iexch==1.and.icorr==4.and.igcx==11.and.igcc==4) then
     shortname = 'WC'
  else if (iexch==7.and.icorr==12.and.igcx==9.and. igcc==7) then
     shortname = 'B3LYP'
  else if (iexch==7.and.icorr==13.and.igcx==9.and. igcc==7) then
     shortname = 'B3LYP-V1R'
  else if (iexch==9.and.icorr==14.and.igcx==28.and. igcc==13) then
     shortname = 'X3LYP'
  else if (iexch==0.and.icorr==3.and.igcx==6.and.igcc==3) then
     shortname = 'OLYP'
  else if (iexch==1.and.icorr==4.and.igcx==17.and.igcc==4) then
     shortname = 'SOGGA'
  else if (iexch==1.and.icorr==4.and.igcx==23.and.igcc==1) then
      shortname = 'OPTBK88'
  else if (iexch==1.and.icorr==4.and.igcx==24.and.igcc==1) then
      shortname = 'OPTB86B'
  else if (iexch==1.and.icorr==4.and.igcx==25.and.igcc==0) then
      shortname = 'EV93'
  else if (iexch==5.and.icorr==0.and.igcx==0.and.igcc==0) then
      shortname = 'HF'
  end if

  if (imeta == 1 ) then
     shortname = 'TPSS'
  else if (imeta == 2) then
     shortname = 'M06L'
  else if (imeta == 3) then
     shortname = 'TB09'
  else if (imeta == 4) then
    if ( iexch == 1 .and. icorr == 1) then
      shortname = 'PZ+META'
    else if (iexch==1.and.icorr==4.and.igcx==3.and.igcc==4) then
      shortname = 'PBE+META'
    end if
  else if (imeta == 5 ) then
    shortname = 'SCAN'
  else if (imeta == 6 ) then
    shortname = 'SCAN0'
  end if

  if ( inlc==1 ) then
     if (iexch==1.and.icorr==4.and.igcx==4.and.igcc==0) then
        shortname = 'VDW-DF'
     else if (iexch==1.and.icorr==4.and.igcx==27.and.igcc==0) then
        shortname = 'VDW-DF-CX'
     else if (iexch==6.and.icorr==4.and.igcx==29.and.igcc==0) then
        shortname = 'VDW-DF-CX0'
     else if (iexch==6.and.icorr==4.and.igcx==31.and.igcc==0) then
        shortname = 'VDW-DF-CX0P'
     else if (iexch==1.and.icorr==4.and.igcx==16.and.igcc==0) then
        shortname = 'VDW-DF-C09'
     else if (iexch==1.and.icorr==4.and.igcx==24.and.igcc==0) then
        shortname = 'VDW-DF-OB86'
     else if (iexch==1.and.icorr==4.and.igcx==23.and.igcc==0) then
        shortname = 'VDW-DF-OBK8'
     else if (iexch==6.and.icorr==4.and.igcx==40.and.igcc==0) then
        shortname = 'VDW-DF-C090'
     end if
  else if ( inlc==2 ) then
     if (iexch==1.and.icorr==4.and.igcx==13.and.igcc==0) then
        shortname = 'VDW-DF2'
     else if (iexch==1.and.icorr==4.and.igcx==16.and.igcc==0) then
        shortname = 'VDW-DF2-C09'
     else if (iexch==1.and.icorr==4.and.igcx==26.and.igcc==0) then
        shortname = 'VDW-DF2-B86R'
     else if (iexch==6.and.icorr==4.and.igcx==30.and.igcc==0) then
        shortname = 'VDW-DF2-0'
     else if (iexch==6.and.icorr==4.and.igcx==38.and.igcc==0) then
        shortname = 'VDW-DF2-BR0'
     end if
  else if ( inlc==3) then
     shortname = 'RVV10'
  end if
  !
  get_dft_short = shortname
  !
  end function get_dft_short
  !
  !---------------------------------------------------------------------
  function get_dft_long( )
  !---------------------------------------------------------------------
  !
  character (len=25) :: get_dft_long
  character (len=25) :: longname
  !
  write(longname,'(4a5)') exc(iexch),corr(icorr),gradx(igcx),gradc(igcc)
  if ( imeta > 0 ) then
     longname = longname(1:20)//trim(meta(imeta))
  else if ( inlc > 0 ) then
     longname = longname(1:20)//trim(nonlocc(inlc))
  end if
  get_dft_long = longname

end function get_dft_long

!-----------------------------------------------------------------------
subroutine write_dft_name
!-----------------------------------------------------------------------
   WRITE( stdout, '(5X,"Exchange-correlation      = ",A, &
        &  " (",I2,3I3,2I2,")")') TRIM( dft ), iexch,icorr,igcx,igcc,inlc,imeta
   IF ( get_exx_fraction() > 0.0_dp ) WRITE( stdout, &
        '(5X,"EXX-fraction              =",F12.2)') get_exx_fraction()
   return
end subroutine write_dft_name
!
!
!-----------------------------------------------------------------------
!------- NONLOCAL CORRECTIONS DRIVERS ----------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
subroutine nlc (rho_valence, rho_core, nspin, enl, vnl, v)
  !-----------------------------------------------------------------------
  !     non local correction for the correlation
  !
  !     input:  rho_valence, rho_core
  !     definition:  E_nl = \int E_nl(rho',grho',rho'',grho'',|r'-r''|) dr
  !     output: enl = E_nl
  !             vnl= D(E_x)/D(rho)
  !             v  = Correction to the potential
  !

  USE vdW_DF, ONLY: xc_vdW_DF, xc_vdW_DF_spin, vdw_type
  USE rVV10,  ONLY: xc_rVV10

  implicit none

  REAL(DP), INTENT(IN) :: rho_valence(:,:), rho_core(:)
  INTEGER, INTENT(IN)  :: nspin
  REAL(DP), INTENT(INOUT) :: v(:,:)
  REAL(DP), INTENT(INOUT) :: enl, vnl

  if ( inlc == 1 .or. inlc == 2 .or. inlc == 4 .or. inlc == 5 .or. inlc == 6 ) then

     vdw_type = inlc
     if( nspin == 1 ) then
        call xc_vdW_DF      (rho_valence, rho_core, enl, vnl, v)
     else if( nspin == 2 ) then
        call xc_vdW_DF_spin (rho_valence, rho_core, enl, vnl, v)
     else
        call errore ('nlc','vdW-DF not available for noncollinear spin case',1)
     end if

  elseif (inlc == 3) then
      if(imeta == 0) then
        call xc_rVV10 (rho_valence(:,1), rho_core, nspin, enl, vnl, v)
      else
        call xc_rVV10 (rho_valence(:,1), rho_core, nspin, enl, vnl, v, 15.7_dp)
      endif
  else
     enl = 0.0_DP
     vnl = 0.0_DP
     v = 0.0_DP
  endif
  !
  return
end subroutine nlc

!
!-----------------------------------------------------------------------
!------- META CORRECTIONS DRIVERS ----------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
subroutine tau_xc (rho, grho, tau, ex, ec, v1x, v2x, v3x, v1c, v2c, v3c)
  !-----------------------------------------------------------------------
  !     gradient corrections for exchange and correlation - Hartree a.u.
  !     See comments at the beginning of module for implemented cases
  !
  !     input:  rho, grho=|\nabla rho|^2
  !
  !     definition:  E_x = \int e_x(rho,grho) dr
  !
  !     output: sx = e_x(rho,grho) = grad corr
  !             v1x= D(E_x)/D(rho)
  !             v2x= D(E_x)/D( D rho/D r_alpha ) / |\nabla rho|
  !             v3x= D(E_x)/D(tau)
  !
  !             sc, v1c, v2c as above for correlation
  !
  implicit none

  real(DP) :: rho, grho, tau, ex, ec, v1x, v2x, v3x, v1c, v2c, v3c

  !_________________________________________________________________________

  if     (imeta == 1) then
     call tpsscxc (rho, grho, tau, ex, ec, v1x, v2x, v3x, v1c, v2c, v3c)
  elseif (imeta == 2) then
     call   m06lxc (rho, grho, tau, ex, ec, v1x, v2x, v3x, v1c, v2c, v3c)
  elseif (imeta == 3) then
     call  tb09cxc (rho, grho, tau, ex, ec, v1x, v2x, v3x, v1c, v2c, v3c)
  elseif (imeta == 4) then
     ! do nothing
  elseif (imeta == 5) then
     call  SCANcxc (rho, grho, tau, ex, ec, v1x, v2x, v3x, v1c, v2c, v3c)
  else
    call errore('tau_xc','wrong igcx and/or igcc',1)
  end if

  return

end subroutine tau_xc

subroutine tau_xc_array (nnr, rho, grho, tau, ex, ec, v1x, v2x, v3x, v1c, v2c, v3c)
  ! HK/MCA : the xc_func_init is slow and is called too many times
  ! HK/MCA : we modify this subroutine so that the overhead could be minimized
  !-----------------------------------------------------------------------
  !     gradient corrections for exchange and correlation - Hartree a.u.
  !     See comments at the beginning of module for implemented cases
  !
  !     input:  rho, grho=|\nabla rho|^2
  !
  !     definition:  E_x = \int e_x(rho,grho) dr
  !
  !     output: sx = e_x(rho,grho) = grad corr
  !             v1x= D(E_x)/D(rho)
  !             v2x= D(E_x)/D( D rho/D r_alpha ) / |\nabla rho|
  !             v3x= D(E_x)/D(tau)
  !
  !             sc, v1c, v2c as above for correlation
  !
  implicit none

  integer, intent(in) :: nnr
  real(DP) :: rho(nnr), grho(nnr), tau(nnr), ex(nnr), ec(nnr)
  real(DP) :: v1x(nnr), v2x(nnr), v3x(nnr), v1c(nnr), v2c(nnr), v3c(nnr)
  !_________________________________________________________________________

  if (imeta == 5) then
     call  scancxc_array (nnr, rho, grho, tau, ex, ec, v1x, v2x, v3x, v1c, v2c, v3c)
  elseif (imeta == 6 ) then ! HK/MCA: SCAN0
     call  scancxc_array (nnr, rho, grho, tau, ex, ec, v1x, v2x, v3x, v1c, v2c, v3c)
     if (exx_started) then
        ex  = (1.0_DP - exx_fraction) * ex
        v1x = (1.0_DP - exx_fraction) * v1x
        v2x = (1.0_DP - exx_fraction) * v2x
        v3x = (1.0_DP - exx_fraction) * v3x
     end if
  else
     call errore('v_xc_meta_array','(CP only) array mode only works for SCAN',1)
  end if

  return

end subroutine tau_xc_array
!
!
!-----------------------------------------------------------------------
subroutine tau_xc_spin (rhoup, rhodw, grhoup, grhodw, tauup, taudw, ex, ec,   &
           &            v1xup, v1xdw, v2xup, v2xdw, v3xup, v3xdw, v1cup, v1cdw,&
           &            v2cup, v2cdw, v3cup, v3cdw)

!-----------------------------------------------------------------------
  !
  !

  implicit none

  real(dp), intent(in)                :: rhoup, rhodw, tauup, taudw
  real(dp), dimension (3), intent(in) :: grhoup, grhodw

  real(dp), intent(out)               :: ex, ec, v1xup, v1xdw, v2xup, v2xdw, v3xup, v3xdw,  &
                                      &  v1cup, v1cdw, v3cup, v3cdw
  real(dp), dimension(3), intent(out) :: v2cup, v2cdw

  !
  !  Local variables
  !
  integer                 :: ipol
  real(dp)                :: rh, zeta, atau, grhoup2, grhodw2
  real(dp), parameter     :: epsr=1.0d-08, zero=0._dp
  !
  !_____________________________

  grhoup2 = zero
  grhodw2 = zero

  v2cup         = zero
  v2cdw         = zero

  ! FIXME: for SCAN, this will be calculated later
  if (imeta /= 4) then

      do ipol=1,3
        grhoup2 = grhoup2 + grhoup(ipol)**2
        grhodw2 = grhodw2 + grhodw(ipol)**2
      end do

  end if

  if (imeta == 1) then

     call tpsscx_spin(rhoup, rhodw, grhoup2, grhodw2, tauup,   &
              &  taudw, ex, v1xup,v1xdw,v2xup,v2xdw,v3xup,v3xdw)

     rh   =  rhoup + rhodw

     zeta = (rhoup - rhodw) / rh
     atau =  tauup + taudw    ! KE-density in Hartree

     call tpsscc_spin(rh,zeta,grhoup,grhodw, atau,ec,              &
     &                v1cup,v1cdw,v2cup,v2cdw,v3cup, v3cdw)


  elseif (imeta == 2) then

     call   m06lxc_spin (rhoup, rhodw, grhoup2, grhodw2, tauup, taudw,      &
            &            ex, ec, v1xup, v1xdw, v2xup, v2xdw, v3xup, v3xdw,  &
            &            v1cup, v1cdw, v2cup(1), v2cdw(1), v3cup, v3cdw)

  elseif (imeta == 5) then

            ! FIXME: not the most efficient use of libxc

            call scanxc_spin(rhoup, rhodw, grhoup, grhodw, tauup, taudw,  &
                     &  ex, v1xup,v1xdw,v2xup,v2xdw,v3xup,v3xdw,          &
                     &  ec, v1cup,v1cdw,v2cup,v2cdw,v3cup,v3cdw )

  else

    call errore('tau_xc_spin','This case not implemented',imeta)

  end if

end subroutine tau_xc_spin

subroutine tau_xc_array_spin (nnr, rho, grho, tau, ex, ec, v1x, v2x, v3x, &
  & v1c, v2c, v3c)
! HK/MCA : the xc_func_init (LIBXC) is slow and is called too many times
! HK/MCA : we modify this subroutine so that the overhead could be minimized
!-----------------------------------------------------------------------
!     gradient corrections for exchange and correlation - Hartree a.u.
!     See comments at the beginning of module for implemented cases
!
!     input:  rho,rho, grho=\nabla rho
!
!     definition:  E_x = \int e_x(rho,grho) dr
!
!     output: sx = e_x(rho,grho) = grad corr
!             v1x= D(E_x)/D(rho)
!             v2x= D(E_x)/D( D rho/D r_alpha ) / |\nabla rho|
!             v3x= D(E_x)/D(tau)
!
!             sc, v1cup, v2cup as above for correlation
!
implicit none

integer, intent(in) :: nnr
real(DP) :: rho(nnr,2), grho(3,nnr,2), tau(nnr,2), ex(nnr), ec(nnr)
real(DP) :: v1x(nnr,2), v2x(nnr,3), v3x(nnr,2), v1c(nnr,2), v2c(nnr,3), v3c(nnr,2)

!Local variables

integer  :: ipol, k, is
real(DP) :: grho2(3,nnr)
!MCA: Libxc format
real(DP) :: rho_(2,nnr), tau_(2,nnr)
real(DP) :: v1x_(2,nnr), v2x_(3,nnr), v3x_(2,nnr), v1c_(2,nnr), v2c_(3,nnr), v3c_(2,nnr)

!_________________________________________________________________________

grho2 = 0.0

!MCA/HK: contracted gradient of density, same format as in libxc
do k=1,nnr

do ipol=1,3
grho2(1,k) = grho2(1,k) + grho(ipol,k,1)**2
grho2(2,k) = grho2(2,k) + grho(ipol,k,1) * grho(ipol,k,2)
grho2(3,k) = grho2(3,k) + grho(ipol,k,2)**2
end do

!MCA: transforming to libxc format (DIRTY HACK)
do is=1,2
rho_(is,k) = rho(k,is)
tau_(is,k) = tau(k,is)
enddo

end do

if (imeta == 5) then

!MCA/HK: using the arrays in libxc format
call  scancxc_array_spin (nnr, rho_, grho2, tau_, ex, ec, &
&                   v1x_, v2x_, v3x_,  &
&                   v1c_, v2c_, v3c_ )

do k=1,nnr

!MCA: from libxc to QE format (DIRTY HACK)
do is=1,2
v1x(k,is) = v1x_(is,k)
v2x(k,is) = v2x_(is,k) !MCA/HK: v2x(:,2) contains the cross terms
v3x(k,is) = v3x_(is,k)
v1c(k,is) = v1c_(is,k)
v2c(k,is) = v2c_(is,k) !MCA/HK: same as v2x
v3c(k,is) = v3c_(is,k)
enddo

v2c(k,3) = v2c_(3,k)
v2x(k,3) = v2x_(3,k)

end do

elseif (imeta == 6 ) then ! HK/MCA: SCAN0
call  scancxc_array (nnr, rho, grho, tau, ex, ec, v1x, v2x, v3x, v1c, v2c, v3c)
if (exx_started) then
ex  = (1.0_DP - exx_fraction) * ex
v1x = (1.0_DP - exx_fraction) * v1x
v2x = (1.0_DP - exx_fraction) * v2x
v3x = (1.0_DP - exx_fraction) * v3x
end if
else
call errore('v_xc_meta_array','(CP only) array mode only works for SCAN',1)
end if

return

end subroutine tau_xc_array_spin
!
!-----------------------------------------------------------------------
SUBROUTINE init_lda_xc()
   !-------------------------------------------------------------------
   !! Gets from inside parameters needed to initialize lda xc-drivers.
   !
   USE kinds,          ONLY: DP
   USE xc_lda_lsda,    ONLY: libxc_switches, iexch_l, icorr_l, &
                             exx_started_l, is_there_finite_size_corr, &
                             exx_fraction_l, finite_size_cell_volume_l
   !
   IMPLICIT NONE
   !
   ! use qe or libxc for the different terms (0: qe, 1: libxc)
   libxc_switches(:) = qe_or_libxc(:)
   !
   ! exchange-correlation indexes
   iexch_l = get_iexch()
   IF (libxc_switches(1)==1) iexch_l = qe_to_libxc_index( iexch, 'exch_LDA' )
   icorr_l = get_icorr()
   IF (libxc_switches(2)==1) icorr_l = qe_to_libxc_index( icorr, 'corr_LDA' )
   !
   IF (iexch_l==-1 .OR. icorr_l==-1) CALL errore( 'init_lda_xc', 'Functional &
                                             & indexes not well defined', 1 )
   !
   ! hybrid exchange vars
   exx_started_l  = exx_is_active()
   exx_fraction_l = 0._DP
   IF ( exx_started ) exx_fraction_l = get_exx_fraction()
   !
   ! finite size correction vars
   CALL get_finite_size_cell_volume( is_there_finite_size_corr, &
                                     finite_size_cell_volume_l )
   !
   RETURN
   !
END SUBROUTINE
!
!
SUBROUTINE init_gga_xc()
   !! Extracts the parameters needed to initialize GGA xc-drivers.
   !
   USE kinds,               ONLY: DP
   USE xc_gga,              ONLY: igcx_l, igcc_l, &
                                  exx_started_g, exx_fraction_g, &
                                  screening_parameter_l, gau_parameter_l
   !
   IMPLICIT NONE
   !
   ! exchange-correlation indexes
   igcx_l = get_igcx()
   igcc_l = get_igcc()
   !
   ! hybrid exchange vars
   exx_started_g  = exx_is_active()
   exx_fraction_g = 0._DP
   IF ( exx_started ) exx_fraction_g = get_exx_fraction()
   !
   screening_parameter_l = get_screening_parameter()
   gau_parameter_l = get_gau_parameter()
   !
   RETURN
   !
END SUBROUTINE
!
!
#if defined(__LIBXC)
  subroutine get_libxc_version
     implicit none
     interface
        subroutine xc_version(major, minor, micro) bind(c)
           use iso_c_binding
           integer(c_int) :: major, minor, micro
        end subroutine xc_version
     end interface
     call xc_version(libxc_major, libxc_minor, libxc_micro)
  end subroutine get_libxc_version
#endif


end module funct

