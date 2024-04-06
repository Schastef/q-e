module atomic_projectors
   ! module for atomic projectors
   USE kinds, ONLY: DP

   ! Atomic projectors for species
   TYPE atproj_type
      CHARACTER(len=3) :: atsym ! atomic symbol
      INTEGER :: ngrid ! number of grid points
      REAL(DP), ALLOCATABLE :: xgrid(:)
      REAL(DP), ALLOCATABLE :: rgrid(:) ! exp(xgrid(:))
      INTEGER :: nproj ! total number of projectors for this species
      INTEGER, ALLOCATABLE :: l(:) ! angular momentum of each wfc
      REAL(DP), ALLOCATABLE :: radial(:, :) ! magnitude at each point, ngrid x nproj
   END TYPE atproj_type

   ! atomic proj input variables, start with atom_proj*
   LOGICAL :: atom_proj
   CHARACTER(LEN=256) :: atom_proj_dir ! directory of external projectors
   LOGICAL :: atom_proj_ext ! switch for using external files instead of orbitals from UPF
   INTEGER, PARAMETER :: nexatproj_max = 2000 ! max allowed number of projectors to be excluded
   INTEGER :: atom_proj_exclude(nexatproj_max) ! index starts from 1
   LOGICAL :: atom_proj_ortho ! whether perform Lowdin orthonormalization
   LOGICAL :: atom_proj_sym ! whether perform symmetrization

   ! atomic proj internal variables, using *atproj*
   INTEGER :: nexatproj ! actual number of excluded projectors
   INTEGER :: natproj ! total number of projectors = n_proj + nexatproj
   LOGICAL, ALLOCATABLE :: atproj_excl(:) ! size = total num of projectors
   INTEGER :: iun_atproj
   TYPE(atproj_type), ALLOCATABLE :: atproj_types(:) ! all atom proj types

   REAL(DP), ALLOCATABLE :: tab_at(:, :, :)
   !! interpolation table for atomic projectors

CONTAINS

   SUBROUTINE skip_comments(file_unit)
      ! read a file_unit and skip lines starting with #
      !
      USE, INTRINSIC :: ISO_FORTRAN_ENV
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(in) :: file_unit
      CHARACTER(len=256) :: line
      INTEGER :: ret_code
      !
      DO
         READ (file_unit, '(A)', iostat=ret_code) line
         IF (ret_code == IOSTAT_END) EXIT
         IF (ret_code /= 0) THEN
            ! read error
            EXIT
         END IF
         IF (INDEX(ADJUSTL(line), "#") == 1) CYCLE
         EXIT
      END DO

      BACKSPACE file_unit

      RETURN
   END SUBROUTINE skip_comments

   SUBROUTINE allocate_atproj_type(typ, ngrid, nproj)
      !
      IMPLICIT NONE
      !
      TYPE(atproj_type), INTENT(INOUT) :: typ
      INTEGER, INTENT(IN) :: ngrid, nproj
      INTEGER :: ierr
      !
      typ%ngrid = ngrid
      typ%nproj = nproj
      ALLOCATE (typ%xgrid(ngrid), stat=ierr)
      IF (ierr /= 0) CALL errore('atomic_projectors', 'Error allocating typ%xgrid', 1)
      ALLOCATE (typ%rgrid(ngrid), stat=ierr)
      IF (ierr /= 0) CALL errore('atomic_projectors', 'Error allocating typ%rgrid', 1)
      ALLOCATE (typ%l(nproj), stat=ierr)
      IF (ierr /= 0) CALL errore('atomic_projectors', 'Error allocating typ%l', 1)
      ALLOCATE (typ%radial(ngrid, nproj), stat=ierr)
      IF (ierr /= 0) CALL errore('atomic_projectors', 'Error allocating typ%radial', 1)

      RETURN
   END SUBROUTINE allocate_atproj_type

   SUBROUTINE read_atomproj(types)
      !
      ! read data files for atom proj
      ! should be called only by one node
      !
      USE kinds, ONLY: dp
      USE io_global, ONLY: stdout, ionode
      USE ions_base, ONLY: nsp, atm
      !
      IMPLICIT NONE
      !
      TYPE(atproj_type), INTENT(INOUT) :: types(nsp)
      !
      INTEGER :: i, j, it
      LOGICAL :: file_exists
      CHARACTER(len=256) :: filename
      INTEGER :: ngrid, nproj

      if (.not. ionode) return

      DO it = 1, nsp
         filename = TRIM(atom_proj_dir)//'/'//TRIM(atm(it))//".dat"
         INQUIRE (FILE=TRIM(filename), EXIST=file_exists)
         IF (.NOT. file_exists) &
            CALL errore('atomic_projectors', 'file not exists: '//TRIM(filename), 1)
         OPEN (NEWUNIT=iun_atproj, file=TRIM(filename), form='formatted')
         CALL skip_comments(iun_atproj)

         READ (iun_atproj, *) ngrid, nproj
         WRITE (stdout, *) " Read from "//TRIM(filename)
         WRITE (stdout, '((A),(I4))') "   number of grid points   = ", ngrid
         WRITE (stdout, '((A),(I4))') "   number of projectors    = ", nproj

         CALL allocate_atproj_type(types(it), ngrid, nproj)
         types(it)%atsym = atm(it)

         READ (iun_atproj, *) (types(it)%l(i), i=1, nproj)
         WRITE (stdout, '((A))', advance='no') "   ang. mom. of projectors = "
         DO i = 1, nproj
            WRITE (stdout, '(I4)', advance='no') types(it)%l(i)
         END DO
         WRITE (stdout, *)
         WRITE (stdout, *)

         DO i = 1, ngrid
            READ (iun_atproj, *) types(it)%xgrid(i), types(it)%rgrid(i), &
               (types(it)%radial(i, j), j=1, nproj)
         END DO

         CLOSE (iun_atproj)
      END DO

      RETURN
   END SUBROUTINE read_atomproj

   SUBROUTINE broadcast_atomproj()
      use ions_base, only: nsp, ityp, nat
      use io_global, only: ionode, ionode_id
      use mp_world, only: world_comm
      use mp, only: mp_bcast
      use basis, only: natomwfc
      use projections, only: fill_nlmchi, nlmchi
      use wannier, only: n_proj

      implicit none

      integer :: it, i, j, ierr, lmax_wfc, iproj, nwfc, n, l, m, iwfc

      ! Broadcast the data
      CALL mp_bcast(n_proj, ionode_id, world_comm)

      if (atom_proj_ext) then
         DO it = 1, nsp
            i = atproj_types(it)%ngrid
            j = atproj_types(it)%nproj
            CALL mp_bcast(i, ionode_id, world_comm)
            CALL mp_bcast(j, ionode_id, world_comm)
            IF (.NOT. ionode) CALL allocate_atproj_type(atproj_types(it), i, j)
            !
            CALL mp_bcast(atproj_types(it)%atsym, ionode_id, world_comm)
            CALL mp_bcast(atproj_types(it)%xgrid, ionode_id, world_comm)
            CALL mp_bcast(atproj_types(it)%rgrid, ionode_id, world_comm)
            CALL mp_bcast(atproj_types(it)%l, ionode_id, world_comm)
            CALL mp_bcast(atproj_types(it)%radial, ionode_id, world_comm)
         END DO

         ! Update nlmchi
         allocate(nlmchi(n_proj))
         iproj = 0
         DO i = 1, nat
            it = ityp(i)
            DO nwfc = 1, atproj_types(it)%nproj
               l = atproj_types(it)%l(nwfc)
               ! Work out n by looking through the previous projectors for this atom
               n = l + 1
               DO iwfc = 1, nwfc - 1
                  if (atproj_types(it)%l(iwfc) == l) n = n + 1
               END DO

               DO m = 1, 2*l + 1
                  iproj = iproj + 1
                  nlmchi(iproj)%na = i
                  nlmchi(iproj)%n = n
                  nlmchi(iproj)%l = l
                  nlmchi(iproj)%m = m
                  nlmchi(iproj)%ind = m + 2*l + 1
                  nlmchi(iproj)%jj = 0.0d0
                  nlmchi(iproj)%els = ' '
               END DO
            END DO
         END DO
      ELSE
         ! need to access nlmchi, natomwfc, lmax_wfc on each core,
         ! the root node has been filled already
         IF (.NOT. ionode) CALL fill_nlmchi(natomwfc, lmax_wfc)
      END IF

      CALL mp_bcast(natproj, ionode_id, world_comm)
      CALL mp_bcast(nexatproj, ionode_id, world_comm)
      IF (.NOT. ionode) THEN
         ALLOCATE (atproj_excl(n_proj + nexatproj), stat=ierr)
         IF (ierr /= 0) CALL errore('atomic_projectors', 'Error allocating atproj_excl', 1)
      END IF
      CALL mp_bcast(atproj_excl, ionode_id, world_comm)

   END SUBROUTINE

   SUBROUTINE read_atomproj_wrapper(types, print_info_arg)

      use io_global, only: stdout
      USE ions_base, ONLY: nsp, atm, nat, ityp
      use wannier, only: n_proj

      implicit none

      TYPE(atproj_type), INTENT(INOUT) :: types(nsp)
      logical, intent(in), optional :: print_info_arg

      logical :: print_info
      integer :: i, it, l, m, nwfc, iproj

      print_info = .false.
      if (present(print_info_arg)) print_info = print_info_arg

      if (print_info) then
         WRITE (stdout, '(a)') '  Using atomic projectors from dir '//TRIM(atom_proj_dir)
         WRITE (stdout, *) ''
      end if
      call read_atomproj(types)
      n_proj = 0
      DO i = 1, nat
         it = ityp(i)
         DO nwfc = 1, types(it)%nproj
            l = types(it)%l(nwfc)
            DO m = 1, 2*l + 1
               n_proj = n_proj + 1
               if (print_info) then
                  WRITE (stdout, 1000, ADVANCE="no") n_proj, i, types(it)%atsym, nwfc, l
                  WRITE (stdout, '(" m=", i2, ")")') m
               end if
            END DO
         END DO
      END DO
      WRITE (stdout, *) ''

1000  FORMAT(5X, "state #", i4, ": atom ", i3, " (", a3, "), wfc ", i2, " (l=", i1)

   END SUBROUTINE read_atomproj_wrapper

   SUBROUTINE load_upfproj_wrapper(print_info_arg)

      use ions_base, only: atm, nat, ityp
      use io_global, only: stdout
      use basis, only: natomwfc
      use projections, only: nlmchi, fill_nlmchi, compute_mj
      use wannier, only: n_proj
      use noncollin_module, only: lspinorb, noncolin

      implicit none

      logical, intent(in), optional :: print_info_arg

      logical :: print_info
      integer :: lmax_wfc, nwfc

      print_info = .false.
      if (present(print_info_arg)) print_info = print_info_arg

      if (print_info) then
         WRITE (stdout, '(a)') '  Use atomic projectors from UPF'
         WRITE (stdout, *) ''
         WRITE (stdout, '( 5x,"(read from pseudopotential files):"/)')
      end if
      CALL fill_nlmchi(natomwfc, lmax_wfc)
      if (print_info) then
         DO nwfc = 1, natomwfc
            WRITE (stdout, 1000, ADVANCE="no") &
               nwfc, nlmchi(nwfc)%na, atm(ityp(nlmchi(nwfc)%na)), &
               nlmchi(nwfc)%n, nlmchi(nwfc)%l
            IF (lspinorb) THEN
               WRITE (stdout, '(" j=", f3.1, " m_j=", f4.1, ")")') &
                  nlmchi(nwfc)%jj, compute_mj(nlmchi(nwfc)%jj, nlmchi(nwfc)%l, nlmchi(nwfc)%m)
            ELSE IF (noncolin) THEN
               WRITE (stdout, '(" m=", i2, " s_z=", f4.1, ")")') &
                  nlmchi(nwfc)%m, 0.5D0 - INT(nlmchi(nwfc)%ind/(2*nlmchi(nwfc)%l + 2))
            ELSE
               WRITE (stdout, '(" m=", i2, ")")') nlmchi(nwfc)%m
            END IF
         END DO
      end if
      WRITE (stdout, *) ''

      n_proj = natomwfc

1000  FORMAT(5X, "state #", i4, ": atom ", i3, " (", a3, "), wfc ", i2, " (l=", i1)

   END SUBROUTINE

   SUBROUTINE count_number_of_excluded_projectors(print_info_arg)

      use wannier, only: n_proj
      use io_global, only: stdout

      implicit none

      logical, intent(in), optional :: print_info_arg

      logical :: print_info
      integer :: ierr, i, j
      logical :: has_excl_proj
      CHARACTER(len=256) :: err_str

      print_info = .false.
      if (present(print_info_arg)) print_info = print_info_arg


      ALLOCATE (atproj_excl(n_proj), stat=ierr)
      IF (ierr /= 0) CALL errore('atomic_projectors', 'Error allocating atproj_excl', 1)
      atproj_excl = .false.
      !
      DO i = 1, nexatproj_max
         IF (atom_proj_exclude(i) > n_proj) THEN
            WRITE (err_str, *) 'atom_proj_exclude(', i, ') = ', atom_proj_exclude(i), &
               '> total number of projectors (', n_proj, ')'
            CALL errore('pw2wannier90', err_str, i)
         ELSEIF (atom_proj_exclude(i) < 0) THEN
            CYCLE
         ELSE
            atproj_excl(atom_proj_exclude(i)) = .true.
         END IF
      END DO
      !
      nexatproj = COUNT(atproj_excl)
      has_excl_proj = (nexatproj > 0)
      !
      natproj = n_proj
      !
      IF (has_excl_proj) THEN
         if (print_info) then
            WRITE (stdout, *) '    excluded projectors: '
            j = 0 ! how many elements have been written
            DO i = 1, n_proj
               IF (atproj_excl(i)) THEN
                  WRITE (stdout, '(i8)', advance='no') i
                  j = j + 1
                  IF (MOD(j, 10) == 0) WRITE (stdout, *)
               END IF
            END DO
            WRITE (stdout, *) ''
         end if
         n_proj = n_proj - nexatproj
      END IF

   END SUBROUTINE

   SUBROUTINE init_atomproj(lmax_proj, print_info_arg)
      use io_global, only: ionode, stdout
      use mp_world, only: world_comm
      use wannier, only: n_proj
      use ions_base, only: nsp

      implicit none

      integer, intent(out) :: lmax_proj
      logical, intent(in), optional :: print_info_arg

      logical :: print_info
      integer :: nt

      print_info = .false.
      if (present(print_info_arg)) print_info = print_info_arg
      
      if (atom_proj_ext) allocate(atproj_types(nsp))

      if (ionode) then
         if (atom_proj_ext) then
            ! Read the file
            call read_atomproj_wrapper(atproj_types, print_info)
         else
            ! Use atomic projectors from UPF
            call load_upfproj_wrapper(print_info)
         end if

         call count_number_of_excluded_projectors()
      end if

      ! Broadcast
      call broadcast_atomproj()

      ! Populate the atomic projectors table
      CALL init_tab_atproj(world_comm)

      ! Calculating lmax_proj
      lmax_proj = 0
      DO nt = 1, nsp
         lmax_proj = MAX(lmax_proj, MAXVAL(atproj_types(nt)%l))
      END DO

   END SUBROUTINE init_atomproj

   !-----------------------------------------------------------------------
   SUBROUTINE atomic_wfc_ext(ik, wfcatom)
      !-----------------------------------------------------------------------
      !! This routine computes the superposition of atomic wavefunctions
      !! for k-point "ik" - output in "wfcatom".
      !
      !  adapted from PW/src/atomic_wfc.f90, PP/src/atomic_wfc_nc_proj.f90
      !  to use external atomic wavefunctions other than UPF ones.
      !
      USE kinds, ONLY: DP
      USE constants, ONLY: tpi, fpi, pi
      USE cell_base, ONLY: omega, tpiba
      USE ions_base, ONLY: nat, ntyp => nsp, ityp, tau
      USE gvect, ONLY: mill, eigts1, eigts2, eigts3, g
      USE klist, ONLY: xk, igk_k, ngk
      USE wvfct, ONLY: npwx
      USE noncollin_module, ONLY: noncolin, npol, angle1, angle2, lspinorb, &
                                  domag, starting_spin_angle
      USE upf_spinorb, ONLY: rot_ylm, lmaxx, fcoef, lmaxx
      USE wannier, ONLY: n_proj
      !
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ik
      !! k-point index
      COMPLEX(DP), INTENT(OUT) :: wfcatom(npwx, npol, natproj)
      !! Superposition of atomic wavefunctions
      !
      ! ... local variables
      !
      INTEGER :: n_starting_wfc, lmax_wfc, nt, l, nb, na, m, lm, ig, iig, &
                 i0, i1, i2, i3, npw
      REAL(DP), ALLOCATABLE :: qg(:), ylm(:, :), chiq(:, :, :), gk(:, :)
      COMPLEX(DP), ALLOCATABLE :: sk(:), aux(:)
      COMPLEX(DP) :: kphase, lphase
      REAL(DP)    :: arg
      INTEGER :: nwfcm
      !! max number of radial atomic projectors across atoms
      INTEGER :: ierr

      CALL start_clock('atomic_wfc_ext')

      ! calculate max angular momentum required in wavefunctions
      lmax_wfc = 0
      nwfcm = 0
      DO nt = 1, ntyp
         lmax_wfc = MAX(lmax_wfc, MAXVAL(atproj_types(nt)%l))
         nwfcm = MAX(nwfcm, atproj_types(nt)%nproj)
      END DO
      !
      npw = ngk(ik)
      !
      ALLOCATE (ylm(npw, (lmax_wfc + 1)**2), stat=ierr)
      IF (ierr /= 0) CALL errore('atomic_projectors', 'Error allocating ylm', 1)
      ALLOCATE (chiq(npw, nwfcm, ntyp), stat=ierr)
      IF (ierr /= 0) CALL errore('atomic_projectors', 'Error allocating chiq', 1)
      ALLOCATE (gk(3, npw), qg(npw), stat=ierr)
      IF (ierr /= 0) CALL errore('atomic_projectors', 'Error allocating gk/qg', 1)
      !
      DO ig = 1, npw
         iig = igk_k(ig, ik)
         gk(1, ig) = xk(1, ik) + g(1, iig)
         gk(2, ig) = xk(2, ik) + g(2, iig)
         gk(3, ig) = xk(3, ik) + g(3, iig)
         qg(ig) = gk(1, ig)**2 + gk(2, ig)**2 + gk(3, ig)**2
      END DO
      !
      !  ylm = spherical harmonics
      !
      CALL ylmr2((lmax_wfc + 1)**2, npw, gk, qg, ylm)
      !
      ! set now q=|k+G| in atomic units
      !
      DO ig = 1, npw
         qg(ig) = SQRT(qg(ig))*tpiba
      END DO
      !
      CALL interp_atproj(npw, qg, nwfcm, chiq)
      !
      DEALLOCATE (qg, gk)
      ALLOCATE (aux(npw), sk(npw), stat=ierr)
      IF (ierr /= 0) CALL errore('atomic_projectors', 'Error allocating aux/sk', 1)
      !
      wfcatom(:, :, :) = (0.0_DP, 0.0_DP)
      n_starting_wfc = 0
      !
      DO na = 1, nat
         arg = (xk(1, ik)*tau(1, na) + xk(2, ik)*tau(2, na) + xk(3, ik)*tau(3, na))*tpi
         kphase = CMPLX(COS(arg), -SIN(arg), KIND=DP)
         !
         !     sk is the structure factor
         !
         DO ig = 1, npw
            iig = igk_k(ig, ik)
            sk(ig) = kphase*eigts1(mill(1, iig), na)* &
                     eigts2(mill(2, iig), na)* &
                     eigts3(mill(3, iig), na)
         END DO
         !
         nt = ityp(na)
         DO nb = 1, atproj_types(nt)%nproj
            l = atproj_types(nt)%l(nb)
            lphase = (0.D0, 1.D0)**l
            !
            !  only support without spin-orbit coupling
            !  if needed in the future, add code according to
            !  PP/src/atomic_wfc_nc_proj.f90
            !
            CALL atomic_wfc___()
            !
            ! END IF
            !
         END DO
         !
      END DO

      DEALLOCATE (aux, sk, chiq, ylm)

      CALL stop_clock('atomic_wfc_ext')

      RETURN

   CONTAINS

      SUBROUTINE atomic_wfc___()
         !
         ! ... LSDA or nonmagnetic case
         !
         DO m = 1, 2*l + 1
            lm = l**2 + m
            n_starting_wfc = n_starting_wfc + 1
            IF (n_starting_wfc > natproj) CALL errore &
               ('atomic_wfc___', 'internal error: too many wfcs', 1)
            !
            DO ig = 1, npw
               wfcatom(ig, 1, n_starting_wfc) = lphase* &
                                                sk(ig)*ylm(ig, lm)*chiq(ig, nb, nt)
            END DO
            !
         END DO
         !
      END SUBROUTINE atomic_wfc___
      !
   END SUBROUTINE atomic_wfc_ext

   SUBROUTINE interp_atproj(npw, qg, nwfcm, chiq)
      !-----------------------------------------------------------------------
      !
      ! computes chiq: radial fourier transform of atomic projector chi
      !
      !! adapted from upflib/interp_atwfc.f90
      !  to support external projectors
      !
      USE kinds, ONLY: dp
      USE ions_base, ONLY: nsp
      USE uspp_data, ONLY: dq
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN)  :: npw
      INTEGER, INTENT(IN)  :: nwfcm
      REAL(dp), INTENT(IN) :: qg(npw)
      REAL(dp), INTENT(OUT):: chiq(npw, nwfcm, nsp)
      !
      INTEGER :: nt, nb, ig
      INTEGER :: i0, i1, i2, i3
      REAL(dp):: qgr, px, ux, vx, wx
      !
      DO nt = 1, nsp
         DO nb = 1, atproj_types(nt)%nproj
            DO ig = 1, npw
               qgr = qg(ig)
               px = qgr/dq - INT(qgr/dq)
               ux = 1.D0 - px
               vx = 2.D0 - px
               wx = 3.D0 - px
               i0 = INT(qgr/dq) + 1
               i1 = i0 + 1
               i2 = i0 + 2
               i3 = i0 + 3
               chiq(ig, nb, nt) = &
                  tab_at(i0, nb, nt)*ux*vx*wx/6.D0 + &
                  tab_at(i1, nb, nt)*px*vx*wx/2.D0 - &
                  tab_at(i2, nb, nt)*px*ux*wx/2.D0 + &
                  tab_at(i3, nb, nt)*px*ux*vx/6.D0
            END DO
         END DO
      END DO

   END SUBROUTINE interp_atproj

   SUBROUTINE init_tab_atproj(intra_bgrp_comm)
      !-----------------------------------------------------------------------
      !! This routine computes a table with the radial Fourier transform
      !! of the atomic wavefunctions.
      !!
      !! adapted from upflib/init_tab_atwfc.f90
      !! to support external projectors
      !
      USE kinds, ONLY: DP
      USE upf_const, ONLY: fpi
      USE uspp_data, ONLY: nqx, dq
      USE ions_base, ONLY: nsp
      USE cell_base, ONLY: omega
      USE mp, ONLY: mp_sum
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: intra_bgrp_comm
      !
      INTEGER :: nt, nb, iq, ir, l, startq, lastq, ndm, nwfcm, ierr
      !
      REAL(DP), ALLOCATABLE :: aux(:), vchi(:), rab(:)
      REAL(DP) :: vqint, pref, q
      !
      ndm = 0
      nwfcm = 0
      DO nt = 1, nsp
         ndm = MAX(ndm, atproj_types(nt)%ngrid)
         nwfcm = MAX(nwfcm, atproj_types(nt)%nproj)
      END DO
      ALLOCATE (aux(ndm), vchi(ndm), rab(ndm), stat=ierr)
      IF (ierr /= 0) CALL errore('atomic_projectors', 'Error allocating aux/vchi/rab', 1)
      !
      ! chiq = radial fourier transform of atomic orbitals chi
      !
      pref = fpi/SQRT(omega)
      ! needed to normalize atomic wfcs (not a bad idea in general)
      CALL divide(intra_bgrp_comm, nqx, startq, lastq)
      !
      ! nqx = INT( (SQRT(ecutwfc) / dq + 4) )
      ALLOCATE (tab_at(nqx, nwfcm, nsp), stat=ierr)
      IF (ierr /= 0) CALL errore('atomic_projectors', 'Error allocating tab_at', 1)
      tab_at(:, :, :) = 0.0_DP
      !
      DO nt = 1, nsp
         rab = (atproj_types(nt)%xgrid(2) - atproj_types(nt)%xgrid(1))* &
               atproj_types(nt)%rgrid
         DO nb = 1, atproj_types(nt)%nproj
            !
            l = atproj_types(nt)%l(nb)
            !
            DO iq = startq, lastq
               q = dq*(iq - 1)
               CALL sph_bes(atproj_types(nt)%ngrid, atproj_types(nt)%rgrid, q, l, aux)
               DO ir = 1, atproj_types(nt)%ngrid
                  vchi(ir) = atproj_types(nt)%radial(ir, nb)*aux(ir)*atproj_types(nt)%rgrid(ir)
               END DO
               CALL simpson(atproj_types(nt)%ngrid, vchi, rab, vqint)
               tab_at(iq, nb, nt) = vqint*pref
            END DO
            !
         END DO
      END DO
      !
      CALL mp_sum(tab_at, intra_bgrp_comm)
      !
      DEALLOCATE (aux, vchi, rab)
      !
      RETURN
      !
   END SUBROUTINE init_tab_atproj

   SUBROUTINE deallocate_atproj
      IMPLICIT NONE

      IF (ALLOCATED(tab_at)) DEALLOCATE (tab_at)

      RETURN
   END SUBROUTINE deallocate_atproj

end module atomic_projectors
