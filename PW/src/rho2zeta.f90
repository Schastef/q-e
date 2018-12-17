!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE rho2zeta( rho, rho_core, nrxx, nspin, iop )
  !---------------------------------------------------------------------------
  ! ... For the collinear case:
  ! ... if ( iopi == 1 )  transform the spin up spin down charge density 
  ! ...                   rho(*,is) into :
  !
  ! ...                      rho(*,1) = ( rho_up + rho_dw ) and
  ! ...                      rho(*,2) = ( rho_up - rho_dw ) / rho_tot = zeta
  !
  ! ... if ( iopi == -1)  do the opposit transformation
  ! 
  ! ... and the analogous operation for the non-collinear case.
  !
  USE constants, ONLY : eps32
  USE io_global, ONLY : stdout
  USE kinds,     ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER, intent(in)    :: iop, nspin, nrxx
  REAL(DP), intent(in)   :: rho_core(nrxx)
  REAL(DP), intent(inout):: rho(nrxx,nspin)
  !
  INTEGER :: is
  !
  IF ( nspin == 1 ) RETURN  
  !
  IF ( iop == -1 ) THEN
     !
     DO is = 2, nspin
        !
        WHERE( rho(:,1)+rho_core(:) > eps32 )
           !
           rho(:,is) = rho(:,is)*(rho(:,1)+rho_core(:))
           !
        ELSEWHERE
           !
           rho(:,is) = 0.D0
              !
        END WHERE
        !
     END DO
     !
  ELSE IF ( iop == 1 ) THEN
        !
     DO is = 2, nspin
        !
        WHERE( rho(:,1)+rho_core(:) > eps32 )
           !
           rho(:,is) = rho(:,is) / ( rho(:,1)+rho_core(:) )
           !
        ELSEWHERE
           !
           rho(:,is) = 0.D0
           !
        END WHERE
        !
     END DO
     !
  ELSE
     !
     WRITE( stdout , '(5X,"iop = ",I5)' ) iop
     !
     CALL errore( 'rho2zeta', 'wrong iop', 1 )
     !
  END IF
  !
END SUBROUTINE rho2zeta
