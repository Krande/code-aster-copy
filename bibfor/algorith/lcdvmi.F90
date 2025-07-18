! --------------------------------------------------------------------
! Copyright (C) 1991 - 2025 - EDF R&D - www.code-aster.org
! This file is part of code_aster.
!
! code_aster is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! code_aster is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with code_aster.  If not, see <http://www.gnu.org/licenses/>.
! --------------------------------------------------------------------
!
subroutine lcdvmi(sigma, y, f, dfds, d2fds, &
                  seq)
    implicit none
    real(kind=8) :: y, sigma(6), f, seq, dfds(6), d2fds(6, 6)
!-----------------------------------------------------------------------
!  ROUTINE D EVALUATION DU CRITERE DE VON-MISES ISOTROPE
!  ET DE SES DERIVEES PARTIELLES PAR RAPPORT A SIGMA
!
!        !!! DONNEES ET RESULTATS PAR PG !!!
!-----------------------------------------------------------------------
!  ENTREES
!    SIGMA : VECTEUR CONTRAINTE
!    Y     : VALEUR COURANTE DU CRITERE
!
!  SORTIES
!    F     : CRITERE (NUL A LA CONVERGENCE)
!    DFDS  : DERIVEE DE F PAR RAPPORT AUX CONTRAINTES
!    D2FDS : DERIVEE SECONDE DE F PAR RAPPORT AUX CONTRAINTES
!    SEQ   : CONTRAINTE EQUIVALENTE (= Y A LA CONVERGENCE)
!-----------------------------------------------------------------------
    integer(kind=8) :: i, j
    real(kind=8) :: s23, s31, s12, t12, t23, t31
!-----------------------------------------------------------------------
!  CALCUL DE F
!  -----------
!
    s23 = (sigma(2)-sigma(3))*(sigma(2)-sigma(3))
    s31 = (sigma(3)-sigma(1))*(sigma(3)-sigma(1))
    s12 = (sigma(1)-sigma(2))*(sigma(1)-sigma(2))
    t12 = sigma(4)*sigma(4)
    t23 = sigma(5)*sigma(5)
    t31 = sigma(6)*sigma(6)
!
    seq = sqrt((s23+s31+s12)/2.d0+1.5d0*(t12+t23+t31))
!
    f = seq-y
!
    if (seq .lt. 1.d-9) goto 20
!
!  CALCUL DE DF/DSIG (DFDS(6))
!  ---------------------------
!
    dfds(1) = (2.d0*sigma(1)-sigma(2)-sigma(3))/(2.d0*seq)
    dfds(2) = (2.d0*sigma(2)-sigma(3)-sigma(1))/(2.d0*seq)
    dfds(3) = (2.d0*sigma(3)-sigma(1)-sigma(2))/(2.d0*seq)
    dfds(4) = 3.d0*sigma(4)/seq/2.d0
    dfds(5) = 3.d0*sigma(5)/seq/2.d0
    dfds(6) = 3.d0*sigma(6)/seq/2.d0
!
    d2fds(1, 1) = (1.0d0-dfds(1)*dfds(1))/seq
    d2fds(1, 2) = (-0.5d0-dfds(1)*dfds(2))/seq
    d2fds(1, 3) = (-0.5d0-dfds(1)*dfds(3))/seq
    d2fds(1, 4) = (-dfds(1)*dfds(4))/seq
    d2fds(1, 5) = (-dfds(1)*dfds(5))/seq
    d2fds(1, 6) = (-dfds(1)*dfds(6))/seq
!
    d2fds(2, 2) = (1.0d0-dfds(2)*dfds(2))/seq
    d2fds(2, 3) = (-0.5d0-dfds(2)*dfds(3))/seq
    d2fds(2, 4) = (-dfds(2)*dfds(4))/seq
    d2fds(2, 5) = (-dfds(2)*dfds(5))/seq
    d2fds(2, 6) = (-dfds(2)*dfds(6))/seq
!
    d2fds(3, 3) = (1.0d0-dfds(3)*dfds(3))/seq
    d2fds(3, 4) = (-dfds(3)*dfds(4))/seq
    d2fds(3, 5) = (-dfds(3)*dfds(5))/seq
    d2fds(3, 6) = (-dfds(3)*dfds(6))/seq
!
    d2fds(4, 4) = (3.d0-dfds(4)*dfds(4))/seq
    d2fds(4, 5) = (-dfds(4)*dfds(5))/seq
    d2fds(4, 6) = (-dfds(4)*dfds(6))/seq
!
    d2fds(5, 5) = (3.d0-dfds(5)*dfds(5))/seq
    d2fds(5, 6) = (-dfds(5)*dfds(6))/seq
!
    d2fds(6, 6) = (3.d0-dfds(6)*dfds(6))/seq
!
    do i = 1, 6
        do j = i, 6
            d2fds(j, i) = d2fds(i, j)
        end do
    end do
!
20  continue
!
end subroutine
