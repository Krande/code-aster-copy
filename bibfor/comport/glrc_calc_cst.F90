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

subroutine glrc_calc_cst(lamf, muf, alf, gf, efp, qff)
!
! person_in_charge: sebastien.fayolle at edf.fr
!
    implicit none
    real(kind=8) :: lamf, muf, alf, gf, trot, trot2
    real(kind=8) :: efp(2)
!
!----------------------------------------------------------------------
!        CALCUL DES CONSTANTES INDEPENDANTES DE DA1, DA2 ET EPS33
!
!     IN :
!         LAMBDA : PARAMETRE D ELASTICITE - MEMBRANE
!         MU     : PARAMETRE D ELASTICITE - MEMBRANE
!         LAMF   : PARAMETRE D ELASTICITE - FLEXION
!         MUF    : PARAMETRE D ELASTICITE - FLEXION
!         ALF    : PARAMETRE DE SEUIL FLEXION
!         GF     : PARAMETRE GAMMA POUR LA FLEXION
!         GMT    : PARAMETRE GAMMA POUR LA MEMBRANE EN TRACTION
!         GMC    : PARAMETRE GAMMA POUR LA MEMBRANE EN COMPRESSION
!         DELTA  : PARAMETRE DE COUPLAGE MEMBRANE-FLEXION
!         TROT   : TRACE DE KAPPA (FLEXION)
!         EMP(2) : VALEURS PROPRES DE EPS 2D
!         EFP(2) : VALEURS PROPRES DE KAPPA
!
!     OUT :
!         QFF(2) : TF
!         QM     : TM
!         COF1   : INTERMEDIAIRE DE CALCUL
!         Q2D    : INTERMEDIAIRE DE CALCUL
!         GI(2)  : PARTIE DE LA DERIVEE DE KSI(EMP) PAR RAPPORT A DA
!         GTR2   : INTERMEDIAIRE DE CALCUL
!         GI(2)  : INTERMEDIAIRE DE CALCUL
!----------------------------------------------------------------------
!
    integer(kind=8) :: k
    real(kind=8) :: qff(2), gf1, gf2
!
!-- ICI ON SUPPOSE QUE GF1=GF2, CE QUI N EST PAS NECESSAIRE
    gf1 = gf
    gf2 = gf
!
    trot = efp(1)+efp(2)
    trot2 = trot**2
!
! -------- CALCUL DE QFF --------------------
!
    if (trot .gt. 0.0d0) then
        qff(1) = 0.5d0*lamf*trot2
        qff(2) = 0.0d0
    else
        qff(1) = 0.0d0
        qff(2) = 0.5d0*lamf*trot2
    end if
!
    do k = 1, 2
        if (efp(k) .gt. 0.0d0) then
            qff(1) = qff(1)+muf*efp(k)**2
        else
            qff(2) = qff(2)+muf*efp(k)**2
        end if
    end do
!
    qff(1) = alf*qff(1)*(1.0d0-gf1)
    qff(2) = alf*qff(2)*(1.0d0-gf2)
!
end subroutine
