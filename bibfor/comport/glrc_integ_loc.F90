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

subroutine glrc_integ_loc(lambda, deuxmu, seuil, alf, &
                          alfmc, gmt, gmc, cof1, &
                          vim, q2d, qff, tr2d, eps33, &
                          de33d1, de33d2, ksi2d, dksi1, dksi2, &
                          da1, da2, kdmax, told, codret, &
                          emp)
!
! person_in_charge: sebastien.fayolle at edf.fr
! aslint: disable=W1504
    implicit none
!
#include "asterf_types.h"
#include "asterfort/glrc_calc_eps33.h"
    integer(kind=8) :: kdmax, codret
    real(kind=8) :: vim(*), gmt, gmc, tr2d, eps33
    real(kind=8) :: lambda, deuxmu, seuil, alf, qff(2), told
    real(kind=8) :: de33d1, de33d2, ksi2d, dksi1, dksi2, da1, da2, treps, treps2
    real(kind=8) :: cof1(2), q2d(2)
!
! ----------------------------------------------------------------------
!
!      CALCUL D'ENDOMMAGEMENT (DA1 ET DA2) ET DE LA COMPOSANTE DE
!       DEFORMATION EPS33, APPELE PAR "LCGLDM" (LOI GLRC_DM)
!
! IN:
!       LAMBDA  : PARAMETRE D ELASTICITE - MEMBRANE
!       DEUXMU  : PARAMETRE D ELASTICITE - MEMBRANE
!       GMT     : PARAMETRE GAMMA POUR LA MEMBRANE EN TRACTION
!       GMC     : PARAMETRE GAMMA POUR LA MEMBRANE EN COMPRESSION
!       SEUIL   : INITIAL MEMBRANE
!       ALF     : PARAMETRE DE SEUIL FLEXION
!       VIM     : VARIABLES INTERNES EN T-
!       Q2D     : PARTIE CONSTANTE DU RESIDU MEMBRANE
!       QFF(2)  : PARTIE CONSTANTE DU RESIDU FLEXION
!       TR2D    : EPS11 + EPS22
!       COF1    : PARAMETRE
! OUT:
!       DA1     : ENDOMMAGEMENT SUR LA PARTIE 1 DE LA PLAQUE
!       DA2     : ENDOMMAGEMENT SUR LA PARTIE 2 DE LA PLAQUE
!       KSI2D   : VALEUR DE LA FONCTION CARACTERISTIQUE DE L'ENDOM.
!       EPS33   : COMPOSANTE 33 DE LA DEFORMATION
!       DE33D1  : DERIVEE DE EPS33 PAR RAPPORT A DA1
!       DE33D2  : DERIVEE DE EPS33 PAR RAPPORT A DA2
!       ELAS    : .TRUE. SI ELASTIQUE
!       ELAS1   : .TRUE. SI DA1 == VIM(1)
!       ELAS2   : .TRUE. SI DA2 == VIM(2)
!       CODRET  : CODE RETOUR DE L'INTEGRATION INTEGRATION DU
!                 0 => PAS DE PROBLEME
!                 1 => ABSENCE DE CONVERGENCE
! ----------------------------------------------------------------------
!
    aster_logical :: lconv1, lconv2
!
    integer(kind=8) :: i
!
    real(kind=8) :: qm1, qm2
    real(kind=8) :: rd1, rd2, dr1d, dr2d, dd1, dd2, seuilr
    real(kind=8) :: alfmc, emp(2), cof2(2), dq2d(2)
!
    call glrc_calc_eps33(lambda, deuxmu, alfmc, gmt, gmc, &
                         tr2d, da1, da2, eps33, de33d1, &
                         de33d2, ksi2d, dksi1, dksi2, cof1, &
                         q2d, emp, cof2, dq2d)
!
    treps = tr2d+eps33
    treps2 = treps**2
!
!----------CONTRIBUTION DE MEMBRANE-----
    qm1 = 0.5d0*cof1(1)*treps2+q2d(1)
    qm2 = 0.5d0*cof1(2)*treps2+q2d(2)
    rd1 = qm1/(1.0d0+da1)**2-seuil
    rd2 = qm2/(1.0d0+da2)**2-seuil
!
!----CONTRIBUTION DES COURBURES---------
    rd1 = rd1+qff(1)/(alf+da1)**2
    rd2 = rd2+qff(2)/(alf+da2)**2
!
    seuilr = max(seuil, 1.0d-6)
!
!-----VERIFIER SI LE SEUIL EST ATTEINT
    lconv1 = rd1 .lt. 0.0d0
    lconv2 = rd2 .lt. 0.0d0
!
    if ((.not. lconv1 .and. da1 .ge. vim(1)) .or. (.not. lconv2 .and. da2 .ge. vim(2))) then
!
        do i = 1, kdmax
!
            dd1 = 0.0d0
            dd2 = 0.0d0
!
            if (rd1 .ge. 0.0d0) then
                dr1d = (cof1(1)*treps*de33d1+0.5d0*cof2(1)*treps2+dq2d(1))/(1.d0+da1)**2 &
                       -2.d0*qm1/(1.d0+da1)**3-2.d0*qff(1)/(alf+da1)**3
                if (abs(dr1d) .lt. 1.0d-14) then
                    dd1 = 0.0d0
                else
                    dd1 = -rd1/dr1d
                end if
            end if
!
            if (rd2 .ge. 0.0d0) then
                dr2d = (cof1(2)*treps*de33d2+0.5d0*cof2(2)*treps2+dq2d(2))/(1.d0+da2)**2 &
                       -2.d0*qm2/(1.d0+da2)**3-2.d0*qff(2)/(alf+da2)**3
!
                if (abs(dr2d) .lt. 1.0d-14) then
                    dd2 = 0.0d0
                else
                    dd2 = -rd2/dr2d
                end if
            end if
!
            if (((abs(dd1*rd1) .lt. told*seuilr) .or. ((rd1 .lt. 0.0d0 .and. da1 .le. vim(1)))) &
                .and. &
                ((abs(dd2*rd2) .lt. told*seuilr) .or. ((rd2 .lt. 0.0d0 .and. da2 .le. vim(2))) &
                 )) goto 10
!
            da1 = da1+dd1
            da2 = da2+dd2
!
            if (da1 .lt. 0.0d0 .and. rd1 .lt. 0.0d0) da1 = vim(1)
            if (da2 .lt. 0.0d0 .and. rd2 .lt. 0.0d0) da2 = vim(2)
!
            call glrc_calc_eps33(lambda, deuxmu, alfmc, gmt, gmc, &
                                 tr2d, da1, da2, eps33, de33d1, &
                                 de33d2, ksi2d, dksi1, dksi2, cof1, &
                                 q2d, emp, cof2, dq2d)
!
            treps = tr2d+eps33
            treps2 = treps**2
!
!----------CONTRIBUTION DE MEMBRANE-----
            qm1 = 0.5d0*cof1(1)*treps2+q2d(1)
            qm2 = 0.5d0*cof1(2)*treps2+q2d(2)
            rd1 = qm1/(1.0d0+da1)**2-seuil
            rd2 = qm2/(1.0d0+da2)**2-seuil
!
!----CONTRIBUTION DES COURBURES---------
            rd1 = rd1+qff(1)/(alf+da1)**2
            rd2 = rd2+qff(2)/(alf+da2)**2
!
        end do
!
!    NON CONVERGENCE POUR LE NOMBRE MAXIMAL D ITERATION PRESCRIT
        codret = 1
!
10      continue
    end if
end subroutine
