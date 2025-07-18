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

subroutine cfadju(alias, ksi1, ksi2, toleou, iproj)
!
! person_in_charge: mickael.abbas at edf.fr
!
    implicit none
#include "asterfort/assert.h"
    character(len=8) :: alias
    real(kind=8) :: ksi1, ksi2, toleou
    integer(kind=8) :: iproj
!
! ----------------------------------------------------------------------
!
! ROUTINE CONTACT (METHODE DISCRETE- APPARIEMENT)
!
! AJUSTE LES COORDONNES PARAMETRIQUES POUR RESTER DANS LA MAILLE
!
! ----------------------------------------------------------------------
!
!
! IN  ALIAS  : TYPE DE L'ELEMENT
!               'SE2','SE3'
!               'TR3','TR6'
!               'QU4','QU8','QU9'
! I/O KSI1   : POINT DE CALCUL SUIVANT KSI1 DES
!               FONCTIONS DE FORME ET LEURS DERIVEES
! I/O KSI2   : POINT DE CALCUL SUIVANT KSI2 DES
!               FONCTIONS DE FORME ET LEURS DERIVEES
! IN  TOLEOU : TOLERANCE POUR PROJECTION HORS SEGMENT
! OUT IPROJ  : VAUT 0 SI POINT PROJETE DANS L'ELEMENT
!                   1 SI POINT PROJETE DANS LA ZONE DEFINIE PAR TOLEOU
!                   2 SI POINT PROJETE EN DEHORS (EXCLUS)
!
! ----------------------------------------------------------------------
!
    real(kind=8) :: ecart, ksi1e, ksi2e, k1pk2, k2mk1
    integer(kind=8) :: izone
!
!   tolerances --- absolue et relative --- pour determiner si deux distances sont egales
    real(kind=8), parameter :: atol = 1.e-12
    real(kind=8), parameter :: rtol = 1.e-12
    aster_logical :: near
!
! ----------------------------------------------------------------------
!
    iproj = 0
    ecart = -1.d0
!
!   calcul de ksi1 + ksi2 et ksi2 - ksi1 (utilises pour le cas du triangle)
    k1pk2 = ksi1+ksi2
    k2mk1 = ksi2-ksi1
!
    if (alias(1:2) .eq. 'SE') then
!
!       premier ajustement : on positionne le point de coordonnees (ksi1, ksi2) sur
!       le bord, s'il est à une distance (normalisee) inférieure a atol du bord
!
!       ksi1 = -1 ?
        if (abs(ksi1+1.d0) .le. atol) ksi1 = -1.d0
!       ksi1 = +1 ?
        if (abs(ksi1-1.d0) .le. atol) ksi1 = +1.d0
!
        if ((ksi1 .ge. -1.d0) .and. (ksi1 .le. 1.d0)) then
            goto 999
        end if
!
! --- CALCUL DE L'ECART
!
        ecart = abs(ksi1)-1.d0
!
! --- RABATTEMENT
!
        iproj = 1
!
        if (ksi1 .lt. -1.d0) then
            ksi1 = -1.d0
        else if (ksi1 .gt. 1.d0) then
            ksi1 = 1.d0
        end if
!
!       ecart est-il egal a toleou ?
        near = abs(ecart-toleou) .le. (atol+toleou*rtol)
!
        if (ecart .gt. toleou .and. .not. near) then
            iproj = 2
        end if
!
    else if (alias(1:2) .eq. 'TR') then
!
!       premier ajustement : on positionne le point de coordonnees (ksi1, ksi2) sur
!       le bord, s'il est à une distance (normalisee) inférieure a atol du bord
!
!       ksi1 = 0 ?
        if (abs(ksi1) .le. atol) ksi1 = 0.d0
!       ksi2 = 0 ?
        if (abs(ksi2) .le. atol) ksi2 = 0.d0
!       ksi1 + ksi2 = 1 ?
        if (abs(k1pk2-1.d0) .le. atol) k1pk2 = +1.d0
!       ksi2 - ksi1 = -1 ?
        if (abs(k2mk1+1.d0) .le. atol) k2mk1 = -1.d0
!       ksi2 - ksi1 = 1 ?
        if (abs(k2mk1-1.d0) .le. atol) k2mk1 = +1.d0
!
        if ((ksi1 .ge. 0.d0) .and. (ksi2 .ge. 0.d0) .and. (k1pk2 .le. 1.d0)) then
            goto 999
        end if
!
! --- SECTEUR CONCERNE
!
        izone = 0
        if (ksi1 .lt. 0.d0) then
            if (ksi2 .lt. 0.d0) then
                izone = 1
            else if ((ksi2 .ge. 0.d0) .and. (ksi2 .le. 1.d0)) then
                izone = 2
            else if (ksi2 .gt. 1.d0) then
                izone = 3
            else
                ASSERT(.false.)
            end if
        end if
        if (ksi2 .lt. 0.d0) then
            if (ksi1 .lt. 0.d0) then
                izone = 1
            else if ((ksi1 .ge. 0.d0) .and. (ksi1 .le. 1.d0)) then
                izone = 8
            else if (ksi1 .gt. 1.d0) then
                izone = 7
            else
                ASSERT(.false.)
            end if
        end if
        if (ksi1 .ge. 0.d0) then
            if (k2mk1 .gt. 1.d0) then
                izone = 4
            elseif ((k1pk2 .gt. 1.d0) .and. (k2mk1 .ge. -1.d0) &
                    .and. (k2mk1 .le. 1.d0)) then
                izone = 5
                ksi1e = 5.d-1*(1.d0+ksi1-ksi2)
                ksi2e = 5.d-1*(1.d0-ksi1+ksi2)
            else if ((ksi2 .ge. 0.d0) .and. (k2mk1 .lt. -1.d0)) then
                izone = 6
            end if
        end if
!
! --- CALCUL DE L'ECART
!
        if (izone .eq. 1) then
            ecart = sqrt(abs(ksi1)*abs(ksi1)+abs(ksi2)*abs(ksi2))
        else if (izone .eq. 2) then
            ecart = sqrt(abs(ksi1)*abs(ksi1))
        else if (izone .eq. 3 .or. izone .eq. 4) then
            ecart = sqrt(abs(ksi1)*abs(ksi1)+(ksi2-1.d0)*(ksi2-1.d0))
        else if (izone .eq. 5) then
            ecart = sqrt((ksi1-ksi1e)*(ksi1-ksi1e)+(ksi2-ksi2e)*(ksi2-ksi2e))
        else if (izone .eq. 6 .or. izone .eq. 7) then
            ecart = sqrt(abs(ksi2)*abs(ksi2)+(ksi1-1.d0)*(ksi1-1.d0))
        else if (izone .eq. 8) then
            ecart = sqrt(abs(ksi2)*abs(ksi2))
        else
            ASSERT(.false.)
        end if
!
! --- RABATTEMENT
!
        iproj = 1
!
        if (izone .eq. 1) then
            ksi1 = 0.d0
            ksi2 = 0.d0
        else if (izone .eq. 2) then
            ksi1 = 0.d0
        else if (izone .eq. 3 .or. izone .eq. 4) then
            ksi1 = 0.d0
            ksi2 = 1.d0
        else if (izone .eq. 5) then
            ksi1 = ksi1e
            ksi2 = ksi2e
        else if (izone .eq. 6 .or. izone .eq. 7) then
            ksi1 = 1.d0
            ksi2 = 0.d0
        else if (izone .eq. 8) then
            ksi2 = 0.d0
        end if
!
!       ecart est-il egal a toleou ?
        near = abs(ecart-toleou) .le. (atol+toleou*rtol)
!
        if (ecart .gt. toleou .and. .not. near) then
            iproj = 2
        end if
!
    else if (alias(1:2) .eq. 'QU') then
!
!       premier ajustement : on positionne le point de coordonnees (ksi1, ksi2) sur
!       le bord, s'il est à une distance (normalisee) inférieure a atol du bord
!
!       ksi1 = -1 ?
        if (abs(ksi1+1.d0) .le. atol) ksi1 = -1.d0
!       ksi1 = +1 ?
        if (abs(ksi1-1.d0) .le. atol) ksi1 = +1.d0
!       ksi2 = -1 ?
        if (abs(ksi2+1.d0) .le. atol) ksi2 = -1.d0
!       ksi2 = +1 ?
        if (abs(ksi2-1.d0) .le. atol) ksi2 = +1.d0
!
!        if ((abs(ksi1).le.1.d0) .and. (abs(ksi2).le.1.d0)) then
        if ((ksi1 .ge. -1.d0) .and. (ksi1 .le. 1.d0) .and. &
            (ksi2 .ge. -1.d0) .and. (ksi2 .le. 1.d0)) then
            goto 999
        end if
!
! --- SECTEUR CONCERNE
!
        izone = 0
        if (ksi1 .lt. -1.d0) then
            if (ksi2 .lt. -1.d0) then
                izone = 1
!
            else if ((ksi2 .ge. -1.d0) .and. (ksi2 .le. 1.d0)) then
                izone = 2
            else if (ksi2 .gt. 1.d0) then
                izone = 3
            else
                ASSERT(.false.)
            end if
        end if
        if (ksi1 .gt. 1.d0) then
            if (ksi2 .lt. -1.d0) then
                izone = 7
            else if ((ksi2 .ge. -1.d0) .and. (ksi2 .le. 1.d0)) then
                izone = 6
            else if (ksi2 .gt. 1.d0) then
                izone = 5
            else
                ASSERT(.false.)
            end if
        end if
        if ((ksi1 .ge. -1.d0) .and. (ksi1 .le. 1.d0)) then
            if (ksi2 .lt. -1.d0) then
                izone = 8
            else if (ksi2 .gt. 1.d0) then
                izone = 4
            end if
        end if
!
! --- CALCUL DE L'ECART
!
        if (izone .eq. 1) then
            ecart = sqrt((abs(ksi1)-1.d0)*(abs(ksi1)-1.d0)+(abs(ksi2)-1.d0)*(abs(ksi2)-1.d0))
        else if (izone .eq. 2) then
            ecart = sqrt((abs(ksi1)-1.d0)*(abs(ksi1)-1.d0))
        else if (izone .eq. 3) then
            ecart = sqrt((abs(ksi1)-1.d0)*(abs(ksi1)-1.d0)+(abs(ksi2)-1.d0)*(abs(ksi2)-1.d0))
        else if (izone .eq. 4) then
            ecart = sqrt((abs(ksi2)-1.d0)*(abs(ksi2)-1.d0))
        else if (izone .eq. 5) then
            ecart = sqrt((abs(ksi1)-1.d0)*(abs(ksi1)-1.d0)+(abs(ksi2)-1.d0)*(abs(ksi2)-1.d0))
        else if (izone .eq. 6) then
            ecart = sqrt((abs(ksi1)-1.d0)*(abs(ksi1)-1.d0))
        else if (izone .eq. 7) then
            ecart = sqrt((abs(ksi1)-1.d0)*(abs(ksi1)-1.d0)+(abs(ksi2)-1.d0)*(abs(ksi2)-1.d0))
        else if (izone .eq. 8) then
            ecart = sqrt((abs(ksi2)-1.d0)*(abs(ksi2)-1.d0))
        else
            ASSERT(.false.)
        end if
!
! --- RABATTEMENT
!
        iproj = 1
!
        if (izone .eq. 1) then
            ksi1 = -1.d0
            ksi2 = -1.d0
        else if (izone .eq. 2) then
            ksi1 = -1.d0
        else if (izone .eq. 3) then
            ksi1 = -1.d0
            ksi2 = 1.d0
        else if (izone .eq. 4) then
            ksi2 = 1.d0
        else if (izone .eq. 5) then
            ksi1 = 1.d0
            ksi2 = 1.d0
        else if (izone .eq. 6) then
            ksi1 = 1.d0
        else if (izone .eq. 7) then
            ksi1 = 1.d0
            ksi2 = -1.d0
        else if (izone .eq. 8) then
            ksi2 = -1.d0
        end if
!
!       ecart est-il egal a toleou ?
        near = abs(ecart-toleou) .le. (atol+toleou*rtol)
!
        if (ecart .gt. toleou .and. .not. near) then
            iproj = 2
        end if
!
    else
        ASSERT(.false.)
    end if
!
999 continue
!
end subroutine
