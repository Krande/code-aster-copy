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
subroutine lcverr(dy, ddy, nr, typ, err)
    implicit none
!       MODULE DE CALCUL DE L'ERREUR DE CONVERGENCE
!       IN  DY     :    VECTEUR SOLUTION
!           DDY    :    VECTEUR CORRECTION SUR LA SOLUTION
!           NR     :    DIMENSION DE DY DDY
!           TYP    :    TYPE D'ERREUR A CALCULER
!                               0 = IDDYI/DYII     < EPS (POUR TOUT I)
!                               1 = IIDDYII/IIDYII < EPS
!                               2 = IIDDYI/DYIII   < EPS
!       OUT ERR    :    VECTEUR ERREUR
!       ----------------------------------------------------------------
#include "asterc/r8prem.h"
    real(kind=8) :: zero
!-----------------------------------------------------------------------
    integer(kind=8) :: i
!-----------------------------------------------------------------------
    parameter(zero=0.d0)
    integer(kind=8) :: n, nd, nr, typ
    real(kind=8) :: dy(*), ddy(*)
    real(kind=8) :: err(*), e(50)
!       ----------------------------------------------------------------
    common/tdim/n, nd
!       ----------------------------------------------------------------
!
!       ERREUR(I) =  !DDYI/DYI! < EPS
!
    if (typ .eq. 0) then
        err(1) = 0.d0
        do i = 1, nr
!                IF(DY(I).EQ.ZERO) THEN
            if (abs(dy(i)) .lt. r8prem()) then
                err(i) = abs(ddy(i))
            else
                err(i) = abs(ddy(i)/dy(i))
            end if
            err(1) = max(err(1), err(i))
        end do
!
!       ERREUR = !!DDY!!/!!DY!! < EPS
!
    else if (typ .eq. 1) then
        e(1) = norm2(ddy(1:nr))
        e(2) = norm2(dy(1:nr))
        if (e(2) .eq. zero) then
            err(1) = e(1)
        else
            err(1) = e(1)/e(2)
        end if
!
!       ERREUR = !!DDYI/DYI!! < EPS
!
    else if (typ .eq. 2) then
        do i = 1, nr
            if (abs(dy(i)) .lt. r8prem()) then
                e(i) = ddy(i)
            else
                e(i) = ddy(i)/dy(i)
            end if
        end do
        err(1) = norm2(e(1:nr))
!
    end if
!
end subroutine
