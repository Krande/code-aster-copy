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
subroutine rvpstd(valee, type, codir, valdir, valeq)
    implicit none
!
!
#include "asterc/r8vide.h"
    real(kind=8) :: valee(*), valeq(*), valdir(*)
    character(len=2) :: type
    integer(kind=8) :: codir, indir1(3), indir2(4), indir3(3)
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i
!-----------------------------------------------------------------------
    data indir1/2, 1, 3/
    data indir2/2, 1, 4, 3/
    data indir3/2, 3, 1/
!
!**********************************************************************
!
!  OPERATION REALISEE
!  ------------------
!
!     CALCUL DE LA TRACE NORMALE EN UN POINT
!
!  ARGUMENTS EN ENTREE
!  -------------------
!
!     VALEE  : TABLE DES VALEUR DES CMP (REPERE GLOBAL)
!     TYPE     VAUT 'V3' POUR LES VECTEURS ET
!              'T2' POUR LES TENSEUR 2X2 ET 'T3' POUR LES 3X3
!     CODIR  : CODE LES DIRECTIONS ACTIVES
!     VALDIR : VALEURS DU VECTEUR DIRECTION (TJS X, Y, Z)
!
!  ARGUMENTS EN SORTIE
!  -------------------
!
!     VALEE  : TABLE DES VALEUR DES TRACES(REPERE GLOBAL)
!
!**********************************************************************
!
!==================== CORPS DE LA ROUTINE =============================
!
    if (codir .eq. 1) then
!
!     /* DIRECTION ACTIVE : X */
!
        if (type .eq. 'V3') then
!
            if (valee(1) .eq. r8vide()) then
                valeq(1) = 0.d0
            else
                valeq(1) = valee(1)*valdir(1)
            end if
!
        else if (type .eq. 'T3') then
!
            do i = 1, 3
                if (valee(i) .eq. r8vide()) then
                    valeq(i) = 0.d0
                else
                    valeq(i) = valee(i)*valdir(1)
                end if
            end do
!
        else
!
            do i = 1, 4
                if (valee(i) .eq. r8vide()) then
                    valeq(i) = 0.d0
                else
                    valeq(i) = valee(i)*valdir(1)
                end if
            end do
!
        end if
!
    else if (codir .eq. 2) then
!
!     /* DIRECTION ACTIVE : Y */
!
        if (type .eq. 'V3') then
!
            if (valee(1) .eq. r8vide()) then
                valeq(1) = 0.d0
            else
                valeq(1) = valee(1)*valdir(2)
            end if
!
        else if (type .eq. 'T3') then
!
            do i = 1, 3
                if (valee(indir1(i)) .eq. r8vide()) then
                    valeq(i) = 0.d0
                else
                    valeq(i) = valee(indir1(i))*valdir(2)
                end if
            end do
!
        else
!
            do i = 1, 4
                if (valee(indir2(i)) .eq. r8vide()) then
                    valeq(i) = 0.d0
                else
                    valeq(i) = valee(indir2(i))*valdir(2)
                end if
            end do
!
        end if
!
    else if (codir .eq. 3) then
!
!     /* DIRECTION ACTIVE : Z */
!
        if (type .eq. 'V3') then
!
            if (valee(1) .eq. r8vide()) then
                valeq(1) = 0.d0
            else
                valeq(1) = valee(1)*valdir(3)
            end if
!
        else
!
            do i = 1, 3
                if (valee(indir3(i)) .eq. r8vide()) then
                    valeq(i) = 0.d0
                else
                    valeq(i) = valee(indir3(i))*valdir(3)
                end if
            end do
!
        end if
!
    else if (codir .eq. 4) then
!
!     /* DIRECTION ACTIVE : X,Y */
!
        if (type .eq. 'V3') then
!
            do i = 1, 2
                if (valee(i) .eq. r8vide()) valee(i) = 0.d0
            end do
            valeq(1) = valee(1)*valdir(1)+valee(2)*valdir(2)
!
        else if (type .eq. 'T3') then
!
            do i = 1, 5
                if (valee(i) .eq. r8vide()) valee(i) = 0.d0
            end do
            valeq(1) = valee(1)*valdir(1)+valee(3)*valdir(2)
            valeq(2) = valee(3)*valdir(1)+valee(2)*valdir(2)
            valeq(3) = valee(4)*valdir(1)+valee(5)*valdir(2)
!
        else
!
            do i = 1, 6
                if (valee(i) .eq. r8vide()) valee(i) = 0.d0
            end do
            valeq(1) = valee(1)*valdir(1)+valee(3)*valdir(2)
            valeq(2) = valee(3)*valdir(1)+valee(2)*valdir(2)
            valeq(3) = valee(4)*valdir(1)+valee(6)*valdir(2)
            valeq(4) = valee(6)*valdir(1)+valee(5)*valdir(2)
!
        end if
!
    else if (codir .eq. 5) then
!
!     /* DIRECTION ACTIVE : X,Z */
!
        if (type .eq. 'V3') then
!
            do i = 1, 2
                if (valee(i) .eq. r8vide()) valee(i) = 0.d0
            end do
            valeq(1) = valee(1)*valdir(1)+valee(2)*valdir(3)
!
        else
!
            do i = 1, 5
                if (valee(i) .eq. r8vide()) valee(i) = 0.d0
            end do
            valeq(1) = valee(1)*valdir(1)+valee(4)*valdir(3)
            valeq(2) = valee(3)*valdir(1)+valee(5)*valdir(3)
            valeq(3) = valee(4)*valdir(1)+valee(2)*valdir(3)
!
        end if
!
    else if (codir .eq. 6) then
!
!     /* DIRECTION ACTIVE : Y,Z */
!
        if (type .eq. 'V3') then
!
            do i = 1, 2
                if (valee(i) .eq. r8vide()) valee(i) = 0.d0
            end do
            valeq(1) = valee(1)*valdir(2)+valee(2)*valdir(3)
!
        else
!
            do i = 1, 5
                if (valee(i) .eq. r8vide()) valee(i) = 0.d0
            end do
            valeq(1) = valee(3)*valdir(2)+valee(4)*valdir(3)
            valeq(2) = valee(1)*valdir(2)+valee(5)*valdir(3)
            valeq(3) = valee(5)*valdir(2)+valee(2)*valdir(3)
!
        end if
!
    else
!
!     /* DIRECTION ACTIVE : X,Y,Z */
!
        if (type .eq. 'V3') then
!
            do i = 1, 3
                if (valee(i) .eq. r8vide()) valee(i) = 0.d0
            end do
            valeq(1) = valee(1)*valdir(1)+valee(2)*valdir(2)+valee(3)*valdir(3)
!
        else
!
            do i = 1, 6
                if (valee(i) .eq. r8vide()) valee(i) = 0.d0
            end do
            valeq(1) = valee(1)*valdir(1)+valee(4)*valdir(2)+valee(5)*valdir(3)
            valeq(2) = valee(4)*valdir(1)+valee(2)*valdir(2)+valee(6)*valdir(3)
            valeq(3) = valee(5)*valdir(1)+valee(6)*valdir(2)+valee(3)*valdir(3)
!
        end if
!
    end if
!
end subroutine
