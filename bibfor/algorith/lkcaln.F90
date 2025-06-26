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
subroutine lkcaln(s, b, vecn, retcom)
!
    implicit none
#include "asterc/r8miem.h"
    integer(kind=8) :: retcom
    real(kind=8) :: b, s(6), vecn(6)
! --- MODELE LETK : LAIGLE VISCOPLASTIQUE--------------------------
! =================================================================
! --- BUT : CALCUL DE N -------------------------------------------
! =================================================================
! IN  : S      : DEVIATEUR DES CONTRAINTES ------------------------
! --- : B      : PARAMETRE DU CALCUL DE LA NORMALE ----------------
! OUT : VECN   : N = (B*S/SII-I)/SQRT(B**2+3) ---------------------
! =================================================================
    integer(kind=8) :: i, ndt, ndi
    real(kind=8) :: sii, racine, un, trois, kron(6), zero
! =================================================================
! --- INITIALISATION DE PARAMETRE ---------------------------------
! =================================================================
    parameter(un=1.0d0)
    parameter(trois=3.0d0)
    parameter(zero=0.0d0)
! =================================================================
    common/tdim/ndt, ndi
! =================================================================
    data kron/un, un, un, zero, zero, zero/
! --- INITIALISATION ----------------------------------------------
! =================================================================
    retcom = 1
    vecn = 0.d0
! =================================================================
! --- CALCUL DE N -------------------------------------------------
! =================================================================

    sii = norm2(s(1:ndt))
    if (sii .gt. r8miem()) then
        racine = sqrt(b*b+trois)
        do i = 1, ndt
            vecn(i) = (b*s(i)/sii-kron(i))/racine
        end do
        retcom = 0
    end if

end subroutine lkcaln
