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
subroutine lkcalg(dfdsig, vecn, g, devgii)
!
    implicit none
#include "asterfort/lcdevi.h"
    real(kind=8) :: dfdsig(6), vecn(6), g(6), devgii
! --- MODELE LETK : LAIGLE VISCOPLASTIQUE--------------------------
! =================================================================
! --- BUT : CALCUL DE G=df/dsig-(df/dsig*n)*n----------------------
! =================================================================
! IN  : DFDSIG : df/dsig ------------------------------------------
!        VECN  : --------------------------------------------------
! OUT : G : G=df/dsig-(df/dsig*n)*n -------------------------------
!      DEVGII : SECOND INVARIANT DE G -----------------------------
! =================================================================
    common/tdim/ndt, ndi
    integer(kind=8) :: ndi, ndt, i
    real(kind=8) :: devg(6), fact1
!
! =================================================================
! --- CALCUL DE G -------------------------------------------------
! =================================================================
    g = 0.d0
!
    fact1 = dot_product(dfdsig(1:ndt), vecn(1:ndt))
!
    do i = 1, ndt
        g(i) = dfdsig(i)-fact1*vecn(i)
    end do
! =================================================================
! --- CALCUL DU DEVIATEUR DE G ET DE SA NORME ---------------------
! =================================================================
    call lcdevi(g, devg)
    devgii = norm2(devg(1:ndt))
! =================================================================
end subroutine
