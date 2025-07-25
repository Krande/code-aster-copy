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

function dpvpeq(x, n, const, fonc1, fonc2, &
                fonc3, fonc4)
    implicit none
    real(kind=8) :: x, n, const
    real(kind=8) :: fonc1, fonc2, fonc3, fonc4
    real(kind=8) :: dpvpeq, fonc, zero
! =====================================================================
! --- LOI DE COMPORTEMENT DE TYPE DRUCKER PRAGER VISCOPLASTIQUE -------
! ---- VISC_DRUC_PRAG -------------------------------------------------
! --- EQUATION NON LINEAIRE EN DP -------------------------------------
! =====================================================================
    parameter(zero=0.0d0)
! =====================================================================
    fonc = fonc1-fonc2*x-fonc3*x**2-fonc4*x**3
!
    if (fonc .lt. zero) then
        fonc = zero
    else
        fonc = fonc
    end if
    dpvpeq = const*fonc**n-x
! =====================================================================
end function
