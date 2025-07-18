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

subroutine ratanp(dpstrs, rprops, ii, jj, mm, edge, apex)
!***********************************************************************
!
! OBJECT:
!
! COMPUTATION OF CONSISTENT TANGENT MODULUS FOR RANKINE TYPE
! ELASTO-PLASTIC MATERIAL WITH ASSOCIATIVE/NON-ASSOCIATIVE FLOW RULE
!
! ----------------------------------------------------------------------
!
!     LOI DE COMPORTEMENT DE RANKINE
!
! IN  PSTRS   : CONTRAINTES PRINCIPALES
! IN  II      : COMPOSANTE DE LA CONTRAINTE PRINCIPALE MINEURE
! IN  JJ      : COMPOSANTE DE LA CONTRAINTE PRINCIPALE MAJEURE
! IN  MM      : COMPOSANTE DE LA CONTRAINTE PRINCIPALE INTERMEDIAIRE
! IN  EDGE    : Y-A-T-IL DEUX  MECANISMES ACTIFS
! IN  APEX    : Y-A-T-IL TROIS MECANISMES ACTIFS
!
! OUT DPSTRS  : MATRICE TANGENTE DANS L'ESPACE DES DIRECTIONS PRINCIPALE
!
!***********************************************************************
    implicit none
! ======================================================================
!
    real(kind=8) :: dpstrs(3, 3)
    real(kind=8) :: rprops(*)
    real(kind=8) :: edge
    real(kind=8) :: apex
!
#include "asterf_types.h"
!
! Declaration of integer type variables
    integer(kind=8) :: ii, jj, mm
!
!     aster_logical :: epflag, epflag0
! Real arrays and variables
    real(kind=8) :: young, poiss, constb, gmodu, bulk
    real(kind=8) :: r2g, r4g, r2bulk, r1d3, r2d3, r2gd3, r4gd3
    real(kind=8) :: consta, denom, b1, b2
    real(kind=8) :: r1ddet, r1, r2, r3, r4, r0
!
    data r0, r1, r2, r3, r4/&
     &    0.0d0, 1.0d0, 2.0d0, 3.0d0, 4.0d0/
!
!
! Set material properties
    young = rprops(1)
    poiss = rprops(2)
!
! Set some constants
    gmodu = young/(r2*(r1+poiss))
    bulk = young/(r3*(r1-r2*poiss))
    r2g = r2*gmodu
    r4g = r4*gmodu
    r2bulk = r2*bulk
    r1d3 = r1/r3
    r2d3 = r2*r1d3
    r2gd3 = r2g*r1d3
    r4gd3 = r4g*r1d3
!
    consta = bulk+r4gd3
    constb = bulk-r2gd3
!
    dpstrs(:, :) = r0
!
! Compute elastoplastic consistent tangent
! ----------------------------------------
    if (edge .eq. r1) then
! -------------------------------------------------------------------------
!
! Tangent consistent with 2-vector return to edge
!
! -------------------------------------------------------------------------
        r1ddet = r3*gmodu*bulk/(r1d3*gmodu+bulk)
        dpstrs(jj, jj) = r1ddet
!
    else if (apex .eq. r0) then
! -------------------------------------------------------------------------
!
! Tangent consistent with 1-vector return to main active plane
!
! -------------------------------------------------------------------------
        denom = consta
!
        b1 = (consta*consta-constb*constb)/denom
        b2 = constb*(consta-constb)/denom
!
        dpstrs(mm, mm) = b1
        dpstrs(mm, jj) = b2
        dpstrs(jj, mm) = b2
        dpstrs(jj, jj) = b1
    end if
!
end subroutine
