! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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

subroutine xminte(ndim, integ, fpg)
!
    implicit none
    integer :: ndim, integ
    character(len=8) :: fpg
!
!
! ROUTINE CONTACT (METHODE XFEM HPP - CALCUL ELEM.)
!
! --- SCHEMA D'INTEGRATION NUMERIQUE SUR LA SURFACE DE CONTACT
!
! ----------------------------------------------------------------------
!
! IN  NDIM   : DIMENSION DU PROBLEME
! IN  INTEG  : NUMERO DU TYPE D'INTEGRATION
! OUT FPG    : NOM (LOCAL) DE LA FAMILLE DE POINTS DE GAUSS
!
! ----------------------------------------------------------------------
!
    if (ndim .eq. 3) then
        if (integ .eq. 1) fpg = 'NOEU'
        if (integ .eq. 62 .or. integ .eq. 72 .or. integ .eq. 82 .or. integ .eq. 92 .or. integ &
            .eq. 102) fpg = 'GAUSS'
        if (integ .eq. 13) fpg = 'SIMP'
        if (integ .eq. 23 .or. integ .eq. 33 .or. integ .eq. 43) fpg = 'SIMP1'
        if (integ .eq. 34 .or. integ .eq. 44 .or. integ .eq. 54) fpg = 'COTES'
        if (integ .eq. 12) fpg = 'FPG3'
        if (integ .eq. 22) fpg = 'FPG3'
        if (integ .eq. 32) fpg = 'FPG4'
        if (integ .eq. 42) fpg = 'FPG6'
        if (integ .eq. 52) fpg = 'FPG7'
!
    else if (ndim .eq. 2) then
        if (integ .eq. 1) fpg = 'NOEU'
        if (integ .eq. 32) fpg = 'GAUSS'
        if (integ .eq. 13) fpg = 'SIMP'
        if (integ .eq. 23 .or. integ .eq. 33 .or. integ .eq. 43) fpg = 'SIMP1'
        if (integ .eq. 34) fpg = 'COTES'
        if (integ .eq. 54 .or. integ .eq. 44) fpg = 'COTES1'
        if (integ .eq. 84 .or. integ .eq. 64 .or. integ .eq. 74 .or. integ .eq. 94 .or. integ &
            .eq. 104) fpg = 'COTES2'
        if (integ .eq. 12) fpg = 'FPG2'
        if (integ .eq. 22) fpg = 'FPG2'
        if (integ .eq. 32) fpg = 'FPG3'
        if (integ .eq. 42 .or. integ .eq. 52 .or. integ .eq. 62 .or. integ .eq. 72 .or. integ &
            .eq. 82 .or. integ .eq. 92 .or. integ .eq. 102) fpg = 'FPG4'
!
    end if
!
end subroutine
