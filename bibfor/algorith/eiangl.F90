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
subroutine eiangl(ndim, nno2, angnau, ang)
    implicit none
#include "asterc/r8dgrd.h"
    integer(kind=8) :: ndim, nno2
    real(kind=8) :: angnau(3), ang(merge(1, 3, ndim .eq. 2), nno2)
!
!--------------------------------------------------
!  DEFINITION DES ANGLES NAUTIQUES AUX NOEUDS
!  EN RADIAN POUR L'ELEMENT D'INTERFACE
!
!  IN  : NDIM,NNO2
!        ANGNAU : ANGLES NAUTIQUES EN DEGRES
!  OUT :
!        ANG : ANGLES NAUTIQUES AUX NOEUDS EN RADIAN
!--------------------------------------------------
    integer(kind=8)::i
    real(kind=8):: deg_2_rad
!--------------------------------------------------
    deg_2_rad = r8dgrd()
    forall (i=1:size(ang, 1)) ang(i, :) = angnau(i)*deg_2_rad
end subroutine
