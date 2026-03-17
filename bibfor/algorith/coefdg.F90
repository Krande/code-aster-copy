! --------------------------------------------------------------------
! Copyright (C) 1991 - 2026 - EDF - www.code-aster.org
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
subroutine coefdg(relaComp, materPara, dpida2)
!
    use MaterialPara_type
    implicit none
!
#include "asterfort/rcvalb.h"
!
    character(len=16), intent(in) :: relaComp
    type(Material_Para), intent(inout) :: materPara
    real(kind=8), intent(out) :: dpida2
!
! --------------------------------------------------------------------------------------------------
!
!     LOIS A GRADIENTS : COEFFICIENT DIAGONAL MATRICE GVNO
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8), parameter :: poum = "+"
    integer(kind=8), parameter :: nbProp = 2
    character(len=8), parameter :: propName(nbProp) = (/"E ", "NU"/)
    real(kind=8) :: propVale(nbProp)
    integer(kind=8) :: propCode(nbProp)
!
! --------------------------------------------------------------------------------------------------
!
    dpida2 = 0.d0
!
    if (relaComp .eq. 'ENDO_CARRE') then
        call rcvalb(materPara%schemePara%fami, &
                    materPara%schemePara%kpg, &
                    materPara%schemePara%ksp, &
                    poum, &
                    materPara%jvMaterCode, &
                    ' ', 'ELAS', &
                    0, ' ', [0.d0], &
                    nbProp, propName, propVale, propCode, 2)
        dpida2 = propVale(1)
    end if
!
end subroutine
