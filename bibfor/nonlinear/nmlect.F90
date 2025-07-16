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
subroutine nmlect(result, model, materField, mateco, caraElem, solver_)
!
    implicit none
!
#include "asterc/getres.h"
#include "asterfort/cresol.h"
#include "asterfort/medomm.h"
!
    character(len=*), intent(out) :: result
    character(len=*), intent(out) :: model
    character(len=*), intent(out) :: materField
    character(len=*), intent(out) :: mateco
    character(len=*), intent(out) :: caraElem
    character(len=*), optional, intent(out) :: solver_
!
! --------------------------------------------------------------------------------------------------
!
! Mechanics - Initializations
!
! Get parameters from command file
!
! --------------------------------------------------------------------------------------------------
!
! Out result           : name of results datastructure
! Out model            : name of model
! Out materField       : name of material
! Out mateco           : name of material characteristics (field)
! Out caraElem         : name of elementary characteristics (field)
! Out solver           : name of datastructure for solver
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16) :: k16dummy
!
! --------------------------------------------------------------------------------------------------
!

! - Get results
    call getres(result, k16dummy, k16dummy)

! - Get parameters from command file
    call medomm(model, materField, mateco, caraElem)

! - Get parameters for solver
    if (present(solver_)) then
        solver_ = '&&OP00XX.SOLVER'
        call cresol(solver_, 'V')
    end if
!
end subroutine
