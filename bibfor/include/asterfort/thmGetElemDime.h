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
#include "asterf_types.h"
!
interface
    subroutine thmGetElemDime(ndim     , nnos   , nnom   ,&
                              mecani   , press1 , press2 , tempe , second, &
                              nddls    , nddlm  , &
                              nddl_meca, nddl_p1, nddl_p2, nddl_2nd, &
                              dimdep   , dimdef , dimcon , dimuel)
        integer(kind=8), intent(in) :: ndim, nnos, nnom
        integer(kind=8), intent(in) :: mecani(5), press1(7), press2(7), tempe(5), second(5)
        integer(kind=8), intent(out) :: nddls, nddlm
        integer(kind=8), intent(out) :: nddl_meca, nddl_p1, nddl_p2, nddl_2nd
        integer(kind=8), intent(out) :: dimdep, dimdef, dimcon, dimuel
    end subroutine thmGetElemDime
end interface
