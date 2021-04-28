! --------------------------------------------------------------------
! Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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
interface
    subroutine cmeGetParameters(option       ,&
                                model        , caraElem    ,&
                                mate         , mateco      , comporMult,&
                                listLoadK8   , nbLoad      ,&
                                rigiMatrElem , massMatrElem,&
                                timeCurr     , timeIncr    , modeFourier,&
                                sigm         , strx        , disp,&
                                calcElemModel)
        character(len=16), intent(out) :: option
        character(len=8), intent(out) :: model, caraElem
        character(len=24), intent(out) :: mate, mateco, comporMult
        character(len=8), pointer :: listLoadK8(:)
        integer, intent(out) :: nbLoad
        character(len=19), intent(out) :: rigiMatrElem, massMatrElem
        real(kind=8), intent(out) :: timeCurr, timeIncr
        integer, intent(out) :: modeFourier
        character(len=8), intent(out) :: sigm, strx, disp
        character(len=8), intent(out) :: calcElemModel
    end subroutine cmeGetParameters
end interface
