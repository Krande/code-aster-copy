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
    subroutine cmePrep(optionz      , modelz       ,&
                       timeCurr     , timeIncr     , chtime     ,&
                       nbLoad       , listLoadK8   , listLoadK24,&
                       calcElemModel, onlyDirichlet,&
                       matrElemz    , listElemCalc)
        character(len=*), intent(in) :: optionz, modelz
        real(kind=8), intent(in) :: timeCurr, timeIncr
        character(len=24), intent(out) :: chtime
        integer, intent(in) :: nbLoad
        character(len=8), pointer :: listLoadK8(:)
        character(len=24), pointer :: listLoadK24(:)
        character(len=8), intent(in) :: calcElemModel
        aster_logical, intent(out) :: onlyDirichlet
        character(len=*), intent(in) :: matrElemz
        character(len=24), intent(out) :: listElemCalc
    end subroutine cmePrep
end interface
