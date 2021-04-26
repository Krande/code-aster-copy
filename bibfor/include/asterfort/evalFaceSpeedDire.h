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
    subroutine evalFaceSpeedDire(cellDime, jvLoad , speedDire,&
                                 ipg     , nx     , ny       ,&
                                 lFunc_  , lReal_ , lCplx_   ,&
                                 lTime_  , time_  ,&
                                 x_      , y_     ,&
                                 z_      , nz_)
        integer, intent(in) :: cellDime, jvLoad
        real(kind=8), intent(out) :: speedDire
        integer, intent(in) :: ipg
        real(kind=8), intent(in) :: nx, ny
        aster_logical, optional, intent(in) :: lFunc_, lReal_, lCplx_, lTime_
        real(kind=8), optional, intent(in) :: time_, x_, y_
        real(kind=8), optional, intent(in) :: z_, nz_
    end subroutine evalFaceSpeedDire
end interface
