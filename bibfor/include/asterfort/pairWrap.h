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
interface
    subroutine pairWrap(method, &
                        mesh, newgeo, mastConxInvName, &
                        mastNeighName, slavNeighName, &
                        pairTole, distRatio, verbosity, &
                        nbCellMast, listCellMast, &
                        nbCellSlav, listCellSlav, &
                        nbNodeMast, listNodeMast, &
                        nbPairZone, baseName)
        integer(kind=8), intent(in) :: method
        character(len=8), intent(in) :: mesh
        character(len=24), intent(in) :: newgeo, mastConxInvName
        character(len=24), intent(in) :: mastNeighName, slavNeighName
        real(kind=8), intent(in) :: pairTole, distRatio
        integer(kind=8), intent(in) :: verbosity
        integer(kind=8), intent(in) :: nbCellMast, listCellMast(nbCellMast)
        integer(kind=8), intent(in) :: nbCellSlav, listCellSlav(nbCellSlav)
        integer(kind=8), intent(in) :: nbNodeMast, listNodeMast(nbNodeMast)
        integer(kind=8), intent(out) :: nbPairZone
        character(len=8), intent(in) :: baseName
    end subroutine pairWrap
end interface
