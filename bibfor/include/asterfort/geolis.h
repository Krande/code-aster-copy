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
!
interface
    subroutine geolis(modgen, sst1, sst2, intf1, intf2,&
                      geom1, geom2, limail, nmga1)
        character(len=8) :: modgen
        character(len=8) :: sst1
        character(len=8) :: sst2
        character(len=8) :: intf1
        character(len=8) :: intf2
        character(len=24) :: geom1
        character(len=24) :: geom2
        character(len=*) :: limail
        integer(kind=8) :: nmga1
    end subroutine geolis
end interface
