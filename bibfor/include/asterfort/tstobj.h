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
    subroutine tstobj(ob, perm, resume, sommi, sommr,&
                      lonuti, lonmax, type, iret, ni)
        character(len=*) :: ob
        character(len=*) :: perm
        integer(kind=8) :: resume
        integer(kind=8) :: sommi
        real(kind=8) :: sommr
        integer(kind=8) :: lonuti
        integer(kind=8) :: lonmax
        character(len=3) :: type
        integer(kind=8) :: iret
        integer(kind=8) :: ni
    end subroutine tstobj
end interface
