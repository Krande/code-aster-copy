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
    subroutine maskmn(nbcmp, nbno, nbec, mcoddl, imask,&
                      numord, nbdef)
        integer(kind=8) :: nbec
        integer(kind=8) :: nbno
        integer(kind=8) :: nbcmp
        integer(kind=8) :: mcoddl(nbno*nbec, 2)
        integer(kind=8) :: imask(nbno*nbec)
        integer(kind=8) :: numord(nbno)
        integer(kind=8) :: nbdef
    end subroutine maskmn
end interface
