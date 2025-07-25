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
    subroutine rcsp01(nbm, adrm, ipt, sp3, sp4,&
                      sp5, alphaa, alphab, nbth, iocs,&
                      sp6)
        integer(kind=8) :: nbm
        integer(kind=8) :: adrm(*)
        integer(kind=8) :: ipt
        real(kind=8) :: sp3
        real(kind=8) :: sp4
        real(kind=8) :: sp5
        real(kind=8) :: alphaa
        real(kind=8) :: alphab
        integer(kind=8) :: nbth
        integer(kind=8) :: iocs
        real(kind=8) :: sp6
    end subroutine rcsp01
end interface
