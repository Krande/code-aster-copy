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
    subroutine pbflu0(rhof, hmoy, rmoy, long, icoq,&
                      imod, nbm, rkip, tcoef, d)
        integer(kind=8) :: nbm
        real(kind=8) :: rhof
        real(kind=8) :: hmoy
        real(kind=8) :: rmoy
        real(kind=8) :: long
        integer(kind=8) :: icoq
        integer(kind=8) :: imod
        real(kind=8) :: rkip
        real(kind=8) :: tcoef(10, nbm)
        real(kind=8) :: d(6)
    end subroutine pbflu0
end interface
