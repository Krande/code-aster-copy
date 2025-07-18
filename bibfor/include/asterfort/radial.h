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
    subroutine radial(nbsig, sigm, sigp, indm, indp,&
                      icine, xm, xp, normdn)
        integer(kind=8) :: nbsig
        real(kind=8) :: sigm(nbsig)
        real(kind=8) :: sigp(nbsig)
        real(kind=8) :: indm
        real(kind=8) :: indp
        integer(kind=8) :: icine
        real(kind=8) :: xm(6)
        real(kind=8) :: xp(6)
        real(kind=8) :: normdn
    end subroutine radial
end interface
