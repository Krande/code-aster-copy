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
    subroutine cplch(icompo, cpiter, ti, tf, numpas,&
                     nomvar, idim, taille, nompal, info)
        integer(kind=8) :: icompo
        integer(kind=4) :: cpiter
        real(kind=4) :: ti
        real(kind=4) :: tf
        integer(kind=4) :: numpas
        character(len=144) :: nomvar
        integer(kind=4) :: idim
        integer(kind=4) :: taille
        character(len=6) :: nompal
        integer(kind=4) :: info
    end subroutine cplch
end interface
