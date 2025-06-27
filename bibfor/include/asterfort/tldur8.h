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
    subroutine tldur8(nommat, hcol, adia, ablo, npivot,&
                      neq, nbbloc, ildeb, ilfin, eps)
        character(len=*) :: nommat
        integer(kind=8) :: hcol(*)
        integer(kind=8) :: adia(*)
        integer(kind=8) :: ablo(*)
        integer(kind=8) :: npivot
        integer(kind=8) :: neq
        integer(kind=8) :: nbbloc
        integer(kind=8) :: ildeb
        integer(kind=8) :: ilfin
        real(kind=8) :: eps
    end subroutine tldur8
end interface
