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
    subroutine sleelt(iunv, maxnod, nbtyma, indic, permut,&
                      nbmail, mint, mant, datset, inum)
        integer(kind=8) :: nbtyma
        integer(kind=8) :: maxnod
        integer(kind=8) :: iunv
        integer(kind=8) :: indic(nbtyma)
        integer(kind=8) :: permut(maxnod, nbtyma)
        integer(kind=8) :: nbmail(nbtyma)
        integer(kind=8) :: mint(nbtyma)
        integer(kind=8) :: mant(nbtyma)
        integer(kind=8) :: datset
        integer(kind=8) :: inum
    end subroutine sleelt
end interface
