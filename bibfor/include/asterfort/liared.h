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
    subroutine liared(nomres, fmli, iblo, liamod, nlilia,&
                      ncolia, promod, nlipro, ncopro, taille,&
                      indcol, nbcol)
        character(len=8) :: nomres
        character(len=24) :: fmli
        integer(kind=8) :: iblo
        character(len=24) :: liamod
        integer(kind=8) :: nlilia
        integer(kind=8) :: ncolia
        character(len=24) :: promod
        integer(kind=8) :: nlipro
        integer(kind=8) :: ncopro
        integer(kind=8) :: taille(2)
        character(len=24) :: indcol
        integer(kind=8) :: nbcol
    end subroutine liared
end interface
