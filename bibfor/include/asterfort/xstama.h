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
    subroutine xstama(noma, nbma, nmafis, jmafis,&
                      ncouch, lisnoe, stano, cnslt, cnsln,&
                      jmafon, jmaen1, jmaen2, jmaen3, nmafon,&
                      nmaen1, nmaen2, nmaen3, typdis)
        character(len=8) :: noma
        integer(kind=8) :: nbma
        integer(kind=8) :: nmafis
        integer(kind=8) :: jmafis
        integer(kind=8) :: ncouch
        character(len=24) :: lisnoe
        integer(kind=8) :: stano(*)
        character(len=19) :: cnslt
        character(len=19) :: cnsln
        integer(kind=8) :: jmafon
        integer(kind=8) :: jmaen1
        integer(kind=8) :: jmaen2
        integer(kind=8) :: jmaen3
        integer(kind=8) :: nmafon
        integer(kind=8) :: nmaen1
        integer(kind=8) :: nmaen2
        integer(kind=8) :: nmaen3
        character(len=16) :: typdis
    end subroutine xstama
end interface
