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
    subroutine fgvdmg(nomsym, nomsd, nommat, nomnap, nomfon,&
                      mexpic, mcompt, mdomag, nbord, nbpt,&
                      ntcmp, nbcmp, numcmp, impr, vdomag)
        character(len=16) :: nomsym
        character(len=19) :: nomsd
        character(len=8) :: nommat
        character(len=8) :: nomnap
        character(len=8) :: nomfon
        character(len=*) :: mexpic
        character(len=*) :: mcompt
        character(len=*) :: mdomag
        integer(kind=8) :: nbord
        integer(kind=8) :: nbpt
        integer(kind=8) :: ntcmp
        integer(kind=8) :: nbcmp
        integer(kind=8) :: numcmp(*)
        integer(kind=8) :: impr
        real(kind=8) :: vdomag(*)
    end subroutine fgvdmg
end interface
