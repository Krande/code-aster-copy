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
    subroutine stkmai(ifl, icl, iv, rv, cv,&
                      cnl, mcl, nbm, nume, numn,&
                      cnx, typ, fmt, irteti, nommai)
        integer(kind=8) :: nbm
        integer(kind=8) :: ifl
        integer(kind=8) :: icl
        integer(kind=8) :: iv
        real(kind=8) :: rv
        character(len=*) :: cv
        character(len=14) :: cnl
        character(len=8) :: mcl(nbm)
        integer(kind=8) :: nume
        integer(kind=8) :: numn
        character(len=24) :: cnx
        character(len=24) :: typ
        integer(kind=8) :: fmt(nbm)
        integer(kind=8) :: irteti
        character(len=24) :: nommai
    end subroutine stkmai
end interface
