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
    subroutine lcloca(coeft, nmat, nbcomm,&
                      nphas, sigi, vini, iphas, granb,&
                      loca, sigg)
        integer(kind=8) :: nmat
        real(kind=8) :: coeft(nmat)
        integer(kind=8) :: nbcomm(nmat, 3)
        integer(kind=8) :: nphas
        real(kind=8) :: sigi(6)
        real(kind=8) :: vini(*)
        integer(kind=8) :: iphas
        real(kind=8) :: granb(6)
        character(len=16) :: loca
        real(kind=8) :: sigg(6)
    end subroutine lcloca
end interface
