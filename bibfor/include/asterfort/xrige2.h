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
    subroutine xrige2(elrefp, elrese, ndim, coorse, igeom,&
                      he, heavn, ddlh, ddlc, nfe,&
                      basloc, nnop, npg, lsn, lst,&
                      sig, matuu, jstno, imate)
        integer(kind=8) :: nnop
        integer(kind=8) :: ndim
        integer(kind=8) :: heavn(27,5)
        character(len=8) :: elrefp
        character(len=8) :: elrese
        real(kind=8) :: coorse(*)
        integer(kind=8) :: igeom
        real(kind=8) :: he
        integer(kind=8) :: ddlh
        integer(kind=8) :: ddlc
        integer(kind=8) :: nfe
        integer(kind=8) :: jstno
        integer(kind=8) :: imate
        real(kind=8) :: basloc(6*nnop)
        integer(kind=8) :: npg
        real(kind=8) :: lsn(nnop)
        real(kind=8) :: lst(nnop)
        real(kind=8) :: sig(48)
        real(kind=8) :: matuu(*)
    end subroutine xrige2
end interface
