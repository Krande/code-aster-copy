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
    subroutine xside3(elrefp, ndim, coorse, elrese, igeom,&
                      he, nfh, ddlc, ddlm, nfe,&
                      basloc, nnop, npg, idecpg, imate,&
                      idepl, lsn, lst, nfiss,&
                      heavn, jstno, sig)
        integer(kind=8) :: nfiss
        integer(kind=8) :: npg
        integer(kind=8) :: nnop
        integer(kind=8) :: ndim
        character(len=8) :: elrefp
        real(kind=8) :: coorse(*)
        character(len=8) :: elrese
        integer(kind=8) :: igeom
        real(kind=8) :: he(nfiss)
        integer(kind=8) :: nfh
        integer(kind=8) :: ddlc
        integer(kind=8) :: ddlm
        integer(kind=8) :: nfe
        real(kind=8) :: basloc(9*nnop)
        integer(kind=8) :: idecpg
        integer(kind=8) :: imate
        integer(kind=8) :: idepl
        real(kind=8) :: lsn(nnop)
        real(kind=8) :: lst(nnop)
        integer(kind=8) :: heavn(nnop, 5)
        integer(kind=8) :: jstno
        real(kind=8) :: sig(6, npg)
    end subroutine xside3
end interface
