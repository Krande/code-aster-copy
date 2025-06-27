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
    subroutine xxbsig(elrefp, elrese, ndim, coorse,&
                      igeom, he, nfh, ddlc, ddlm,&
                      nfe, basloc, nnop, npg, sigma,&
                      lsn, lst, nfiss,&
                      heavn, jstno, codopt, ivectu, imate)
        integer(kind=8) :: codopt
        integer(kind=8) :: nfiss
        integer(kind=8) :: npg
        integer(kind=8) :: nnop
        integer(kind=8) :: nfe
        integer(kind=8) :: nfh
        integer(kind=8), optional :: imate
        integer(kind=8) :: ndim
        character(len=8) :: elrefp
        character(len=8) :: elrese
        real(kind=8) :: coorse(*)
        integer(kind=8) :: igeom
        real(kind=8) :: he(nfiss)
        integer(kind=8) :: ddlc
        integer(kind=8) :: ddlm
        real(kind=8) :: basloc(3*ndim*nnop)
        real(kind=8) :: sigma(codopt*(2*ndim-1)+1, codopt*(npg-1)+1)
        real(kind=8) :: lsn(nnop)
        real(kind=8) :: lst(nnop)
        integer(kind=8) :: heavn(nnop, 5)
        integer(kind=8) :: ivectu
        integer(kind=8) :: jstno
    end subroutine xxbsig
end interface
