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
    subroutine xpesro(elrefp, ndim, coorse, igeom, jheavt, ncomp,&
                      heavn, nfh, ddlc, nfe, nfiss,&
                      ise, nnop, jlsn, jlst, ivectu,&
                      fno, imate, jbaslo, jstno)
        integer(kind=8) :: nnop
        integer(kind=8) :: ndim
        character(len=8) :: elrefp
        real(kind=8) :: coorse(*)
        integer(kind=8) :: igeom
        integer(kind=8) :: jheavt
        integer(kind=8) :: ncomp
        integer(kind=8) :: heavn(27,5)
        integer(kind=8) :: nfh
        integer(kind=8) :: ddlc
        integer(kind=8) :: nfe
        integer(kind=8) :: nfiss
        integer(kind=8) :: ise
        integer(kind=8) :: jlsn
        integer(kind=8) :: jlst
        integer(kind=8) :: imate
        integer(kind=8) :: jbaslo
        integer(kind=8) :: jstno
        integer(kind=8) :: ivectu
        real(kind=8) :: fno(ndim*nnop)
    end subroutine xpesro
end interface
