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
    subroutine xsifle(ndim, ifa, jptint, cface,&
                      igeom, nfh, jheavn, singu, nfe, ddlc,&
                      ddlm, jlsn, jlst, jstno, ipres, ipref, itemps,&
                      idepl, nnop, valres, basloc, ithet,&
                      nompar, option, igthet, jbasec,&
                      contac)
        integer(kind=8) :: nnop
        integer(kind=8) :: ndim
        integer(kind=8) :: ifa
        integer(kind=8) :: jptint
        integer(kind=8) :: cface(30, 6)
        integer(kind=8) :: igeom
        integer(kind=8) :: nfh
        integer(kind=8) :: singu
        integer(kind=8) :: nfe
        integer(kind=8) :: ddlc
        integer(kind=8) :: ddlm
        integer(kind=8) :: jlst
        integer(kind=8) :: jheavn
        integer(kind=8) :: ipres
        integer(kind=8) :: ipref
        integer(kind=8) :: itemps
        integer(kind=8) :: idepl
        integer(kind=8) :: jstno
        integer(kind=8) :: jlsn
        real(kind=8) :: valres(3)
        real(kind=8) :: basloc(9*nnop)
        integer(kind=8) :: ithet
        character(len=8) :: nompar(4)
        character(len=16) :: option
        integer(kind=8) :: igthet
        integer(kind=8) :: jbasec
        integer(kind=8) :: contac
    end subroutine xsifle
end interface
