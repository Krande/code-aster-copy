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
interface 
    subroutine xsifl1(elrefp, angl, basloc, coeff, coeff3, ddlm,&
                      ddls, dfdi, ff, he, heavn, idepl,&
                      igthet, ipref, ipres, ithet, jac,&
                      jlsn, jlst, jstno, ka, mu, nd,&
                      ndim, nfh, nnop, nnops, itemps,&
                      nompar, option, singu, xg, igeom)
        character(len=8) :: elrefp
        integer(kind=8) :: nnop
        integer(kind=8) :: ndim
        real(kind=8) :: angl(2)
        real(kind=8) :: basloc(9*nnop)
        real(kind=8) :: coeff
        real(kind=8) :: coeff3
        integer(kind=8) :: ddlm
        integer(kind=8) :: ddls
        integer(kind=8) :: heavn(nnop,5)
        real(kind=8) :: dfdi(nnop, ndim)
        real(kind=8) :: ff(27)
        real(kind=8) :: he(2)
        integer(kind=8) :: idepl
        integer(kind=8) :: igthet
        integer(kind=8) :: ipref
        integer(kind=8) :: ipres
        integer(kind=8) :: ithet
        real(kind=8) :: jac
        integer(kind=8) :: jlst
        integer(kind=8) :: jlsn
        integer(kind=8) :: jstno
        real(kind=8) :: ka
        real(kind=8) :: mu
        real(kind=8) :: nd(3)
        integer(kind=8) :: nfh
        integer(kind=8) :: nnops
        integer(kind=8) :: itemps
        character(len=8) :: nompar(4)
        character(len=16) :: option
        integer(kind=8) :: singu
        real(kind=8) :: xg(3)
        integer(kind=8) :: igeom
    end subroutine xsifl1
end interface 
