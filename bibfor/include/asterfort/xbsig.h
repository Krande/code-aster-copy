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
    subroutine xbsig(ndim, nnop, nfh, nfe,&
                     ddlc, ddlm, igeom, jpintt,&
                     cnset, heavt, lonch, basloc, sigma,&
                     nbsig, lsn, lst, ivectu,&
                     jpmilt, nfiss, jheavn, jstno, imate)
        integer(kind=8) :: nfiss
        integer(kind=8) :: nnop
        integer(kind=8) :: ndim
        integer(kind=8) :: nfh
        integer(kind=8) :: nfe
        integer(kind=8) :: ddlc
        integer(kind=8) :: ddlm
        integer(kind=8) :: igeom
        integer(kind=8) :: imate
        integer(kind=8) :: jpintt
        integer(kind=8) :: cnset(128)
        integer(kind=8) :: heavt(*)
        integer(kind=8) :: lonch(10)
        real(kind=8) :: basloc(*)
        real(kind=8) :: sigma(*)
        integer(kind=8) :: nbsig
        real(kind=8) :: lsn(nnop)
        real(kind=8) :: lst(nnop)
        integer(kind=8) :: ivectu
        integer(kind=8) :: jpmilt
        integer(kind=8) :: jheavn
        integer(kind=8) :: jstno
    end subroutine xbsig
end interface
