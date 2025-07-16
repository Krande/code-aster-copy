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
    subroutine xdecqu(nnose, it, ndim, cnset, jlsn,&
                      igeom, pinter, ninter, npts, ainter,&
                      pmilie, nmilie, mfis, tx, txlsn,&
                      pintt, pmitt, ifiss, nfiss, fisco,&
                      nfisc, cut, coupe, exit, joncno, condition_joncno)
        integer(kind=8) :: ndim
        integer(kind=8) :: nnose
        integer(kind=8) :: it
        integer(kind=8) :: cnset(*)
        integer(kind=8) :: jlsn
        integer(kind=8) :: igeom
        real(kind=8) :: pinter(*)
        integer(kind=8) :: ninter
        integer(kind=8) :: npts
        real(kind=8) :: ainter(*)
        real(kind=8) :: pmilie(*)
        integer(kind=8) :: nmilie
        integer(kind=8) :: mfis
        real(kind=8) :: tx(3, 7)
        real(kind=8) :: txlsn(28)
        real(kind=8) :: pintt(*)
        real(kind=8) :: pmitt(*)
        integer(kind=8) :: ifiss
        integer(kind=8) :: nfiss
        integer(kind=8) :: fisco(*)
        integer(kind=8) :: nfisc
        aster_logical :: cut
        integer(kind=8) :: coupe(nfiss)
        integer(kind=8) :: exit(2)
        integer(kind=8) :: joncno
        aster_logical :: condition_joncno
    end subroutine xdecqu
end interface
