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
    subroutine xdecqv(nnose, it, cnset, heavt, lsn, igeom,&
                      ninter, npts, ndim, ainter, nse, cnse,&
                      heav, nsemax, pinter, pmilie, pintt, pmitt, cut,&
                      ncomp, nfisc, nfiss, ifiss, elp, fisco,&
                      lonref, txlsn, tx)
        integer(kind=8) :: nnose
        integer(kind=8) :: it
        integer(kind=8) :: cnset(*)
        real(kind=8) :: lsn(*)
        integer(kind=8) :: igeom
        integer(kind=8) :: ninter
        integer(kind=8) :: npts
        integer(kind=8) :: ndim
        real(kind=8) :: ainter(*)
        real(kind=8) :: pinter(*)
        real(kind=8) :: pmilie(*)
        real(kind=8) :: pintt(*)
        real(kind=8) :: pmitt(*)
        integer(kind=8) :: nse
        integer(kind=8) :: cnse(6, 10)
        real(kind=8) :: heav(*)
        integer(kind=8) :: nsemax
        aster_logical :: cut
        integer(kind=8) :: heavt(*)
        integer(kind=8) :: ncomp
        integer(kind=8) :: nfisc
        integer(kind=8) :: nfiss
        integer(kind=8) :: ifiss
        integer(kind=8) :: fisco(*)
        character(len=8) :: elp
        real(kind=8) :: lonref
        real(kind=8) :: txlsn(28)
        real(kind=8) :: tx(3, 7)
    end subroutine xdecqv
end interface
