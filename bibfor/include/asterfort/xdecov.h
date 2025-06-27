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
    subroutine xdecov(ndim, elp, nnop, nnose, it,&
                      pintt, cnset, heavt, ncomp, lsn,&
                      fisco, igeom, nfiss, ifiss, pinter,&
                      ninter, npts, ainter, nse, cnse,&
                      heav, nfisc, nsemax)
        integer(kind=8) :: nfisc
        integer(kind=8) :: nnop
        integer(kind=8) :: ndim
        character(len=8) :: elp
        integer(kind=8) :: nnose
        integer(kind=8) :: it
        real(kind=8) :: pintt(*)
        integer(kind=8) :: cnset(*)
        integer(kind=8) :: heavt(*)
        integer(kind=8) :: ncomp
        real(kind=8) :: lsn(*)
        integer(kind=8) :: fisco(*)
        integer(kind=8) :: igeom
        integer(kind=8) :: nfiss
        integer(kind=8) :: ifiss
        real(kind=8) :: pinter(*)
        integer(kind=8) :: ninter
        integer(kind=8) :: npts
        real(kind=8) :: ainter(*)
        integer(kind=8) :: nse
        integer(kind=8) :: cnse(6, 10)
        real(kind=8) :: heav(*)
        integer(kind=8) :: nsemax
    end subroutine xdecov
end interface
