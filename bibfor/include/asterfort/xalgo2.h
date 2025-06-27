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
    subroutine xalgo2(ndim, elrefp, it, nnose, cnset, typma, ndime,&
                      geom, lsnelp, pmilie, ninter, ainter, ar, npts, nptm, &
                      pmmax, nmilie, mfis, lonref, pinref, pintt, pmitt, jonc, exit)
        integer(kind=8) :: ndim
        integer(kind=8) :: it
        integer(kind=8) :: nnose
        integer(kind=8) :: cnset(*)
        integer(kind=8) :: ndime
        real(kind=8) :: geom(81)
        integer(kind=8) :: ninter
        integer(kind=8) ::  ar(12, 3)
        integer(kind=8) :: npts
        integer(kind=8) :: nptm
        integer(kind=8) :: nbar
        integer(kind=8) :: pmmax
        integer(kind=8) :: nmilie
        integer(kind=8) :: mfis
        character(len=8) :: typma
        character(len=8) :: elrefp
        real(kind=8) :: lonref
        real(kind=8) :: ainter(*)
        real(kind=8) :: pmilie(*)
        real(kind=8) :: pinref(*) 
        real(kind=8) :: lsnelp(27)
        real(kind=8) :: pintt(*)
        real(kind=8) :: pmitt(*)
        aster_logical :: jonc
        integer(kind=8) :: exit(2)
    end subroutine xalgo2
end interface 
