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
    subroutine poibij(npv, vabs, geom, fsvr, nbm,&
                      vicoq, torco, tcoef, freq, imasse,&
                      maj, vecpr)
        integer(kind=8) :: nbm
        integer(kind=8) :: npv
        real(kind=8) :: vabs(npv)
        real(kind=8) :: geom(9)
        real(kind=8) :: fsvr(7)
        integer(kind=8) :: vicoq(nbm)
        real(kind=8) :: torco(4, nbm)
        real(kind=8) :: tcoef(10, nbm)
        real(kind=8) :: freq(2*nbm*npv)
        integer(kind=8) :: imasse
        real(kind=8) :: maj(nbm)
        real(kind=8) :: vecpr(*)
    end subroutine poibij
end interface
