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
    subroutine meobl1(eps, b, d, deltab, deltad,&
                      mult, lambda, mu, ecrob, ecrod,&
                      alpha, k1, k2, bdim, dsidep)
        real(kind=8) :: eps(6)
        real(kind=8) :: b(6)
        real(kind=8) :: d
        real(kind=8) :: deltab(6)
        real(kind=8) :: deltad
        real(kind=8) :: mult
        real(kind=8) :: lambda
        real(kind=8) :: mu
        real(kind=8) :: ecrob
        real(kind=8) :: ecrod
        real(kind=8) :: alpha
        real(kind=8) :: k1
        real(kind=8) :: k2
        integer(kind=8) :: bdim
        real(kind=8) :: dsidep(6, 6)
    end subroutine meobl1
end interface
