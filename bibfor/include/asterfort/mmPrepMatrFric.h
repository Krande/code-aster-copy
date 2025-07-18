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
    subroutine mmPrepMatrFric(ndim  , nbcps ,&
                              tau1  , tau2  , mprojt,&
                              rese  , nrese ,&
                              dlagrf, djeut ,&
                              e     , a     ,&
                              b     , d     ,&
                              r     , tt    ,&
                              dlagft, pdlaft,&
                              pdjeut, prese)
        integer(kind=8), intent(in) :: ndim, nbcps
        real(kind=8), intent(in) :: tau1(3), tau2(3), mprojt(3, 3)
        real(kind=8), intent(in) :: rese(3), nrese
        real(kind=8), intent(in) :: dlagrf(2), djeut(3)
        real(kind=8), intent(out) :: e(3, 3), a(2, 3), b(2, 3)
        real(kind=8), intent(out) :: d(3, 3), r(2, 2), tt(3, 3)
        real(kind=8), intent(out) :: dlagft(3), pdlaft(3), pdjeut(3), prese(3)
    end subroutine mmPrepMatrFric
end interface
