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
    subroutine mmdepm(nbdm  , ndim  ,&
                      nne   , nnm   ,&
                      jdepm , jdepde,&
                      ffe   , ffm   ,&
                      ddeple, ddeplm,&
                      deplme, deplmm)
        integer(kind=8), intent(in) :: nbdm, ndim, nne, nnm
        integer(kind=8), intent(in) :: jdepde, jdepm
        real(kind=8), intent(in) :: ffe(9), ffm(9)
        real(kind=8), intent(out) :: ddeple(3), deplme(3)
        real(kind=8), intent(out) :: ddeplm(3), deplmm(3)
    end subroutine mmdepm
end interface
