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
    subroutine pk2sig(ndim, f, jac, pk2, sig,&
                      ind)
        integer(kind=8), intent(in) :: ndim, ind
        real(kind=8), intent(in) :: f(3,3), jac
        real(kind=8), intent(inout) :: pk2(2*ndim), sig(2*ndim)
    end subroutine pk2sig
end interface
