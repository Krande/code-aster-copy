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
    subroutine zerof2(func, funcp, para, nb_para, x0, xap, epsi, nitmax,&
                      solu, iret, n)
#include "asterf_types.h"
        real(kind=8) :: x0
        real(kind=8) :: para(nb_para)
        integer(kind=8), intent(in) :: nb_para
        real(kind=8) :: xap
        real(kind=8) :: epsi
        integer(kind=8) :: nitmax
        real(kind=8) :: solu
        integer(kind=8) :: iret
        integer(kind=8) :: n
        interface
        function funcp(x, param)
            real(kind=8), intent(in) :: x
            real(kind=8), intent(in) :: param(*)
            real(kind=8) :: funcp
        end function
        function func(x)
            real(kind=8) :: x
            real(kind=8) :: func
        end function
        end interface
    end subroutine zerof2
end interface
