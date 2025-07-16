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

subroutine matrRotLGTher(aniso, ndim, coorpg, matr)
!.
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/matrot.h"
#include "asterfort/rcangm.h"
!
    integer(kind=8), intent(in) :: ndim
    real(kind=8), intent(in) :: coorpg(3)
    aster_logical, intent(in) :: aniso
    real(kind=8), intent(out) :: matr(3, 3)
!
    integer(kind=8) :: j
    real(kind=8) :: alpha, angl(3), p(3, 3)
!
    matr = 0.d0
    if (aniso) then
        p = 0.d0
        call rcangm(ndim, coorpg, angl)
        if (ndim == 3) then
            call matrot(angl, p)
            p = transpose(p)
        else
            alpha = angl(1)
            p(1, 1) = cos(alpha)
            p(2, 1) = sin(alpha)
            p(1, 2) = -sin(alpha)
            p(2, 2) = cos(alpha)
        end if
        matr = p
    else
        do j = 1, ndim
            matr(j, j) = 1.d0
        end do
    end if
!
end subroutine
