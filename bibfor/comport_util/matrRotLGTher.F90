! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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

subroutine matrRotLGTher(aniso, icamas, ndim, coorpg, matr)
!.
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "jeveux.h"
#include "asterfort/matrot.h"
#include "asterfort/utrcyl.h"
#include "asterc/r8dgrd.h"
!
    integer, intent(in) :: icamas, ndim
    real(kind=8), intent(in) :: coorpg(3)
    aster_logical, intent(in) :: aniso
    real(kind=8), intent(out) :: matr(3, 3)
!
    integer ::  j
    real(kind=8) :: orig(3), xu, yu, xnorm
    real(kind=8) :: alpha, angl(3), p(3, 3), dire(3)
    real(kind=8) :: aalpha, abeta
!
    matr = 0.d0
    if (aniso) then
        p = 0.d0
        if (zr(icamas) .gt. 0.d0) then
            if (ndim == 3) then
                angl(1) = zr(icamas+1)*r8dgrd()
                angl(2) = zr(icamas+2)*r8dgrd()
                angl(3) = zr(icamas+3)*r8dgrd()
                call matrot(angl, p)
                p = transpose(p)
            else
                alpha = zr(icamas+1)*r8dgrd()
                p(1, 1) = cos(alpha)
                p(2, 1) = sin(alpha)
                p(1, 2) = -sin(alpha)
                p(2, 2) = cos(alpha)
            end if
        else
            orig(1:ndim) = zr(icamas+3+1:icamas+3+ndim)
            if (ndim == 3) then
                aalpha = zr(icamas+1)*r8dgrd()
                abeta = zr(icamas+2)*r8dgrd()
                dire(1) = cos(aalpha)*cos(abeta)
                dire(2) = sin(aalpha)*cos(abeta)
                dire(3) = -sin(abeta)
                call utrcyl(coorpg, dire, orig, p)
            else
                xu = orig(1)-coorpg(1)
                yu = orig(2)-coorpg(2)
                xnorm = sqrt(xu**2+yu**2)
                xu = xu/xnorm
                yu = yu/xnorm
                p(1, 1) = xu
                p(2, 1) = yu
                p(1, 2) = -yu
                p(2, 2) = xu
            end if
        end if
        matr = p
    else
        do j = 1, ndim
            matr(j, j) = 1.d0
        end do
    end if
!
end subroutine
