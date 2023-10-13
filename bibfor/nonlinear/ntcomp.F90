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

subroutine ntcomp(icomp, icamas, ndim, temp, dtemp, coorpg, aniso, ifon, fluxglo, Kglo)
!.
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "jeveux.h"
#include "asterfort/rcfode.h"
#include "asterfort/matrRotLGTher.h"
!
    integer, intent(in) :: icomp, icamas, ndim, ifon(6)
    real(kind=8), intent(in) :: temp, dtemp(3), coorpg(3)
    aster_logical, intent(in) :: aniso
    real(kind=8), intent(out) :: fluxglo(3)
    real(kind=8), intent(out) :: Kglo(3, 3)
!
    integer :: j
    real(kind=8) :: lambor(3), lambda, r8bid
    real(kind=8) ::  p(3, 3), Kloc(3, 3)
!
    fluxglo = 0.d0
    Kglo = 0.d0
!
    if (zk16(icomp) (1:5) .eq. 'THER_') then
!
! ------- EVALUATION DE LA CONDUCTIVITE LAMBDA
!
        lambor = 0.d0
        if (aniso) then
            call rcfode(ifon(4), temp, lambor(1), r8bid)
            call rcfode(ifon(5), temp, lambor(2), r8bid)
            if (ndim == 3) then
                call rcfode(ifon(6), temp, lambor(3), r8bid)
            end if
        else
            call rcfode(ifon(2), temp, lambda, r8bid)
        end if
!
! ------- TRAITEMENT DE L ANISOTROPIE
!
        if (aniso) then
            call matrRotLGTher(aniso, icamas, ndim, coorpg, p)
            Kloc = transpose(p)
            do j = 1, ndim
                Kloc(j, 1:3) = lambor(j)*Kloc(j, 1:3)
            end do
            Kglo = matmul(p, Kloc)
        else
            do j = 1, ndim
                Kglo(j, j) = lambda
            end do
        end if
        fluxglo = matmul(Kglo, dtemp)
!
    else
        ASSERT(ASTER_FALSE)
    end if
!
end subroutine
