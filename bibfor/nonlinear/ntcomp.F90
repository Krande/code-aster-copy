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

subroutine ntcomp(rela_name, ndim, temp, dtemp, coorpg, aniso, ifon, fluxglo, Kglo, dfluxglo)
!.
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "jeveux.h"
#include "asterfort/rcfode.h"
#include "asterfort/matrRotLGTher.h"
!
    character(len=16), intent(in) :: rela_name
    integer(kind=8), intent(in) :: ndim, ifon(6)
    real(kind=8), intent(in) :: temp, dtemp(3), coorpg(3)
    aster_logical, intent(in) :: aniso
    real(kind=8), intent(out) :: fluxglo(3)
    real(kind=8), intent(out) :: Kglo(3, 3), dfluxglo(3)
!
    integer(kind=8) :: j
    real(kind=8) :: lambor(3), dlambor(3), lambda, dlambda
    real(kind=8) ::  p(3, 3), Kloc(3, 3), Kloc2(3, 3)
!
    fluxglo = 0.d0
    Kglo = 0.d0
    dfluxglo = 0.d0
!
    if (rela_name(1:5) .eq. 'THER_') then
!
! ------- EVALUATION DE LA CONDUCTIVITE LAMBDA
!
        lambor = 0.d0
        dlambor = 0.d0
        if (aniso) then
            call rcfode(ifon(4), temp, lambor(1), dlambor(1))
            call rcfode(ifon(5), temp, lambor(2), dlambor(2))
            if (ndim == 3) then
                call rcfode(ifon(6), temp, lambor(3), dlambor(3))
            end if
        else
            call rcfode(ifon(2), temp, lambda, dlambda)
        end if
!
! ------- TRAITEMENT DE L ANISOTROPIE
!
        if (aniso) then
            call matrRotLGTher(aniso, ndim, coorpg, p)
            Kloc = transpose(p)
            Kloc2 = transpose(p)
            do j = 1, ndim
                Kloc(j, 1:3) = lambor(j)*Kloc(j, 1:3)
                Kloc2(j, 1:3) = dlambor(j)*Kloc2(j, 1:3)
            end do
            Kglo = matmul(p, Kloc)
            fluxglo = matmul(Kglo, dtemp)
!
            Kloc2 = matmul(p, Kloc2)
            dfluxglo = matmul(Kloc2, dtemp)
        else
            if (ndim == 3) then
                Kglo(1, 1) = lambda
                Kglo(2, 2) = lambda
                Kglo(3, 3) = lambda
!
                dfluxglo = dlambda*dtemp
!
                fluxglo(1) = lambda*dtemp(1)
                fluxglo(2) = lambda*dtemp(2)
                fluxglo(3) = lambda*dtemp(3)
            else
                Kglo(1, 1) = lambda
                Kglo(2, 2) = lambda
!
                dfluxglo = dlambda*dtemp
!
                fluxglo(1) = lambda*dtemp(1)
                fluxglo(2) = lambda*dtemp(2)
            end if
        end if
!
    else
        ASSERT(ASTER_FALSE)
    end if
!
end subroutine
