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

subroutine ntcomp(icomp, icamas, ndim, temp, dtemp, coorpg, aniso, ifon, fluxglo)
!.
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "jeveux.h"
#include "asterfort/rcfode.h"
#include "asterfort/matrot.h"
#include "asterfort/utrcyl.h"
#include "asterfort/utpvlg.h"
#include "asterfort/utpvgl.h"
#include "asterc/r8dgrd.h"
!
    integer, intent(in) :: icomp, icamas, ndim, ifon(6)
    real(kind=8), intent(in) :: temp, dtemp(3), coorpg(3)
    aster_logical, intent(in) :: aniso
    real(kind=8), intent(out) :: fluxglo(3)
!
    aster_logical :: global
    integer :: n1, n2
    real(kind=8) :: lambor(3), orig(3), lambda, r8bid, xu, yu, xnorm
    real(kind=8) :: alpha, angl(3), p(3, 3), dire(3), fluloc(3)
    real(kind=8) :: aalpha, abeta
!
    fluxglo = 0.d0
    global = ASTER_FALSE
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
            fluxglo(1:ndim) = lambda*dtemp(1:ndim)
        end if
!
! ------- TRAITEMENT DE L ANISOTROPIE
!
        if (aniso) then
            if (zr(icamas) .gt. 0.d0) then
                global = ASTER_TRUE
                if (ndim == 3) then
                    angl(1) = zr(icamas+1)*r8dgrd()
                    angl(2) = zr(icamas+2)*r8dgrd()
                    angl(3) = zr(icamas+3)*r8dgrd()
                    call matrot(angl, p)
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
                end if
            end if
            if (ndim == 3) then
                if (.not. global) then
                    call utrcyl(coorpg, dire, orig, p)
                end if
                fluxglo(1) = dtemp(1)
                fluxglo(2) = dtemp(2)
                fluxglo(3) = dtemp(3)
                n1 = 1
                n2 = 3
                call utpvgl(n1, n2, p, fluxglo, fluloc)
                fluloc(1) = lambor(1)*fluloc(1)
                fluloc(2) = lambor(2)*fluloc(2)
                fluloc(3) = lambor(3)*fluloc(3)
                n1 = 1
                n2 = 3
                call utpvlg(n1, n2, p, fluloc, fluxglo)
            else
                if (.not. global) then
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
                fluxglo(1) = dtemp(1)
                fluxglo(2) = dtemp(2)
                fluloc(1) = p(1, 1)*dtemp(1)+p(2, 1)*dtemp(2)
                fluloc(2) = p(1, 2)*dtemp(1)+p(2, 2)*dtemp(2)
                fluloc(1) = lambor(1)*fluloc(1)
                fluloc(2) = lambor(2)*fluloc(2)
                fluxglo(1) = p(1, 1)*fluloc(1)+p(1, 2)*fluloc(2)
                fluxglo(2) = p(2, 1)*fluloc(1)+p(2, 2)*fluloc(2)
            end if
        end if
!
    else
        ASSERT(ASTER_FALSE)
    end if
!
end subroutine
