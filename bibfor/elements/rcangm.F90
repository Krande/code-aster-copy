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
subroutine rcangm(ndim, coor, angl_naut)
    implicit none
#include "jeveux.h"
#include "asterc/r8dgrd.h"
#include "asterfort/angvx.h"
#include "asterfort/angvxy.h"
#include "asterfort/tecach.h"
#include "asterfort/utrcyl.h"
    integer(kind=8) :: ndim
    real(kind=8) :: angl_naut(3), coor(3)
! ......................................................................
!    - ORIENTATION DU MASSIF
!
!   IN      NDIM    I      : DIMENSION DU PROBLEME
!   IN      COOR    R        COORDONNEE DU POINT
!                            (CAS CYLINDRIQUE)
!   OUT     ANGL_NAUT R    : ANGLE NAUTIQUE
! ......................................................................
    integer(kind=8) :: icamas, iret, i
    real(kind=8) :: p(3, 3), xg(3), yg(3), orig(3), dire(3)
    real(kind=8) :: alpha, beta, xu, yu, xnorm
!     ------------------------------------------------------------------
!
    call tecach('NNO', 'PCAMASS', 'L', iret, iad=icamas)
    angl_naut(:) = 0.d0
!
    if (iret .eq. 0) then
        if (zr(icamas) .gt. 0.d0) then
            angl_naut(1) = zr(icamas+1)*r8dgrd()
            if (ndim .eq. 3) then
                angl_naut(2) = zr(icamas+2)*r8dgrd()
                angl_naut(3) = zr(icamas+3)*r8dgrd()
            end if
!
        else if (abs(zr(icamas)+1.d0) .lt. 1.d-3) then
!
! ON TRANSFORME LA DONNEE DU REPERE CYLINDRIQUE EN ANGLE NAUTIQUE
!
            orig(1:ndim) = zr(icamas+3+1:icamas+3+ndim)
            if (ndim .eq. 3) then
                alpha = zr(icamas+1)*r8dgrd()
                beta = zr(icamas+2)*r8dgrd()
                dire(1) = cos(alpha)*cos(beta)
                dire(2) = sin(alpha)*cos(beta)
                dire(3) = -sin(beta)
                call utrcyl(coor, dire, orig, p)
                do i = 1, 3
                    xg(i) = p(1, i)
                    yg(i) = p(2, i)
                end do
                call angvxy(xg, yg, angl_naut)
            else
                xu = coor(1)-orig(1)
                yu = coor(2)-orig(2)
                xnorm = sqrt(xu**2+yu**2)
                xu = xu/xnorm
                yu = yu/xnorm
                p(1, 1) = xu
                p(2, 1) = yu
                p(1, 2) = -yu
                p(2, 2) = xu
                xg(1) = xu
                xg(2) = yu
                xg(3) = 0.d0
                call angvx(xg, alpha, beta)
                angl_naut(1) = alpha
            end if
        end if
    end if
!
end subroutine
