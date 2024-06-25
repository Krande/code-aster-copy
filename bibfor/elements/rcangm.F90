! --------------------------------------------------------------------
! Copyright (C) 1991 - 2024 - EDF R&D - www.code-aster.org
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
#include "asterc/r8nnem.h"
#include "asterfort/angvxy.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "asterfort/utrcyl.h"
    integer :: ndim
    real(kind=8) :: angl_naut(3), coor(3)
! ......................................................................
!    - ORIENTATION DU MASSIF
!
!   IN      NDIM    I      : DIMENSION DU PROBLEME
!   IN      COOR    R        COORDONNEE DU POINT
!                            (CAS CYLINDRIQUE)
!   OUT     ANGL_NAUT R    : ANGLE NAUTIQUE
! ......................................................................
    integer :: icamas, iret, i
    real(kind=8) :: p(3, 3), xg(3), yg(3), orig(3), dire(3)
    real(kind=8) :: alpha, beta
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
! (EN 3D, EN 2D ON MET A 0)
!
            if (ndim .eq. 3) then
                alpha = zr(icamas+1)*r8dgrd()
                beta = zr(icamas+2)*r8dgrd()
                dire(1) = cos(alpha)*cos(beta)
                dire(2) = sin(alpha)*cos(beta)
                dire(3) = -sin(beta)
                orig(1) = zr(icamas+4)
                orig(2) = zr(icamas+5)
                orig(3) = zr(icamas+6)
                call utrcyl(coor, dire, orig, p)
                do i = 1, 3
                    xg(i) = p(1, i)
                    yg(i) = p(2, i)
                end do
                call angvxy(xg, yg, angl_naut)
            else
                call utmess('F', 'ELEMENTS2_38')
            end if
        end if
    end if
!
end subroutine
