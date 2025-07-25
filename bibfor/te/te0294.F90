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
subroutine te0294(option, nomte)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dfdm2d.h"
#include "asterfort/dfdm3d.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
!
    character(len=16) :: option, nomte
! person_in_charge: josselin.delmas at edf.fr
!
!     BUT:
!         CALCUL DES VECTEURS ELEMENTAIRES
!         OPTION : 'SECM_ZZ1'
!
! ......................................................................
!
!
!
!
    integer(kind=8) :: i, k, kp, nno, nnos, npg, ndim, nbcmp
    integer(kind=8) :: ipoids, ivf, idfde, jgano, igeom, isief
    integer(kind=8) :: ivect1, ivect2, ivect3, ivect4, ivect5, ivect6
!
    real(kind=8) :: dfdx(27), dfdy(27), dfdz(27), poids, r
!
    aster_logical :: laxi
!
! ----------------------------------------------------------------------
!
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
!
    if (ndim .eq. 2) then
        nbcmp = 4
    else if (ndim .eq. 3) then
        nbcmp = 6
    else
        ASSERT(.false.)
    end if
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PSIEF_R', 'L', isief)
    call jevech('PVECTR1', 'E', ivect1)
    call jevech('PVECTR2', 'E', ivect2)
    call jevech('PVECTR3', 'E', ivect3)
    call jevech('PVECTR4', 'E', ivect4)
    if (ndim .eq. 3) then
        call jevech('PVECTR5', 'E', ivect5)
        call jevech('PVECTR6', 'E', ivect6)
    end if
!
    laxi = .false.
    if (lteatt('AXIS', 'OUI')) laxi = .true.
!
    do kp = 1, npg
        k = (kp-1)*nno
        if (ndim .eq. 2) then
            call dfdm2d(nno, kp, ipoids, idfde, zr(igeom), &
                        poids, dfdx, dfdy)
        else
            call dfdm3d(nno, kp, ipoids, idfde, zr(igeom), &
                        poids, dfdx, dfdy, dfdz)
        end if
!
        if (laxi) then
            r = 0.d0
            do i = 1, nno
                r = r+zr(igeom+2*(i-1))*zr(ivf+k+i-1)
            end do
            poids = poids*r
        end if
!
        do i = 1, nno
            k = (kp-1)*nno
            zr(ivect1+i-1) = zr(ivect1+i-1)+poids*zr(ivf+k+i-1)*zr(isief+nbcmp*(kp-1))
            zr(ivect2+i-1) = zr(ivect2+i-1)+poids*zr(ivf+k+i-1)*zr(isief+nbcmp*(kp-1)+1)
            zr(ivect3+i-1) = zr(ivect3+i-1)+poids*zr(ivf+k+i-1)*zr(isief+nbcmp*(kp-1)+2)
            zr(ivect4+i-1) = zr(ivect4+i-1)+poids*zr(ivf+k+i-1)*zr(isief+nbcmp*(kp-1)+3)
            if (ndim .eq. 3) then
                zr(ivect5+i-1) = zr(ivect5+i-1)+poids*zr(ivf+k+i-1)*zr(isief+nbcmp*(kp-1)+&
                                 &4)
                zr(ivect6+i-1) = zr(ivect6+i-1)+poids*zr(ivf+k+i-1)*zr(isief+nbcmp*(kp-1)+&
                                 &5)
            end if
        end do
    end do
!
end subroutine
