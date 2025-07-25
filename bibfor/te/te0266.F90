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
subroutine te0266(option, nomte)
    implicit none
#include "jeveux.h"
#include "asterfort/dfdm2d.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/rcvalb.h"
    character(len=16) :: option, nomte
!
!     BUT:
!       CALCUL DES FLUX DE TEMPERATURE AUX POINTS DE GAUSS
!       ELEMENTS 2D AXI
!       OPTION : 'FLUX_ELGA'
!
! ---------------------------------------------------------------------
!
!
!
    integer(kind=8) :: icodre(1)
    integer(kind=8) :: nno, kp, i, k, itempe, itemp, iflux, iharm, nh
    integer(kind=8) :: ipoids, ivf, idfde, igeom, imate
    integer(kind=8) :: npg, nnos, jgano, ndim, kpg, spt, j, nbcmp
!
    real(kind=8) :: valres(1), fluxr, fluxz, fluxt
    real(kind=8) :: dfdr(9), dfdz(9), poids, xh, r
!
    character(len=8) :: fami, poum
!
!-----------------------------------------------------------------------
!
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PHARMON', 'L', iharm)
    nh = zi(iharm)
    xh = dble(nh)
    call jevech('PMATERC', 'L', imate)
    call jevech('PINSTR', 'L', itemp)
    call jevech('PTEMPER', 'L', itempe)
    call jevech('PFLUXPG', 'E', iflux)
!
    fami = 'FPG1'
    kpg = 1
    spt = 1
    poum = '+'
    nbcmp = 3
    call rcvalb(fami, kpg, spt, poum, zi(imate), &
                ' ', 'THER', 1, 'INST', [zr(itemp)], &
                1, 'LAMBDA', valres, icodre, 1)
!
    do kp = 1, npg
        k = (kp-1)*nno
        call dfdm2d(nno, kp, ipoids, idfde, zr(igeom), &
                    poids, dfdr, dfdz)
!
        r = 0.d0
        do i = 1, nno
            r = r+zr(igeom+2*i-2)*zr(ivf+k+i-1)
        end do
!
        fluxr = 0.0d0
        fluxz = 0.0d0
        fluxt = 0.0d0
        do j = 1, nno
            fluxr = fluxr+zr(itempe+j-1)*dfdr(j)
            fluxz = fluxz+zr(itempe+j-1)*dfdz(j)
            fluxt = fluxt-zr(itempe+j-1)*zr(ivf+k+j-1)*xh/r
        end do
!
        zr(iflux+(kp-1)*nbcmp-1+1) = -valres(1)*fluxr
        zr(iflux+(kp-1)*nbcmp-1+2) = -valres(1)*fluxz
        zr(iflux+(kp-1)*nbcmp-1+3) = -valres(1)*fluxt
!
    end do
!
end subroutine
