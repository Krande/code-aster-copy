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
subroutine te0260(option, nomte)
    implicit none
#include "jeveux.h"
#include "asterfort/dfdm2d.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utmess.h"
    character(len=16) :: option, nomte
! ......................................................................
!    - FONCTION REALISEE:  CALCUL DES MATRICES ELEMENTAIRES
!                          OPTION : 'RIGI_THER'
!                          ELEMENTS FOURIER
!
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
! ......................................................................
!
    integer(kind=8) :: icodre(1)
    character(len=8) :: fami, poum
    real(kind=8) :: dfdr(9), dfdz(9), poids, r, valres(1)
    integer(kind=8) :: nno, kp, npg1, i, j, k, itemps, imattt, ndim, nnos, jgano
    integer(kind=8) :: ipoids, ivf, idfde, igeom, imate, kpg, spt
!
!
!-----------------------------------------------------------------------
    integer(kind=8) :: iharm, ij, nh
    real(kind=8) :: r2, wij, xh, xh2
!-----------------------------------------------------------------------
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg1, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PHARMON', 'L', iharm)
    nh = zi(iharm)
    if (nh .eq. -1) then
        call utmess('F', 'ELEMENTS3_63')
    end if
    xh = dble(nh)
    xh2 = xh*xh
    call jevech('PMATERC', 'L', imate)
    call jevech('PMATTTR', 'E', imattt)
    call jevech('PINSTR', 'L', itemps)
!
    fami = 'FPG1'
    kpg = 1
    spt = 1
    poum = '+'
    call rcvalb(fami, kpg, spt, poum, zi(imate), &
                ' ', 'THER', 1, 'INST', [zr(itemps)], &
                1, 'LAMBDA', valres, icodre, 1)
!
    do kp = 1, npg1
        k = (kp-1)*nno
        call dfdm2d(nno, kp, ipoids, idfde, zr(igeom), &
                    poids, dfdr, dfdz)
!
        r = 0.d0
        do i = 1, nno
            r = r+zr(igeom+2*(i-1))*zr(ivf+k+i-1)
        end do
        r2 = r*r
        poids = poids*r
!
        ij = imattt-1
        do i = 1, nno
!
            do j = 1, i
                wij = zr(ivf+k+i-1)*zr(ivf+k+j-1)
                ij = ij+1
                zr(ij) = zr(ij)+poids*valres(1)*(dfdr(i)*dfdr(j)+dfdz(i)*dfdz(j)&
                         & +xh2*wij/r2)
            end do
        end do
    end do
end subroutine
