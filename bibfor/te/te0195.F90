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
subroutine te0195(option, nomte)
    implicit none
#include "jeveux.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/fointe.h"
#include "asterfort/jevech.h"
#include "asterfort/vff2dn.h"
!
    character(len=16) :: option, nomte
! ......................................................................
!    - FONCTION REALISEE:  CALCUL DES VECTEURS ELEMENTAIRES
!                          OPTION : 'CHAR_ME_PF1D2D  '
!                          ELEMENTS AXI FOURIER
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
! ......................................................................
!
!-----------------------------------------------------------------------
    integer(kind=8) :: icode, jgano, nbres, nddl, ndim, nnos
!-----------------------------------------------------------------------
    parameter(nbres=3)
    character(len=8) :: nompar(nbres)
    real(kind=8) :: valpar(nbres), poids, r, tx, ty, tz, z, nx, ny
    integer(kind=8) :: nno, kp, npg, ipoids, ivf, idfde, igeom
    integer(kind=8) :: itemps, ivectu, k, i, l, iforc
!
!
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PINSTR', 'L', itemps)
    call jevech('PVECTUR', 'E', ivectu)
!
    nompar(1) = 'X'
    nompar(2) = 'Y'
    nompar(3) = 'INST'
    valpar(3) = zr(itemps)
    call jevech('PFF1D2D', 'L', iforc)
    nddl = 3
!
    do kp = 1, npg
        k = (kp-1)*nno
        call vff2dn(ndim, nno, kp, ipoids, idfde, &
                    zr(igeom), nx, ny, poids)
        r = 0.d0
        z = 0.d0
        do i = 1, nno
            l = (kp-1)*nno+i
            r = r+zr(igeom+2*i-2)*zr(ivf+l-1)
            z = z+zr(igeom+2*i-1)*zr(ivf+l-1)
        end do
        poids = poids*r
        valpar(1) = r
        valpar(2) = z
        call fointe('FM', zk8(iforc), 3, nompar, valpar, &
                    tx, icode)
        call fointe('FM', zk8(iforc+1), 3, nompar, valpar, &
                    ty, icode)
        call fointe('FM', zk8(iforc+2), 3, nompar, valpar, &
                    tz, icode)
        do i = 1, nno
            zr(ivectu+nddl*(i-1)) = zr(ivectu+nddl*(i-1))+tx*zr(ivf+k+i-1)*poids
            zr(ivectu+nddl*(i-1)+1) = zr(ivectu+nddl*(i-1)+1)+ty*zr(ivf+k+i-1)*poids
            zr(ivectu+nddl*(i-1)+2) = zr(ivectu+nddl*(i-1)+2)+tz*zr(ivf+k+i-1)*poids
        end do
    end do
end subroutine
