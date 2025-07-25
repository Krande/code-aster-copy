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

subroutine te0017(option, nomte)
    implicit none
#include "jeveux.h"
#include "asterfort/dfdm3d.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/fointe.h"
#include "asterfort/jevech.h"
!
    character(len=16) :: option, nomte
!.......................................................................
!
!     BUT: CALCUL DES VECTEURS ELEMENTAIRES EN MECANIQUE
!          ELEMENTS ISOPARAMETRIQUES 3D
!
!          OPTION : 'CHAR_MECA_FORC_F '
!
!     ENTREES  ---> OPTION : OPTION DE CALCUL
!              ---> NOMTE  : NOM DU TYPE ELEMENT
!.......................................................................
!
!
    integer(kind=8) :: ipoids, ivf, idfde, igeom, itemps, iforc, ier
    integer(kind=8) :: jgano, nno, ndl, kp, npg1, ii, i, l, ivectu, ndim, nnos
    real(kind=8) :: poids, fx, fy, fz
    real(kind=8) :: xx, yy, zz, valpar(4)
    character(len=8) :: nompar(4)
!     ------------------------------------------------------------------
!
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, &
                     npg=npg1, jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PVECTUR', 'E', ivectu)
    call jevech('PINSTR', 'L', itemps)
    call jevech('PFF3D3D', 'L', iforc)
!
    valpar(4) = zr(itemps)
    nompar(1) = 'X'
    nompar(2) = 'Y'
    nompar(3) = 'Z'
    nompar(4) = 'INST'
!
    ndl = 3*nno
    do i = 1, ndl
        zr(ivectu+i-1) = 0.0d0
    end do
!
!    BOUCLE SUR LES POINTS DE GAUSS
!
    do kp = 1, npg1
!
        l = (kp-1)*nno
        call dfdm3d(nno, kp, ipoids, idfde, zr(igeom), &
                    poids)
!
        xx = 0.d0
        yy = 0.d0
        zz = 0.d0
        do i = 1, nno
            xx = xx+zr(igeom+3*i-3)*zr(ivf+l+i-1)
            yy = yy+zr(igeom+3*i-2)*zr(ivf+l+i-1)
            zz = zz+zr(igeom+3*i-1)*zr(ivf+l+i-1)
        end do
        valpar(1) = xx
        valpar(2) = yy
        valpar(3) = zz
        call fointe('FM', zk8(iforc), 4, nompar, valpar, &
                    fx, ier)
        call fointe('FM', zk8(iforc+1), 4, nompar, valpar, &
                    fy, ier)
        call fointe('FM', zk8(iforc+2), 4, nompar, valpar, &
                    fz, ier)
!
        do i = 1, nno
            ii = 3*(i-1)
            zr(ivectu+ii) = zr(ivectu+ii)+poids*zr(ivf+l+i-1)*fx
            zr(ivectu+ii+1) = zr(ivectu+ii+1)+poids*zr(ivf+l+i-1)*fy
            zr(ivectu+ii+2) = zr(ivectu+ii+2)+poids*zr(ivf+l+i-1)*fz
!
        end do
!
    end do
!
end subroutine
