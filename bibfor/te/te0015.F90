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

subroutine te0015(option, nomte)
!
    implicit none
#include "jeveux.h"
#include "asterfort/dfdm3d.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcvalb.h"
!
    character(len=16) :: nomte, option
!.......................................................................
!
!     BUT: CALCUL DES VECTEURS ELEMENTAIRES EN MECANIQUE
!          ELEMENTS ISOPARAMETRIQUES 3D
!
!          OPTION : 'CHAR_MECA_PESA_R '
!
!     ENTREES  ---> OPTION : OPTION DE CALCUL
!          ---> NOMTE  : NOM DU TYPE ELEMENT
!.......................................................................
!
!
    integer(kind=8) :: icodre(1)
    character(len=16) :: phenom
    real(kind=8) :: rho(1), coef
    real(kind=8) :: poids
    integer(kind=8) :: ipoids, ivf, idfde, igeom
    integer(kind=8) :: jgano, imate, ipesa, ivectu, nnos
    integer(kind=8) :: ndim, nno, npg, ndl, kp, l, i, ii, j
!
!
!
!
!
!
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, &
                     npg=npg, jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PMATERC', 'L', imate)
    call jevech('PPESANR', 'L', ipesa)
    call jevech('PVECTUR', 'E', ivectu)
!
    call rccoma(zi(imate), 'ELAS', 1, phenom, icodre(1))
    call rcvalb('FPG1', 1, 1, '+', zi(imate), &
                ' ', phenom, 0, ' ', [0.d0], &
                1, 'RHO', rho, icodre(1), 1)
!
    ndl = 3*nno
    do i = 1, ndl
        zr(ivectu+i-1) = 0.0d0
    end do
!
!    BOUCLE SUR LES POINTS DE GAUSS
!
    do kp = 1, npg
!
        l = (kp-1)*nno
        call dfdm3d(nno, kp, ipoids, idfde, zr(igeom), &
                    poids)
!
        coef = rho(1)*poids*zr(ipesa)
!
        do i = 1, nno
            ii = 3*(i-1)
!
            do j = 1, 3
                zr(ivectu+ii+j-1) = zr(ivectu+ii+j-1)+coef*zr(ivf+l+i-1)*zr(ipesa+j)
            end do
!
        end do
!
    end do
!
!
end subroutine
