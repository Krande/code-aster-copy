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

subroutine te0196(option, nomte)
    implicit none
#include "jeveux.h"
#include "asterfort/dfdm2d.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/rcvalb.h"
!
    character(len=16) :: option, nomte
! ......................................................................
!    - FONCTION REALISEE:  CALCUL DES TERMES ELEMENTAIRES EN MECANIQUE
!                          OPTION : 'CHAR_ME_PESANR  '
!                          ELEMENTS AXISYMETRIQUES FOURIER
!
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
! ......................................................................
!
    integer(kind=8) :: icodre(1)
    real(kind=8) :: poids, r
    character(len=8) :: fami, poum
    integer(kind=8) :: nno, kp, npg1, i, ivectu, ipesa, ndim, nnos, jgano
    integer(kind=8) :: ipoids, ivf, idfde, igeom, imate, kpg, spt
!
!
!-----------------------------------------------------------------------
    integer(kind=8) :: k
    real(kind=8) :: rho(1)
!-----------------------------------------------------------------------
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, &
                     npg=npg1, jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PMATERC', 'L', imate)
    call jevech('PPESANR', 'L', ipesa)
    call jevech('PVECTUR', 'E', ivectu)
!
    fami = 'FPG1'
    kpg = 1
    spt = 1
    poum = '+'
    call rcvalb(fami, kpg, spt, poum, zi(imate), &
                ' ', 'ELAS', 0, ' ', [0.d0], &
                1, 'RHO', rho, icodre, 1)
!
    do kp = 1, npg1
        k = (kp-1)*nno
        call dfdm2d(nno, kp, ipoids, idfde, zr(igeom), &
                    poids)
!
        r = 0.d0
        do i = 1, nno
            r = r+zr(igeom+2*(i-1))*zr(ivf+k+i-1)
        end do
        poids = poids*r
!
        poids = poids*rho(1)
        do i = 1, nno
            k = (kp-1)*nno
            zr(ivectu+3*i-3) = zr(ivectu+3*i-3)+poids*zr(ipesa)*zr(ipesa+1)*zr(ivf+k+i-1)
            zr(ivectu+3*i-2) = zr(ivectu+3*i-2)+poids*zr(ipesa)*zr(ipesa+2)*zr(ivf+k+i-1)
            zr(ivectu+3*i-1) = zr(ivectu+3*i-1)+poids*zr(ipesa)*zr(ipesa+3)*zr(ivf+k+i-1)
        end do
    end do
end subroutine
