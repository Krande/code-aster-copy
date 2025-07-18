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

subroutine te0264(option, nomte)
    implicit none
#include "jeveux.h"
#include "asterfort/dfdm2d.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/fointe.h"
#include "asterfort/jevech.h"
!
    character(len=16) :: option, nomte
! ......................................................................
!    - FONCTION REALISEE:  CALCUL DES MATRICES ELEMENTAIRES
!                          OPTION : 'CHAR_THER_SOUR_F'
!                          ELEMENTS FOURIER
!
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
! ......................................................................
!
!-----------------------------------------------------------------------
    integer(kind=8) :: icode, nbres
    real(kind=8) :: soun, sounp1, theta
!-----------------------------------------------------------------------
    parameter(nbres=3)
    character(len=8) :: nompar(nbres)
    real(kind=8) :: valpar(nbres)
    real(kind=8) :: poids, r, z, sour
    integer(kind=8) :: nno, kp, npg1, i, k, itemps, ivectt, isour, nnos, jgano
    integer(kind=8) :: ipoids, ivf, idfde, igeom, ndim
!
!
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, &
                     npg=npg1, jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PINSTR', 'L', itemps)
    call jevech('PSOURCF', 'L', isour)
    call jevech('PVECTTR', 'E', ivectt)
    theta = zr(itemps+2)
    nompar(1) = 'X'
    nompar(2) = 'Y'
    nompar(3) = 'INST'
!
    do kp = 1, npg1
        k = (kp-1)*nno
        call dfdm2d(nno, kp, ipoids, idfde, zr(igeom), &
                    poids)
        r = 0.d0
        z = 0.d0
        do i = 1, nno
            r = r+zr(igeom+2*(i-1))*zr(ivf+k+i-1)
            z = z+zr(igeom+2*(i-1)+1)*zr(ivf+k+i-1)
        end do
        poids = poids*r
        valpar(1) = r
        valpar(2) = z
        valpar(3) = zr(itemps)
        call fointe('FM', zk8(isour), 3, nompar, valpar, &
                    sounp1, icode)
        valpar(3) = zr(itemps)-zr(itemps+1)
        call fointe('FM', zk8(isour), 3, nompar, valpar, &
                    soun, icode)
        if (theta < -0.5) then
            sour = sounp1
        else
            sour = theta*sounp1+(1.0d0-theta)*soun
        end if
        do i = 1, nno
            zr(ivectt+i-1) = zr(ivectt+i-1)+poids*zr(ivf+k+i-1)*sour
        end do
    end do
end subroutine
