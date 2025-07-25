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
subroutine te0269(option, nomte)
    implicit none
#include "jeveux.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/vff2dn.h"
!
    character(len=16) :: option, nomte
! ......................................................................
!    - FONCTION REALISEE:  CALCUL DES VECTEURS ELEMENTAIRES
!                          OPTION : 'CHAR_THER_ECHA_R'
!                          ELEMENTS FOURIER
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
! ......................................................................
!
    real(kind=8) :: poids, r, nx, ny, tpg, theta
    integer(kind=8) :: nno, kp, npg, ipoids, ivf, idfde, igeom
    integer(kind=8) :: ivectt, k, i, icoefh, itex
!
!-----------------------------------------------------------------------
    integer(kind=8) :: itemp, itemps, jgano, ndim, nnos
!-----------------------------------------------------------------------
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PCOEFHR', 'L', icoefh)
    call jevech('PT_EXTR', 'L', itex)
    call jevech('PTEMPER', 'L', itemp)
    call jevech('PINSTR', 'L', itemps)
    call jevech('PVECTTR', 'E', ivectt)
!
    theta = zr(itemps+2)
!
    do kp = 1, npg
        k = (kp-1)*nno
        call vff2dn(ndim, nno, kp, ipoids, idfde, &
                    zr(igeom), nx, ny, poids)
        r = 0.d0
        tpg = 0.d0
        do i = 1, nno
            r = r+zr(igeom+2*i-2)*zr(ivf+k+i-1)
            tpg = tpg+zr(itemp+i-1)*zr(ivf+k+i-1)
        end do
        poids = poids*r
        if (theta < -0.5) then
            do i = 1, nno
                zr(ivectt+i-1) = zr(ivectt+i-1)+poids*zr(ivf+k+i-1)*zr(icoefh)*(zr(itex)-tpg)
            end do
        else
            do i = 1, nno
                zr(ivectt+i-1) = zr(ivectt+i-1)+poids*zr(ivf+k+i-1)*zr(icoefh)*(zr(itex)-(1.0d0-&
                                &theta)*tpg)
            end do
        end if
    end do
end subroutine
