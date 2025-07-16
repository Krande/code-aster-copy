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
subroutine te0127(option, nomte)
    implicit none
#include "jeveux.h"
#include "asterc/r8t0.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
!
    character(len=16) :: option, nomte
! ......................................................................
!    - FONCTION REALISEE:  CALCUL DES VECTEURS RESIDUS
!                          OPTION : 'RESI_THER_COEF_R'
!                          OPTION : 'RESI_THER_RAYO_R'
!                          ELEMENTS 3D
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
!
    real(kind=8) :: nx, ny, nz, sx(9, 9), sy(9, 9), sz(9, 9), jac, tpg
    integer(kind=8) :: ipoids, ivf, idfdx, idfdy, igeom, jgano
    integer(kind=8) :: ndim, nno, ipg, npg1, iveres, iech, iray, nnos
!
    integer(kind=8) :: idec, jdec, kdec, ldec
    real(kind=8) :: hech, sigma, epsil, tz0
! DEB ------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ino, itemp, itemps, j, jno
!-----------------------------------------------------------------------
    tz0 = r8t0()
!
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg1, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfdx, jgano=jgano)
    idfdy = idfdx+1
!
    if (option(11:14) .eq. 'COEF') then
        call jevech('PCOEFHR', 'L', iech)
        hech = zr(iech)
    else if (option(11:14) .eq. 'RAYO') then
        call jevech('PRAYONR', 'L', iray)
        sigma = zr(iray)
        epsil = zr(iray+1)
    end if
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PINSTR', 'L', itemps)
    call jevech('PTEMPEI', 'L', itemp)
    call jevech('PRESIDU', 'E', iveres)
!
!    CALCUL DES PRODUITS VECTORIELS OMI   OMJ
!
    do ino = 1, nno
        i = igeom+3*(ino-1)-1
        do jno = 1, nno
            j = igeom+3*(jno-1)-1
            sx(ino, jno) = zr(i+2)*zr(j+3)-zr(i+3)*zr(j+2)
            sy(ino, jno) = zr(i+3)*zr(j+1)-zr(i+1)*zr(j+3)
            sz(ino, jno) = zr(i+1)*zr(j+2)-zr(i+2)*zr(j+1)
        end do
    end do
!
    do ipg = 1, npg1
        kdec = (ipg-1)*nno*ndim
        ldec = (ipg-1)*nno
!
        nx = 0.0d0
        ny = 0.0d0
        nz = 0.0d0
        do i = 1, nno
            idec = (i-1)*ndim
            do j = 1, nno
                jdec = (j-1)*ndim
                nx = nx+zr(idfdx+kdec+idec)*zr(idfdy+kdec+jdec)*sx(i, j)
                ny = ny+zr(idfdx+kdec+idec)*zr(idfdy+kdec+jdec)*sy(i, j)
                nz = nz+zr(idfdx+kdec+idec)*zr(idfdy+kdec+jdec)*sz(i, j)
            end do
        end do
        jac = sqrt(nx*nx+ny*ny+nz*nz)
!
        tpg = 0.d0
        do i = 1, nno
            tpg = tpg+zr(itemp+i-1)*zr(ivf+ldec+i-1)
        end do
        if (option(11:14) .eq. 'COEF') then
            do i = 1, nno
                zr(iveres+i-1) = zr(iveres+i-1)+jac*zr(ipoids+ipg-1)*zr(ivf+ldec+i-1)&
                                &*hech*tpg
            end do
        else if (option(11:14) .eq. 'RAYO') then
            do i = 1, nno
                zr(iveres+i-1) = zr(iveres+i-1)+jac*zr(ipoids+ipg-1)*zr(ivf+ldec+i-1)&
                                &*sigma*epsil*(tpg+tz0)**4
            end do
        end if
!
    end do
! FIN ------------------------------------------------------------------
end subroutine
