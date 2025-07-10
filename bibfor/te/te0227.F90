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
subroutine te0227(option, nomte)
!
    implicit none
!
#include "jeveux.h"
#include "asterc/r8depi.h"
#include "asterfort/dfdm1d.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/getDensity.h"
!
    character(len=16) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
!    - FONCTION REALISEE:  CALCUL DES MATRICES ELEMENTAIRES
!                          COQUE 1D
!                          OPTION : 'MASS_INER       '
!                          ELEMENT: MECXSE3
!
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: dfdx(3), r, rm, poids, cour, nx, ny, yg
    real(kind=8) :: rho, x(3), y(3), xxi, xyi, yyi
    real(kind=8) :: matine(6), volume, depi
    integer(kind=8) :: nno, ipoids, ivf, idfdk, igeom, imate, icaco
    integer(kind=8) :: kp, npg, i, j, k, lcastr
!
! --------------------------------------------------------------------------------------------------
!
    depi = r8depi()
!
    call elrefe_info(fami='RIGI', nno=nno, &
                     npg=npg, jpoids=ipoids, jvf=ivf, jdfde=idfdk)
!
!
    call jevech('PGEOMER', 'L', igeom)
    do i = 1, nno
        x(i) = zr(igeom-2+2*i)
        y(i) = zr(igeom-1+2*i)
    end do
!
    call jevech('PMATERC', 'L', imate)
    call jevech('PCACOQU', 'L', icaco)

    call getDensity(zi(imate), rho)
    rm = rho*zr(icaco)
!
    call jevech('PMASSINE', 'E', lcastr)
!
    volume = 0.d0
    matine = 0.d0
!
!     --- BOUCLE SUR LES POINTS DE GAUSS ---
!
    do kp = 1, npg
        k = (kp-1)*nno
        call dfdm1d(nno, zr(ipoids+kp-1), zr(idfdk+k), zr(igeom), dfdx, &
                    cour, poids, nx, ny)
        r = 0.d0
        do i = 1, nno
            r = r+zr(igeom+2*(i-1))*zr(ivf+k+i-1)
        end do
        poids = poids*r
        volume = volume+poids
!
        do i = 1, nno
!           --- CDG ---
            zr(lcastr+1) = zr(lcastr+1)+poids*x(i)*zr(ivf+k+i-1)
            zr(lcastr+2) = zr(lcastr+2)+poids*y(i)*zr(ivf+k+i-1)
!           --- INERTIE ---
            xxi = 0.d0
            xyi = 0.d0
            yyi = 0.d0
            do j = 1, nno
                xxi = xxi+x(i)*zr(ivf+k+i-1)*x(j)*zr(ivf+k+j-1)
                xyi = xyi+x(i)*zr(ivf+k+i-1)*y(j)*zr(ivf+k+j-1)
                yyi = yyi+y(i)*zr(ivf+k+i-1)*y(j)*zr(ivf+k+j-1)
            end do
            matine(1) = matine(1)+poids*yyi
            matine(2) = matine(2)+poids*xyi
            matine(3) = matine(3)+poids*xxi
        end do
    end do
!
    yg = zr(lcastr+2)/volume
    zr(lcastr) = depi*volume*rm
    zr(lcastr+3) = yg
    zr(lcastr+1) = 0.d0
    zr(lcastr+2) = 0.d0
!
!    --- ON DONNE LES INERTIES AU CDG ---
    matine(6) = matine(3)*rm*depi
    matine(1) = matine(1)*rm*depi+matine(6)/2.d0-zr(lcastr)*yg*yg
    matine(2) = 0.d0
    matine(3) = matine(1)
    zr(lcastr+4) = matine(1)
    zr(lcastr+5) = matine(3)
    zr(lcastr+6) = matine(6)
    zr(lcastr+7) = matine(2)
    zr(lcastr+8) = matine(4)
    zr(lcastr+9) = matine(5)
!
end subroutine
