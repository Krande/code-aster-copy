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
subroutine te0285(option, nomte)
!
    implicit none
!
#include "jeveux.h"
#include "asterc/r8depi.h"
#include "asterfort/assert.h"
#include "asterfort/dfdm2d.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/getDensity.h"
!
    character(len=16) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
!     CALCUL DES OPTIONS: 'MASS_INER' ELEMENTS 2-D AXI D-PLAN, C-PLAN
!                         'CARA_GEOM' ELEMENTS 2-D D-PLAN
!
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nno, kp, npg, i, j, k, lcastr
    integer(kind=8) :: ipoids, ivf, idfde, igeom, imate
    real(kind=8) :: rho, xg, yg, depi
    real(kind=8) :: dfdx(9), dfdy(9), poids, r, x(9), y(9)
    real(kind=8) :: matine(6), xxi, xyi, yyi, volume
    real(kind=8) :: ixrp2, iyrp2, xp(9), yp(9), xpg, ypg
!
! --------------------------------------------------------------------------------------------------
!
    depi = r8depi()
!
    call elrefe_info(fami='MASS', nno=nno, &
                     npg=npg, jpoids=ipoids, jvf=ivf, jdfde=idfde)
!
    if (option .eq. 'MASS_INER') then
        call jevech('PMATERC', 'L', imate)
        call jevech('PMASSINE', 'E', lcastr)
        call getDensity(zi(imate), rho)

    else if (option .eq. 'CARA_GEOM') then
        rho = 1.d0
        call jevech('PCARAGE', 'E', lcastr)

    else
        ASSERT(ASTER_FALSE)
    end if
!
    call jevech('PGEOMER', 'L', igeom)
    do i = 1, nno
        x(i) = zr(igeom-2+2*i)
        y(i) = zr(igeom-1+2*i)
        xp(i) = 0.d0
        yp(i) = 0.d0
    end do
!
    do i = 0, 3
        zr(lcastr+i) = 0.d0
    end do
    matine = 0.d0
!
!     --- BOUCLE SUR LES POINTS DE GAUSS ---
    volume = 0.d0
    do kp = 1, npg
        k = (kp-1)*nno
        call dfdm2d(nno, kp, ipoids, idfde, zr(igeom), &
                    poids, dfdx, dfdy)
        if (lteatt('AXIS', 'OUI')) then
            r = 0.d0
            do i = 1, nno
                r = r+zr(igeom-2+2*i)*zr(ivf+k+i-1)
            end do
            poids = poids*r
        end if
        volume = volume+poids
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
    if (lteatt('AXIS', 'OUI')) then
        xg = 0.d0
        yg = zr(lcastr+2)/volume
        zr(lcastr) = depi*volume*rho
        zr(lcastr+3) = yg
        zr(lcastr+1) = 0.d0
        zr(lcastr+2) = 0.d0
!
!        --- ON DONNE LES INERTIES AU CDG ---
        matine(6) = matine(3)*rho*depi
        matine(1) = matine(1)*rho*depi+matine(6)/2.d0-zr(lcastr)*yg*yg
        matine(2) = 0.d0
        matine(3) = matine(1)
!
    else
        zr(lcastr) = volume*rho
        zr(lcastr+1) = zr(lcastr+1)/volume
        zr(lcastr+2) = zr(lcastr+2)/volume
        zr(lcastr+3) = 0.d0
!
!        --- ON DONNE LES INERTIES AU CDG ---
        xg = zr(lcastr+1)
        yg = zr(lcastr+2)
        matine(1) = matine(1)*rho-zr(lcastr)*yg*yg
        matine(2) = matine(2)*rho-zr(lcastr)*xg*yg
        matine(3) = matine(3)*rho-zr(lcastr)*xg*xg
        matine(6) = matine(1)+matine(3)
    end if
    zr(lcastr+4) = matine(1)
    zr(lcastr+5) = matine(3)
    zr(lcastr+6) = matine(6)
    zr(lcastr+7) = matine(2)
    zr(lcastr+8) = matine(4)
    zr(lcastr+9) = matine(5)
!
    if (option .eq. 'CARA_GEOM') then
!
! --- CALCUL DE IXRP2 = SOMME((Y*(X**2 + Y**2).DS) ET
! --- CALCUL DE IYRP2 = SOMME((X*(X**2 + Y**2).DS) :
!     --------------------------------------------
        do i = 1, nno
            xp(i) = x(i)-xg
            yp(i) = y(i)-yg
        end do
!
        ixrp2 = 0.d0
        iyrp2 = 0.d0
!
        do kp = 1, npg
            k = (kp-1)*nno
            call dfdm2d(nno, kp, ipoids, idfde, zr(igeom), &
                        poids, dfdx, dfdy)
!
            xpg = 0.d0
            ypg = 0.d0
            do i = 1, nno
                xpg = xpg+xp(i)*zr(ivf+k+i-1)
                ypg = ypg+yp(i)*zr(ivf+k+i-1)
            end do
!
            ixrp2 = ixrp2+xpg*(xpg*xpg+ypg*ypg)*poids
            iyrp2 = iyrp2+ypg*(xpg*xpg+ypg*ypg)*poids
!
        end do
!
        zr(lcastr+10) = ixrp2
        zr(lcastr+11) = iyrp2
!
    end if
!
end subroutine
