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
subroutine te0109(option, nomte)
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/cq3d2d.h"
#include "asterfort/dfdm2d.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utmess.h"
    character(len=16) :: option, nomte
!
!     BUT:
!       CALCUL DES FLUX DE TEMPERATURE AUX POINTS DE GAUSS
!       ELEMENTS COQUE
!       OPTION : 'FLUX_ELGA'
!
! ---------------------------------------------------------------------
!
!
!
    integer(kind=8) :: nbres
    parameter(nbres=3)
!
    integer(kind=8) :: icodre(nbres)
    integer(kind=8) :: i, kp, itempe, icacoq, imate, iflupg, inbspi
    integer(kind=8) :: ivf, igeom, idfde, ipoids, ndim
    integer(kind=8) :: nno, nnos, npg, jgano, kpg, spt
    integer(kind=8) :: itemps, k, mater, nbcmp, cdec, nbcou, nivc
!
    real(kind=8) :: valres(nbres), conduc, h, ord
    real(kind=8) :: coor2d(14), dfdx(7), dfdy(7), poids, dtdx, dtdy, dtdz
    real(kind=8) :: ts, tm, ti, dtsdx, dtmdx, dtidx, dtsdy, dtmdy, dtidy, px3
    real(kind=8) :: va1a2(3), na1a2, x1, y1, z1, x2, y2, z2, x3, y3, z3
    real(kind=8) :: pvec1(3), pvec2(3), npvec1, fx, fy, fz
    real(kind=8) :: ep, fac1, fac2, fac3
!
    character(len=8) :: fami, poum
    character(len=16) :: nomres(nbres)
    character(len=32) :: phenom
!
!-----------------------------------------------------------------------
!
    valres(1) = 0.d0
    valres(2) = 0.d0
    valres(3) = 0.d0
    fami = 'FPG1'
    kpg = 1
    spt = 1
    poum = '+'
    nbcmp = 9
!
    call jevech('PMATERC', 'L', imate)
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PNBSP_I', 'L', inbspi)
    call jevech('PCACOQU', 'L', icacoq)
    call jevech('PTEMPER', 'L', itempe)
    call jevech('PINSTR', 'L', itemps)
    call jevech('PFLUXPG', 'E', iflupg)
!
    nbcou = zi(inbspi)
    ASSERT(nbcou .eq. 1)
!
! --- RECUPERATION DE LA NATURE DU MATERIAU DANS PHENOM
!     -------------------------------------------------
    mater = zi(imate)
    call rccoma(mater, 'THER', 1, phenom, icodre(1))
!
!       --------------------------
! ----- CAS DES COQUES ISOTROPES :
!       --------------------------
    if (phenom .eq. 'THER') then
!
        nomres(1) = 'LAMBDA'
        call rcvalb(fami, kpg, spt, poum, mater, &
                    ' ', 'THER', 1, 'INST', [zr(itemps)], &
                    1, nomres, valres, icodre, 1)
        conduc = valres(1)
        h = zr(icacoq)/2.d0
        ord = 0.d0
        ep = 2.d0*h
    else
        call utmess('F', 'ELEMENTS3_18', sk=phenom)
    end if
!
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
!
    do i = 1, 3
        va1a2(i) = zr(igeom+i+2)-zr(igeom+i-1)
    end do
    na1a2 = sqrt(va1a2(1)**2+va1a2(2)**2+va1a2(3)**2)
    do i = 1, 3
        va1a2(i) = va1a2(i)/na1a2
    end do
!
    x1 = zr(igeom)
    y1 = zr(igeom+1)
    z1 = zr(igeom+2)
    x2 = zr(igeom+3)
    y2 = zr(igeom+4)
    z2 = zr(igeom+5)
    x3 = zr(igeom+6)
    y3 = zr(igeom+7)
    z3 = zr(igeom+8)
    pvec1(1) = (y2-y1)*(z3-z1)-(z2-z1)*(y3-y1)
    pvec1(2) = (z2-z1)*(x3-x1)-(z3-z1)*(x2-x1)
    pvec1(3) = (x2-x1)*(y3-y1)-(x3-x1)*(y2-y1)
    npvec1 = sqrt(pvec1(1)**2+pvec1(2)**2+pvec1(3)**2)
    do i = 1, 3
        pvec1(i) = pvec1(i)/npvec1
    end do
!
    pvec2(1) = (pvec1(2)*va1a2(3)-pvec1(3)*va1a2(2))
    pvec2(2) = (pvec1(3)*va1a2(1)-pvec1(1)*va1a2(3))
    pvec2(3) = (pvec1(1)*va1a2(2)-pvec1(2)*va1a2(1))
!
    call cq3d2d(nno, zr(igeom), 1.d0, 0.d0, coor2d)
!
    do nivc = -1, 1
!
        if (nivc .lt. 0) then
            px3 = ord-ep/2.d0
            cdec = 3
        else if (nivc .gt. 0) then
            px3 = ord+ep/2.d0
            cdec = 6
        else
            px3 = ord
            cdec = 0
        end if
!
        do kp = 1, npg
            k = (kp-1)*nno
            call dfdm2d(nno, kp, ipoids, idfde, coor2d, &
                        poids, dfdx, dfdy)
            dtmdx = 0.d0
            dtidx = 0.d0
            dtsdx = 0.d0
            dtmdy = 0.d0
            dtidy = 0.d0
            dtsdy = 0.d0
            tm = 0.d0
            ti = 0.d0
            ts = 0.d0
!
            do i = 1, nno
                dtmdx = dtmdx+zr(itempe+3*i-3)*dfdx(i)
                dtmdy = dtmdy+zr(itempe+3*i-3)*dfdy(i)
                dtidx = dtidx+zr(itempe+3*i-2)*dfdx(i)
                dtidy = dtidy+zr(itempe+3*i-2)*dfdy(i)
                dtsdx = dtsdx+zr(itempe+3*i-1)*dfdx(i)
                dtsdy = dtsdy+zr(itempe+3*i-1)*dfdy(i)
                tm = tm+zr(itempe+3*i-3)*zr(ivf+k+i-1)
                ti = ti+zr(itempe+3*i-2)*zr(ivf+k+i-1)
                ts = ts+zr(itempe+3*i-1)*zr(ivf+k+i-1)
            end do
            fac1 = (1.d0-(px3/h)**2)
            fac2 = -px3*(1.d0-px3/h)/(2.d0*h)
            fac3 = px3*(1.d0+px3/h)/(2.d0*h)
            dtdx = dtmdx*fac1+dtidx*fac2+dtsdx*fac3
            dtdy = dtmdy*fac1+dtidy*fac2+dtsdy*fac3
            dtdz = ts*(.5d0+px3/h)/h-2.d0*tm*px3/h**2-ti*(.5d0-px3/h)/h
!
            fx = -conduc*dtdx
            fy = -conduc*dtdy
            fz = -conduc*dtdz
!
            zr(iflupg+(kp-1)*nbcmp-1+cdec+1) = fx*va1a2(1)+fy*pvec2(1)+fz*pvec1(1)
            zr(iflupg+(kp-1)*nbcmp-1+cdec+2) = fx*va1a2(2)+fy*pvec2(2)+fz*pvec1(2)
            zr(iflupg+(kp-1)*nbcmp-1+cdec+3) = fx*va1a2(3)+fy*pvec2(3)+fz*pvec1(3)
!
        end do
!
    end do
!
end subroutine
