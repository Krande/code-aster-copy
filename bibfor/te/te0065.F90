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
subroutine te0065(option, nomte)
!
    implicit none
!
#include "MeshTypes_type.h"
#include "jeveux.h"
#include "asterfort/dfdm3d.h"
#include "asterfort/elref2.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/poutre_modloc.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utmess.h"
!
    character(len=16) :: option, nomte
!     CALCULE DES TERMES PROPRES A UN STRUCTURE
!     OPTION : 'MASS_INER'              (ELEMENTS ISOPARAMETRIQUES 3D)
!     ------------------------------------------------------------------
    integer(kind=8), parameter :: nbres = 3
!-----------------------------------------------------------------------
    integer(kind=8) :: l, lcastr, ndim, nnos
    real(kind=8) :: rho(1), xxi, yyi, zero, zzi
!-----------------------------------------------------------------------
!
    integer(kind=8) :: icodre(nbres)
    character(len=8) :: lielrf(MT_NBFAMX)
    character(len=16) :: nomres(nbres)
    character(len=32) :: phenom
    real(kind=8) :: poids, volume
    real(kind=8) :: x(MT_NNOMAX), y(MT_NNOMAX), z(MT_NNOMAX), xg, yg, zg, matine(6)
    real(kind=8) :: rhopou, rhoflu, tpg, valres(nbres), ayz, ycell, rapp, yf
    integer(kind=8) :: ipoids, ivf, idfde, igeom, nbv, lcorr
    integer(kind=8) :: jgano, nno, kp, npg, i, j, imate, ntrou
!     ------------------------------------------------------------------
    zero = 0.d0
!
    call elref2(nomte, MT_NBFAMX, lielrf, ntrou)
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PMATERC', 'L', imate)
!
    call rccoma(zi(imate), 'ELAS', 1, phenom, icodre(1))
!
    if (lielrf(2) (1:4) .eq. 'POHO') then
!
!        POUR LES ELEMENTS DE LA MODELISATION '3D_FAISCEAU':
!        ===================================================
!
!        - DETERMINATION DU RHO 'POUTRE': RHOPOU
        if (phenom .eq. 'ELAS') then
            nomres(1) = 'RHO'
            nbv = 1
        else
            call utmess('F', 'ELEMENTS3_98')
        end if
        tpg = 0.d0
        call rcvalb('FPG1', 1, 1, '+', zi(imate), &
                    ' ', phenom, 0, ' ', [tpg], &
                    nbv, nomres, valres, icodre, 1)
        rhopou = valres(1)
!
!        - DETERMINATION DU RHO 'FLUIDE': RHOFLU
        call rccoma(zi(imate), 'FLUIDE', 1, phenom, icodre(1))
        if (phenom .eq. 'FLUIDE') then
            nomres(1) = 'RHO'
            nbv = 1
        else
            call utmess('F', 'ELEMENTS3_98')
        end if
        tpg = 0.d0
        call rcvalb('FPG1', 1, 1, '+', zi(imate), &
                    ' ', phenom, 0, ' ', [tpg], &
                    nbv, nomres, valres, icodre, 1)
        rhoflu = valres(1)
!
!        - DETERMINATION DU RHO 'EQUIVALENT' : RHO
!          RHO = ( RHOPOU * AYZ * RAPP ) +  ( RHOFLUI * YF )
!                 RAPP :=  COEF_ECH **2 / A_CELL
!                 YF   :=  A_FLUI  / A_CELL
!                 AYZ  := AIRE_SECTION_POUTRE
        call poutre_modloc('CAGNPO', ['A1'], 1, valeur=ayz)
        call jevech('PCAPOUF', 'L', lcorr)
        ycell = zr(lcorr+4)
        rapp = zr(lcorr+5)
        rapp = rapp*rapp/ycell
        yf = zr(lcorr+3)/ycell
        rho(1) = (rhopou*ayz*rapp)+(rhoflu*yf)
        call elrefe_info(elrefe=lielrf(1), fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, &
                         npg=npg, jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
!
    else
        call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, &
                         npg=npg, jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
        if (phenom .eq. 'ELAS' .or. phenom .eq. 'ELAS_ISTR' .or. phenom .eq. 'ELAS_ORTH') then
            call rcvalb('FPG1', 1, 1, '+', zi(imate), &
                        ' ', phenom, 0, ' ', [0.d0], &
                        1, 'RHO', rho, icodre, 1)
        else
            call utmess('F', 'ELEMENTS_50')
        end if
    end if
!
!
    do i = 1, nno
        x(i) = zr(igeom+3*(i-1))
        y(i) = zr(igeom+3*i-2)
        z(i) = zr(igeom+3*i-1)
    end do
!
    call jevech('PMASSINE', 'E', lcastr)
    do i = 0, 3
        zr(lcastr+i) = zero
    end do
    do i = 1, 6
        matine(i) = zero
    end do
!
!     --- BOUCLE SUR LES POINTS DE GAUSS
    volume = 0.d0
    do kp = 1, npg
        l = (kp-1)*nno
        call dfdm3d(nno, kp, ipoids, idfde, zr(igeom), &
                    poids)
!
        volume = volume+poids
        do i = 1, nno
!           --- CDG ---
            zr(lcastr+1) = zr(lcastr+1)+poids*x(i)*zr(ivf+l+i-1)
            zr(lcastr+2) = zr(lcastr+2)+poids*y(i)*zr(ivf+l+i-1)
            zr(lcastr+3) = zr(lcastr+3)+poids*z(i)*zr(ivf+l+i-1)
!           --- INERTIE ---
            xxi = 0.d0
            yyi = 0.d0
            zzi = 0.d0
            do j = 1, nno
                xxi = xxi+x(i)*zr(ivf+l+i-1)*x(j)*zr(ivf+l+j-1)
                yyi = yyi+y(i)*zr(ivf+l+i-1)*y(j)*zr(ivf+l+j-1)
                zzi = zzi+z(i)*zr(ivf+l+i-1)*z(j)*zr(ivf+l+j-1)
                matine(2) = matine(2)+poids*x(i)*zr(ivf+l+i-1)*y(j)*zr(ivf+l+j-1)
                matine(4) = matine(4)+poids*x(i)*zr(ivf+l+i-1)*z(j)*zr(ivf+l+j-1)
                matine(5) = matine(5)+poids*y(i)*zr(ivf+l+i-1)*z(j)*zr(ivf+l+j-1)
            end do
            matine(1) = matine(1)+poids*(yyi+zzi)
            matine(3) = matine(3)+poids*(xxi+zzi)
            matine(6) = matine(6)+poids*(xxi+yyi)
        end do
    end do
!
    xg = zr(lcastr+1)/volume
    yg = zr(lcastr+2)/volume
    zg = zr(lcastr+3)/volume
    zr(lcastr) = volume*rho(1)
    zr(lcastr+1) = xg
    zr(lcastr+2) = yg
    zr(lcastr+3) = zg
!
!     ---ON DONNE LES INERTIES EN G ---
    zr(lcastr+4) = matine(1)*rho(1)-zr(lcastr)*(yg*yg+zg*zg)
    zr(lcastr+5) = matine(3)*rho(1)-zr(lcastr)*(xg*xg+zg*zg)
    zr(lcastr+6) = matine(6)*rho(1)-zr(lcastr)*(xg*xg+yg*yg)
    zr(lcastr+7) = matine(2)*rho(1)-zr(lcastr)*(xg*yg)
    zr(lcastr+8) = matine(4)*rho(1)-zr(lcastr)*(xg*zg)
    zr(lcastr+9) = matine(5)*rho(1)-zr(lcastr)*(yg*zg)
!
end subroutine
