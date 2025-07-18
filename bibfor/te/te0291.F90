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

subroutine te0291(option, nomte)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dfdm2d.h"
#include "asterfort/dfdm3d.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/rcvalb.h"
#include "asterfort/uthk.h"
    character(len=16) :: option, nomte
! person_in_charge: josselin.delmas at edf.fr
!
!     BUT:
!         CALCUL DE L'INDICATEUR D'ERREUR EN ENERGIE
!         SUR UN ELEMENT AVEC LA METHODE DE ZHU-ZIENKIEWICZ.
!         OPTION : 'CALC_ESTI_ERRE'
!
! ......................................................................
!
!
!
!
    integer(kind=8) :: nno, kp, npg1, i, k, nnos, jgano, ndim
    integer(kind=8) :: ipoids, ivf, idfde, igeom, niv, nbcmp
    integer(kind=8) :: ierr, imate, isigm, isigno, mater
!
    real(kind=8) :: dfdx(27), dfdy(27), dfdz(27), poids, valres(2)
    real(kind=8) :: sigl11, sigl22, sigl33, sigl12, sigl13, sigl23
    real(kind=8) :: sigc11, sigc22, sigc33, sigc12, sigc13, sigc23
    real(kind=8) :: esig11, esig22, esig33, esig12, esig13, esig23
    real(kind=8) :: e, nu, eest, nor, norsig, nu0, he, r
!
    integer(kind=8) :: icodre(2)
    character(len=4) :: fami
    character(len=16) :: nomres(2)
!
    aster_logical :: laxi
!
! ----------------------------------------------------------------------
!
    fami = 'RIGI'
    call elrefe_info(fami=fami, ndim=ndim, nno=nno, nnos=nnos, npg=npg1, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PMATERC', 'L', imate)
    mater = zi(imate)
    nomres(1) = 'E'
    nomres(2) = 'NU'
!
!     CHAMP DE CONTRAINTE CALCULE
    call jevech('PSIEF_R', 'L', isigm)
!     CHAMP DE CONTRAINTE LISSE
    call jevech('PSIGMA', 'L', isigno)
!
    call jevech('PERREUR', 'E', ierr)
!
    norsig = 0.d0
    zr(ierr) = 0.d0
!
!    BOUCLE SUR LES POINTS DE GAUSS
!
    laxi = .false.
    if (lteatt('AXIS', 'OUI')) laxi = .true.
!
    do kp = 1, npg1
        k = (kp-1)*nno
        if (ndim .eq. 2) then
            call dfdm2d(nno, kp, ipoids, idfde, zr(igeom), &
                        poids, dfdx, dfdy)
            nbcmp = 4
        else if (ndim .eq. 3) then
            call dfdm3d(nno, kp, ipoids, idfde, zr(igeom), &
                        poids, dfdx, dfdy, dfdz)
            nbcmp = 6
        end if
!
        if (laxi) then
            r = 0.d0
            do i = 1, nno
                r = r+zr(igeom+2*(i-1))*zr(ivf+k+i-1)
            end do
            poids = poids*r
        end if
!
        sigl11 = 0.d0
        sigl22 = 0.d0
        sigl33 = 0.d0
        sigl12 = 0.d0
        sigl13 = 0.d0
        sigl23 = 0.d0
        do i = 1, nno
            sigl11 = sigl11+zr(isigno-1+nbcmp*(i-1)+1)*zr(ivf+k+i-1)
            sigl22 = sigl22+zr(isigno-1+nbcmp*(i-1)+2)*zr(ivf+k+i-1)
            sigl33 = sigl33+zr(isigno-1+nbcmp*(i-1)+3)*zr(ivf+k+i-1)
            sigl12 = sigl12+zr(isigno-1+nbcmp*(i-1)+4)*zr(ivf+k+i-1)
            if (ndim .eq. 3) then
                sigl13 = sigl13+zr(isigno-1+nbcmp*(i-1)+5)*zr(ivf+k+i-1)
                sigl23 = sigl23+zr(isigno-1+nbcmp*(i-1)+6)*zr(ivf+k+i-1)
            end if
!
        end do
!
        call rcvalb(fami, kp, 1, '+', mater, &
                    ' ', 'ELAS', 0, ' ', [0.d0], &
                    2, nomres, valres, icodre, 1)
        e = valres(1)
        nu = valres(2)
!
!    ESTIMATION DE L'ERREUR EN NORME DE L' ENERGIE
!
        sigc11 = zr(isigm-1+nbcmp*(kp-1)+1)
        sigc22 = zr(isigm-1+nbcmp*(kp-1)+2)
        sigc33 = zr(isigm-1+nbcmp*(kp-1)+3)
        sigc12 = zr(isigm-1+nbcmp*(kp-1)+4)
!
        esig11 = sigl11-sigc11
        esig22 = sigl22-sigc22
        esig33 = sigl33-sigc33
        esig12 = sigl12-sigc12
!
        if (ndim .eq. 2) then
            eest = esig11**2+esig22**2+esig33**2+2*(1.d0+nu)*(esig12)**2-2*nu*esig11*esig22-2*&
                   &nu*esig11*esig33-2*nu*esig22*esig33
            zr(ierr) = zr(ierr)+eest*poids/e
!
!    NORME DE L' ENERGIE DE LA SOLUTION CALCULEE
!
            nor = sigc11**2+sigc22**2+sigc33**2+2*(1.d0+nu)*(sigc12)**2-2*nu*sigc11*sigc22-2*n&
                  &u*sigc11*sigc33-2*nu*sigc22*sigc33
            norsig = norsig+nor*poids/e
!
        else if (ndim .eq. 3) then
!
            sigc13 = zr(isigm-1+nbcmp*(kp-1)+5)
            sigc23 = zr(isigm-1+nbcmp*(kp-1)+6)
!
            esig13 = sigl13-sigc13
            esig23 = sigl23-sigc23
!
            eest = esig11**2+esig22**2+esig33**2+2*(1.d0+nu)*(esig12)**2+2*(1.d0+nu)*(esig13)*&
                   &*2+2*(1.d0+nu)*(esig23)**2-2*nu*esig11*esig22-2*nu*esig11*esig33-2*nu*esig2&
                   &2*esig33
            zr(ierr) = zr(ierr)+eest*poids/e
!
!    NORME DE L' ENERGIE DE LA SOLUTION CALCULEE
!
            nor = sigc11**2+sigc22**2+sigc33**2+2*(1.d0+nu)*(sigc12)**2+2*(1.d0+nu)*(sigc13)**&
                  &2+2*(1.d0+nu)*(sigc23)**2-2*nu*sigc11*sigc22-2*nu*sigc11*sigc33-2*nu*sigc22*&
                  &sigc33
            norsig = norsig+nor*poids/e
        else
            ASSERT(.false.)
        end if
!
    end do
!
    niv = 1
    call uthk(nomte, zr(igeom), he, ndim, niv)
!
    if ((zr(ierr)+norsig) .ne. 0.d0) then
        nu0 = 100.d0*sqrt(zr(ierr)/(zr(ierr)+norsig))
    else
        nu0 = 0.d0
    end if
!
    zr(ierr) = sqrt(zr(ierr))
    zr(ierr+1) = nu0
    zr(ierr+2) = sqrt(norsig)
    zr(ierr-1+10) = he
!
!
end subroutine
