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

subroutine te0292(option, nomte)
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
!
    character(len=16) :: option, nomte
! person_in_charge: josselin.delmas at edf.fr
!
!     BUT:
!         CALCUL DE L'INDICATEUR D'ERREUR EN QUANTITE D'INTERET
!         SUR UN ELEMENT 2D AVEC LA METHODE DE ZHU-ZIENKIEWICZ.
!         OPTION : 'ERRE_QIZZ'
!
! ......................................................................
!
!
!
!
    integer(kind=8) :: nno, kp, npg1, i, k, nnos, jgano, ndim, mater
    integer(kind=8) :: ipoids, ivf, idfde, igeom, imate, isiefp, isiefd
    integer(kind=8) :: isigp, isigd, ierr, niv, nbcmp
!
    real(kind=8) :: dfdx(9), dfdy(9), dfdz(9), poids, valres(2)
    real(kind=8) :: sigp11, sigp22, sigp33, sigp12, sigp13, sigp23
    real(kind=8) :: sigd11, sigd22, sigd33, sigd12, sigd13, sigd23
    real(kind=8) :: r
    real(kind=8) :: e, nu, nor, norsig, nu0, eest, he
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
    laxi = .false.
    if (lteatt('AXIS', 'OUI')) laxi = .true.
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PMATERC', 'L', imate)
    mater = zi(imate)
    nomres(1) = 'E'
    nomres(2) = 'NU'
    call jevech('PSIEFP_R', 'L', isiefp)
    call jevech('PSIEFD_R', 'L', isiefd)
    call jevech('PSIGMAP', 'L', isigp)
    call jevech('PSIGMAD', 'L', isigd)
    call jevech('PERREUR', 'E', ierr)
!
    norsig = 0.d0
    zr(ierr) = 0.d0
!
!    BOUCLE SUR LES POINTS DE GAUSS
!
    do kp = 1, npg1
        k = (kp-1)*nno
!
        if (ndim .eq. 2) then
            call dfdm2d(nno, kp, ipoids, idfde, zr(igeom), &
                        poids, dfdx, dfdy)
            nbcmp = 2
        else if (ndim .eq. 3) then
            call dfdm3d(nno, kp, ipoids, idfde, zr(igeom), &
                        poids, dfdx, dfdy, dfdz)
            nbcmp = 3
        else
            ASSERT(.false.)
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
        sigp11 = 0.d0
        sigp22 = 0.d0
        sigp33 = 0.d0
        sigp12 = 0.d0
        sigp13 = 0.d0
        sigp23 = 0.d0
        sigd11 = 0.d0
        sigd22 = 0.d0
        sigd33 = 0.d0
        sigd12 = 0.d0
        sigd13 = 0.d0
        sigd23 = 0.d0
!
        do i = 1, nno
            sigp11 = sigp11+zr(isigp-1+nbcmp*(i-1)+1)*zr(ivf+k+i-1)
            sigp22 = sigp22+zr(isigp-1+nbcmp*(i-1)+2)*zr(ivf+k+i-1)
            sigp33 = sigp33+zr(isigp-1+nbcmp*(i-1)+3)*zr(ivf+k+i-1)
            sigp12 = sigp12+zr(isigp-1+nbcmp*(i-1)+4)*zr(ivf+k+i-1)
            sigd11 = sigd11+zr(isigd-1+nbcmp*(i-1)+1)*zr(ivf+k+i-1)
            sigd22 = sigd22+zr(isigd-1+nbcmp*(i-1)+2)*zr(ivf+k+i-1)
            sigd33 = sigd33+zr(isigd-1+nbcmp*(i-1)+3)*zr(ivf+k+i-1)
            sigd12 = sigd12+zr(isigd-1+nbcmp*(i-1)+4)*zr(ivf+k+i-1)
            if (ndim .eq. 3) then
                sigp13 = sigp13+zr(isigp-1+nbcmp*(i-1)+5)*zr(ivf+k+i-1)
                sigp23 = sigp23+zr(isigp-1+nbcmp*(i-1)+6)*zr(ivf+k+i-1)
                sigd13 = sigd13+zr(isigd-1+nbcmp*(i-1)+5)*zr(ivf+k+i-1)
                sigd23 = sigd23+zr(isigd-1+nbcmp*(i-1)+6)*zr(ivf+k+i-1)
            end if
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
        if (ndim .eq. 2) then
            eest = abs( &
                   ( &
                   sigp11-zr( &
                   isiefp-1+4*(kp-1)+1))*(sigd11-zr(isiefd-1+4*(kp-1)+1)))+abs((sigp22-zr(isie&
                   &fp-1+4*(kp-1)+2))*(sigd22-zr(isiefd-1+4*(kp-1)+2)))+abs((sigp33-zr(isiefp&
                   &-1+4*(kp-1)+3))*(sigd33-zr(isiefd-1+4*(kp-1)+3)))+abs((sigp12-zr(isiefp-1+4&
                   &*(kp-1)+4))*(sigd12-zr(isiefd-1+4*(kp-1)+4) &
                   ) &
                   )*(1.d0+nu &
                   )
!
            zr(ierr) = zr(ierr)+eest*poids/e
!
!    NORME DE L' ENERGIE DE LA SOLUTION CALCULEE
!
            nor = zr( &
                  isiefp-1+4*(kp-1)+1)**2+zr(isiefp-1+4*(kp-1)+2)**2+zr(isiefp-1+4*(kp-1)+3)&
                 &**2+(1.d0+nu)*zr(isiefp-1+4*(kp-1)+4 &
                                   )**2
            norsig = norsig+nor*poids/e
!
        else if (ndim .eq. 3) then
            eest = abs( &
                   ( &
                   sigp11-zr( &
                   isiefp-1+6*(kp-1)+1))*(sigd11-zr(isiefd-1+6*(kp-1)+1)))+abs((sigp22-zr(isie&
                   &fp-1+6*(kp-1)+2))*(sigd22-zr(isiefd-1+6*(kp-1)+2)))+abs((sigp33-zr(isiefp&
                   &-1+6*(kp-1)+3))*(sigd33-zr(isiefd-1+6*(kp-1)+3)))+abs((sigp12-zr(isiefp-1+6&
                   &*(kp-1)+4))*(sigd12-zr(isiefd-1+6*(kp-1)+4)))*(1.d0+nu)+abs((sigp13-zr(isi&
                   &efp-1+6*(kp-1)+5))*(sigd13-zr(isiefd-1+6*(kp-1)+5)))*(1.d0+nu)+abs((sigp23&
                   &-zr(isiefp-1+6*(kp-1)+6))*(sigd23-zr(isiefd-1+6*(kp-1)+6) &
                   ) &
                   )*(1.d0+nu &
                   )
!
            zr(ierr) = zr(ierr)+eest*poids/e
!
!    NORME DE L' ENERGIE DE LA SOLUTION CALCULEE
!
            nor = zr( &
                  isiefp-1+6*(kp-1)+1)**2+zr(isiefp-1+6*(kp-1)+2)**2+zr(isiefp-1+6*(kp-1)+3)&
                 &**2+(1.d0+nu)*zr(isiefp-1+6*(kp-1)+4)**2+(1.d0+nu)*zr(isiefp-1+6*(kp-1)+5)*&
                 &*2+(1.d0+nu)*zr(isiefp-1+6*(kp-1)+6 &
                                  )**2
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
