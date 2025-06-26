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
subroutine mefrec(ndim, nbcyl, nbgrp, numgrp, xint, &
                  yint, rint, sgn, orig, beta, &
                  ppxx, ppxy, ppyx, ppyy, vnxx, &
                  vnxy, vnyx, vnyy, tmp)
    implicit none
!
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/mefasr.h"
#include "asterfort/mtcrog.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
    integer(kind=8) :: nbcyl, nbgrp, ndim(14), numgrp(*), sgn(*), orig(*)
    real(kind=8) :: rint(*), xint(*), yint(*), beta(*)
    real(kind=8) :: ppxx(nbcyl, nbgrp), ppxy(nbcyl, nbgrp)
    real(kind=8) :: ppyx(nbcyl, nbgrp), ppyy(nbcyl, nbgrp)
    real(kind=8) :: vnxx(nbcyl, nbgrp), vnxy(nbcyl, nbgrp)
    real(kind=8) :: vnyx(nbcyl, nbgrp), vnyy(nbcyl, nbgrp)
!     ASSEMBLAGE ET CALCUL DES COEFFICIENTS INTERVENANT DANS
!     L EXPRESSION DES FORCES DE PRESSION PERTURBEE, ET ET DES FORCES
!     NORMALES DE FROTTEMENTS SUR CHAQUE CYLINDRE DANS LE CAS D UNE
!     ENCEINTE RECTANGULAIRE
!     OPERATEUR APPELANT : OP0144 , FLUST3, MEFIST
! ----------------------------------------------------------------------
!     OPTION DE CALCUL   : CALC_FLUI_STRU , CALCUL DES PARAMETRES DE
!     COUPLAGE FLUIDE-STRUCTURE POUR UNE CONFIGURATION DE TYPE "FAISCEAU
!     DE TUBES SOUS ECOULEMENT AXIAL"
! ----------------------------------------------------------------------
! IN  : NDIM   : TABLEAU DES DIMENSIONS
! IN  : NBCYL  : NOMBRE DE CYLINDRES
! IN  : NUMGRP : INDICES DES GROUPES D EQUIVALENCE
! IN  : XINT   : COORDONNEES 'X' DES CENTRES DES CYLINDRES DANS
!                LE REPERE AXIAL
! IN  : YINT   : COORDONNEES 'Y' DES CENTRES DES CYLINDRES DANS
!                LE REPERE AXIAL
! IN  : RINT   : RAYONS DES CYLINDRES
! IN  : SGN    : -1 OU +1, COEFFICIENT INTERVENANT DANS LA DECOMPOSITION
!                EN SERIE DE LAURENT, SELON LE NIVEAU D IMAGE
! IN  : ORIG   : NUMERO DU CYLINDRE D ORIGINE DES CYLINDRES REELS OU
!                IMAGES
! IN  : BETA   : ANGLE CUMULE INTERVENANT DANS LA DECOMPOSITION EN
!                SERIE DE LAURENT, POUR LES CYLINDRES IMAGES
! OUT : PPXX   : COEFFICIENT DE MASSES AJOUTEES INTERVENANT DANS LES
!                EFFORTS DE PRESSION PERTURBES SUIVANT XX
! OUT : PPXY   : COEFFICIENT DE MASSES AJOUTEES INTERVENANT DANS LES
!                EFFORTS DE PRESSION PERTURBES SUIVANT XY
! OUT : PPYX   : COEFFICIENT DE MASSES AJOUTEES INTERVENANT DANS LES
!                EFFORTS DE PRESSION PERTURBES SUIVANT YX
! OUT : PPYY   : COEFFICIENT DE MASSES AJOUTEES INTERVENANT DANS LES
!                EFFORTS DE PRESSION PERTURBES SUIVANT YY
! OUT : VNXX   : COEFFICIENT INTERVENANT DANS L EXPRESSION DES EFFORTS
!                VISQUEUX NORMAUX SUIVANT XX
! OUT : VNXY   : COEFFICIENT INTERVENANT DANS L EXPRESSION DES EFFORTS
!                VISQUEUX NORMAUX SUIVANT XY
! OUT : VNYX   : COEFFICIENT INTERVENANT DANS L EXPRESSION DES EFFORTS
!                VISQUEUX NORMAUX SUIVANT YX
! OUT : VNYY   : COEFFICIENT INTERVENANT DANS L EXPRESSION DES EFFORTS
!                VISQUEUX NORMAUX SUIVANT YY
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
    integer(kind=8) :: i, j, k
    integer(kind=8) :: ncyl
    real(kind=8) :: tmp(4, *), rayoi, rayoj
! ----------------------------------------------------------------------
!
! --- LECTURE DES DIMENSIONS
!-----------------------------------------------------------------------
    integer(kind=8) :: ia, ib, idir, ier, igrp, itrav, ix
    integer(kind=8) :: ixx, nbtron, nmax, nv
!-----------------------------------------------------------------------
    nbcyl = ndim(3)
    nbgrp = ndim(4)
    nbtron = ndim(5)
!
    call jemarq()
!
!
! --- TABLEAUX DE TRAVAIL, ALLOCATION MEMOIRE
    nv = 2*nbtron*nbcyl
    call wkvect('&&MEFREC.TMP.AB', 'V V R', nv*(3+2*nbgrp+nv), ia)
    ib = ia+nv*nv
    ixx = ib+nv
    ix = ixx+nv*2*nbgrp
    itrav = ix+nv
!
! --- INITIALISATIONS
!
    nmax = 4*nbtron*nbcyl*nbgrp
    do i = 1, nmax
        zr(ixx+i-1) = 0.d0
    end do
!
    nmax = 2*nbtron*nbcyl
    do igrp = 1, nbgrp
        do idir = 1, 2
            do i = 1, nmax*nmax
                zr(ia+i-1) = 0.d0
            end do
            do i = 1, nmax
                zr(ib+i-1) = 0.d0
                zr(ix+i-1) = 0.d0
            end do
!
! ---       ASSEMBLAGE
            call mefasr(ndim, nbcyl, nbgrp, nbtron, numgrp, &
                        idir, igrp, xint, yint, rint, &
                        sgn, orig, beta, zr(ia), zr(ib))
!
! ---       RESOLUTION DU SYSTEME A.X = B PAR LA METHODE DE CROUT
            ier = 1
            call mtcrog(zr(ia), zr(ib), nmax, nmax, 1, &
                        zr(ix), zr(itrav), ier)
            if (ier .eq. 1) then
                call utmess('F', 'ALGELINE_76')
            end if
!
            do i = 1, nmax
                zr(ixx+i-1+nmax*(2*igrp-idir)) = zr(ix+i-1)
            end do
        end do
    end do
!
!
! --- CALCUL DES COEFFICIENTS PPXX, PPXY, PPYX, PPYY,
! --- ET VNXX, VNXY, VNYX, VNYY
    do i = 1, nbgrp
        do j = 1, nbcyl
            ppxx(j, i) = 2.d0*zr(ixx-1+2*nbtron*(nbcyl*(2*i-1)+j-1)+1)
            ppxy(j, i) = 2.d0*zr(ixx-1+2*nbtron*(nbcyl*(2*i-2)+j-1)+1)
            ppyx(j, i) = 2.d0*zr(ixx-1+2*nbtron*(nbcyl*(2*i-1)+j-1)+2)
            ppyy(j, i) = 2.d0*zr(ixx-1+2*nbtron*(nbcyl*(2*i-2)+j-1)+2)
!
            if (i .eq. numgrp(j)) then
                ppxx(j, i) = ppxx(j, i)+1.d0
                ppyy(j, i) = ppyy(j, i)+1.d0
            end if
        end do
    end do
!
    do i = 1, nbgrp
        do j = 1, nbgrp
            tmp(1, j) = 0.d0
            tmp(2, j) = 0.d0
            tmp(3, j) = 0.d0
            tmp(4, j) = 0.d0
        end do
        do j = 1, nbcyl
            tmp(1, numgrp(j)) = tmp(1, numgrp(j))+ppxx(j, i)
            tmp(2, numgrp(j)) = tmp(2, numgrp(j))+ppxy(j, i)
            tmp(3, numgrp(j)) = tmp(3, numgrp(j))+ppyx(j, i)
            tmp(4, numgrp(j)) = tmp(4, numgrp(j))+ppyy(j, i)
        end do
        do j = 1, nbgrp
            ppxx(j, i) = tmp(1, j)
            ppxy(j, i) = tmp(2, j)
            ppyx(j, i) = tmp(3, j)
            ppyy(j, i) = tmp(4, j)
        end do
    end do
!
! --- ON FORCE LA SYMETRIE DES COEFFICIENTS
!
    do i = 1, nbgrp
        do k = 1, nbcyl
            if (numgrp(k) .eq. i) then
                rayoi = rint(k)
            end if
        end do
        ppxy(i, i) = ppxy(i, i)+ppyx(i, i)
        ppxy(i, i) = ppxy(i, i)/2.d0
        ppyx(i, i) = ppxy(i, i)
        do j = 1, i-1
            do k = 1, nbcyl
                if (numgrp(k) .eq. j) then
                    rayoj = rint(k)
                end if
            end do
            ppxx(i, j) = rayoi*rayoi*ppxx(i, j)+rayoj*rayoj*ppxx(j, i)
            ppxx(i, j) = ppxx(i, j)/2.d0
            ppxx(j, i) = ppxx(i, j)
            ppxx(i, j) = ppxx(i, j)/rayoi/rayoi
            ppxx(j, i) = ppxx(j, i)/rayoj/rayoj
!
            ppxy(i, j) = rayoi*rayoi*ppxy(i, j)+rayoj*rayoj*ppyx(j, i)
            ppxy(i, j) = ppxy(i, j)/2.d0
            ppyx(j, i) = ppxy(i, j)
            ppxy(i, j) = ppxy(i, j)/rayoi/rayoi
            ppyx(j, i) = ppyx(j, i)/rayoj/rayoj
!
            ppyx(i, j) = rayoi*rayoi*ppyx(i, j)+rayoj*rayoj*ppxy(j, i)
            ppyx(i, j) = ppyx(i, j)/2.d0
            ppxy(j, i) = ppyx(i, j)
            ppyx(i, j) = ppyx(i, j)/rayoi/rayoi
            ppxy(j, i) = ppxy(j, i)/rayoj/rayoj
!
            ppyy(i, j) = rayoi*rayoi*ppyy(i, j)+rayoj*rayoj*ppyy(j, i)
            ppyy(i, j) = ppyy(i, j)/2.d0
            ppyy(j, i) = ppyy(i, j)
            ppyy(i, j) = ppyy(i, j)/rayoi/rayoi
            ppyy(j, i) = ppyy(j, i)/rayoj/rayoj
        end do
    end do
!
    do i = 1, nbgrp
        do j = 1, nbgrp
            vnxx(j, i) = 0.5d0*ppxx(j, i)
            vnxy(j, i) = 0.5d0*ppxy(j, i)
            vnyx(j, i) = 0.5d0*ppyx(j, i)
            vnyy(j, i) = 0.5d0*ppyy(j, i)
            if (j .eq. i) then
                ncyl = 0
                do k = 1, nbcyl
                    if (numgrp(k) .eq. i) ncyl = ncyl+1
                end do
                vnxx(j, i) = vnxx(j, i)+0.5d0*ncyl
                vnyy(j, i) = vnyy(j, i)+0.5d0*ncyl
            end if
        end do
    end do
!
! --- MENAGE
    call jedetr('&&MEFREC.TMP.AB')
    call jedema()
end subroutine
