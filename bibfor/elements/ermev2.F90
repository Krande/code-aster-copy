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

subroutine ermev2(nno, igeom, ff, sig, nbcmp, &
                  dfdx, dfdy, poids, poiaxi, dsx, &
                  dsy, norme)
    implicit none
#include "jeveux.h"
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/lteatt.h"
    integer(kind=8) :: nno, igeom, nbcmp
    integer(kind=8) :: poiaxi
    real(kind=8) :: ff(nno), dfdx(nno), dfdy(nno), poids, dsx, dsy, norme, sig(*)
!  ERREUR EN MECANIQUE - TERME VOLUMIQUE - DIMENSION 2
!  **        **                *                     *
! ======================================================================
!
!     BUT:
!         PREMIER TERME DE L'ESTIMATEUR D'ERREUR EN RESIDU EXPLICITE :
!         CALCUL DE LA DIVERGENCE ET DE LA NORME DE SIGMA EN UN POINT
!         DE GAUSS EN 2D. UTILISE POUR UN ELEMENT "CLASSIQUE" OU UN
!         SOUS-ELEMENT ISSU DU DECOUPAGE X-FEM
!
!
!     ARGUMENTS:
!     ----------
!
!      ENTREE :
!-------------
! IN   NNO    : NOMBRE DE NOEUDS DU TYPE_MAILLE
! IN   IGEOM  : ADRESSE DANS ZR DU TABLEAU DES COORDONNEES
! IN   FF     : TABLEAU DES VALEURS DES FONCTIONS DE FORME AU POINT DE
!               GAUSS COURANT
! IN   SIG    : TABLEAU DES CONTRAINTES AUX NOEUDS
! IN   NBCMP  : NOMBRE DE COMPOSANTES
! IN   DFDX   : DERIVEES DES FONCTIONS DE FORME / X
! IN   DFDY   : DERIVEES DES FONCTIONS DE FORME / Y
! IN   POIAXI : EN MODELISATION AXI :
!               =0 : ON NE MULTIPLIE PAS LE POIDS PAR LE RAYON
!               =1 : ON MULTIPLIE LE POIDS PAR LE RAYON
!               EN MODELISATION AUTRE : SANS OBJET
!
!      SORTIE :
!-------------
! OUT  DSX    : PREMIERE COMPOSANTE DE DIVERGENCE SIGMA
! OUT  DSY    : SECONDE COMPOSANTE DE DIVERGENCE SIGMA
! OUT  NORME  : NORME DE SIGMA AU POINT DE GAUSS
!
!
!     ENTREE ET SORTIE :
!----------------------
! IN/OUT   POIDS  : NE SERT QU'A ETRE MULTIPLIE PAR LE RAYON DANS LE
!                   CAS DES MODELISATIONS AXI
!
! ......................................................................
!
!
!
!
    integer(kind=8) :: i
!
    real(kind=8) :: dsig11, dsig12, dsig22, dsig21, spg11, spg22, spg33, spg12
    real(kind=8) :: r, sig11, sig22, sig33, sig12
!
! ----------------------------------------------------------------------
!
    dsig11 = 0.d0
    dsig12 = 0.d0
    dsig22 = 0.d0
    dsig21 = 0.d0
!
    spg11 = 0.d0
    spg22 = 0.d0
    spg33 = 0.d0
    spg12 = 0.d0
!
!====
! 1. MODELISATION AXI
!====
!
    if (lteatt('AXIS', 'OUI')) then
!
        r = 0.d0
        do i = 1, nno
            r = r+zr(igeom-1+2*(i-1)+1)*ff(i)
!
            sig11 = sig(nbcmp*(i-1)+1)
            sig22 = sig(nbcmp*(i-1)+2)
            sig33 = sig(nbcmp*(i-1)+3)
            sig12 = sig(nbcmp*(i-1)+4)
!
            dsig11 = dsig11+sig11*dfdx(i)
            dsig12 = dsig12+sig12*dfdy(i)
            dsig22 = dsig22+sig22*dfdy(i)
            dsig21 = dsig21+sig12*dfdx(i)
!
            spg11 = spg11+sig11*ff(i)
            spg22 = spg22+sig22*ff(i)
            spg33 = spg33+sig33*ff(i)
            spg12 = spg12+sig12*ff(i)
!
        end do
!
!
        ASSERT(abs(r) .gt. r8prem())
!
        dsx = dsig11+dsig12+(1.d0/r)*(spg11-spg33)
        dsy = dsig21+dsig22+(1.d0/r)*spg12
        if (poiaxi .eq. 1) then
            poids = poids*r
        end if
!
!====
! 2. AUTRE MODELISATION
!====
!
    else
!
        do i = 1, nno
            sig11 = sig(nbcmp*(i-1)+1)
            sig22 = sig(nbcmp*(i-1)+2)
            sig33 = sig(nbcmp*(i-1)+3)
            sig12 = sig(nbcmp*(i-1)+4)
!
            dsig11 = dsig11+sig11*dfdx(i)
            dsig12 = dsig12+sig12*dfdy(i)
            dsig22 = dsig22+sig22*dfdy(i)
            dsig21 = dsig21+sig12*dfdx(i)
!
            spg11 = spg11+sig11*ff(i)
            spg22 = spg22+sig22*ff(i)
            spg33 = spg33+sig33*ff(i)
            spg12 = spg12+sig12*ff(i)
!
        end do
!
        dsx = dsig11+dsig12
        dsy = dsig21+dsig22
!
    end if
!
!====
! 3.
!====
!
    norme = spg11**2+spg22**2+spg33**2+spg12**2
!
end subroutine
