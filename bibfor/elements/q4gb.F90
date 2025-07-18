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

subroutine q4gb(caraq4, xyzl, igau, jacgau, bmat)
    implicit none
#include "jeveux.h"
#include "asterfort/bcoqaf.h"
#include "asterfort/dsqbfb.h"
#include "asterfort/dxqbm.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jquad4.h"
#include "asterfort/q4gbc.h"
    integer(kind=8) :: igau
    real(kind=8) :: xyzl(3, 1), caraq4(*), bmat(8, 1), jacgau
! --- CALCUL DE LA MATRICE (B) RELIANT LES DEFORMATIONS DU PREMIER
! --- ORDRE AUX DEPLACEMENTS AU POINT D'INTEGRATION D'INDICE IGAU
! --- POUR UN ELEMENT DE TYPE Q4G
! --- (I.E. (EPS_1) = (B)*(UN))
! --- D'AUTRE_PART, ON CALCULE LE PRODUIT NOTE JACGAU = JACOBIEN*POIDS
!     ------------------------------------------------------------------
!     IN  XYZL(3,NNO)   : COORDONNEES DES CONNECTIVITES DE L'ELEMENT
!                         DANS LE REPERE LOCAL DE L'ELEMENT
!     IN  IGAU          : INDICE DU POINT D'INTEGRATION
!     OUT JACGAU        : PRODUIT JACOBIEN*POIDS AU POINT D'INTEGRATION
!                         COURANT
!     OUT BMAT(8,1)     : MATRICE (B) AU POINT D'INTEGRATION COURANT
    integer(kind=8) :: ndim, nno, nnos, npg, ipoids, icoopg, ivf, idfdx, idfd2, jgano
    real(kind=8) :: bm(3, 8), bf(3, 12), bc(2, 12), qsi, eta, jacob(5)
! ------------------------------------------------------------------
!
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, &
                     npg=npg, jpoids=ipoids, jcoopg=icoopg, jvf=ivf, jdfde=idfdx, &
                     jdfd2=idfd2, jgano=jgano)
!
! --- COORDONNEES DU POINT D'INTEGRATION COURANT :
!     ------------------------------------------
    qsi = zr(icoopg-1+ndim*(igau-1)+1)
    eta = zr(icoopg-1+ndim*(igau-1)+2)
!
! --- CALCUL DE LA MATRICE JACOBIENNE ET DE SON DETERMINANT AU POINT
! --- D'INTEGRATION D'INDICE IGAU
!     ---------------------------
    call jquad4(xyzl, qsi, eta, jacob)
!
! --- PRODUIT JACOBIEN*POIDS
!     ----------------------
    jacgau = zr(ipoids+igau-1)*jacob(1)
!
! --- CALCUL DE LA MATRICE B_MEMBRANE NOTEE, ICI, (BM)
!     ------------------------------------------------
    call dxqbm(qsi, eta, jacob(2), bm)
!
! --- CALCUL DE LA MATRICE B_FLEXION RELATIVE AUX INCONNUES W, BETAX
! --- ET BETAY.
! --- CETTE MATRICE EST IDENTIQUE A LA PARTIE DE LA MATRICE B_FLEXION
! --- DU DSQ OU L'ON NE TIENT COMPTE QUE DE L'INTERPOLATION LINEAIRE
! --- DES INCONNUES W, BETAX ET BETAY, CETTE MATRICE EST NOTEE (BFB)
! --- POUR LE DSQ.
! --- ON RAPPELLE QUE LE Q4G EST UN ELEMENT ISOPARAMETRIQUE EXPRIME
! --- EN COORDONNEES CARTESIENNES.
!     ----------------------------
    call dsqbfb(qsi, eta, jacob(2), bf)
!
! --- CALCUL DE LA MATRICE B_CISAILLEMENT, NOTEE (BC)
!     -----------------------------------------------
    call q4gbc(qsi, eta, jacob(2), caraq4, bc)
!
! --- AFFECTATION DE LA MATRICE B COMPLETE, NOTEE (BMAT)
! --- AVEC LES MATRICES (BM), (BF) ET (BC)
!     ------------------------------------
    call bcoqaf(bm, bf, bc, nno, bmat)
!
end subroutine
