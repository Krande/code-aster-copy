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
subroutine mbgnlr(lVect, lMatr, &
                  nno, ncomp, imate, icompo, &
                  dff, alpha, beta, h, &
                  preten, igeom, ideplm, ideplp, &
                  kpg, fami, ipoids, icontp, &
                  ivectu, imatuu)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "jeveux.h"
#include "asterfort/mbhenh.h"
#include "asterfort/mbhesv.h"
#include "asterfort/mbpk2c.h"
#include "asterfort/mbtgin.h"
#include "asterfort/mbvfie.h"
#include "asterfort/subaco.h"
#include "asterfort/sumetr.h"
#include "asterfort/subacv.h"
!
    aster_logical, intent(in) :: lVect, lMatr
    character(len=4) :: fami
    integer(kind=8) :: nno, ncomp, kpg
    integer(kind=8) :: imate, icompo, igeom, ideplm, ideplp, ipoids, icontp, ivectu
    integer(kind=8) :: imatuu
    real(kind=8) :: dff(2, nno), alpha, beta, h, preten
!
! ----------------------------------------------------------------------
!    - FONCTION REALISEE:  CALCUL DES OPTIONS DE COMPORTEMENT :
!                            - FULL_MECA
!                            - RAPH_MECA
!                            - RIGI_MECA_TANG
!                          POUR LES MEMBRANES EN GRANDES DEFORMATIONS
! ----------------------------------------------------------------------
! IN  NOMTE             NOM DU TYPE ELEMENT
!     VECTEU            BOOL: 1 SI FULL_MECA OU RAPH_MECA
!     MATRIC            BOOL: 1 SI FULL_MECA OU RIGI_MECA
!     NNO               NOMBRE DE NOEUDS
!     NCOMP             NOMBRE DE COMPOSANTS DANS LES VECTEURS COLONNES
!                       DE CONTRAINTE ET DEFORMATION
!     IMATE             ADRESSE DANS ZI DU TABLEAU PMATERC
!     ICOMPO            ADRESSE DANS ZK16 DU TABLEAU PCOMPOR
!     DFF               DERIVEE DES F. DE FORME
!     ALPHA, BETA       ANGLES DEF. LA BASE DE L'ÉCRITURE DES CONTRAINTES
!     H                 EPAISSEUR DE LA MEMBRANE
!     PRETEN            PRECONTRAINTES
!     IGEOM             ADRESSE DANS ZR DU TABLEAU PGEOMER
!     IDEPLM            ADRESSE DANS ZR DU TABLEAU PDEPLMR
!     IDEPLP            ADRESSE DANS ZR DU TABLEAU PDEPLPR
!     KPG               NUMERO DU POINT DE GAUSS DANS LA BOUCLE
!     FAMI              NOM DE LA FAMILLE DE POINTS DE GAUSS :
!                       'RIGI','MASS',..
!     IPOIDS            ADRESSE DANS ZR DU TABLEAU POIDS
!     ICONTP            ADRESSE DANS ZR DU TABLEAU PCONTPR
!     IVECTU            ADRESSE DANS ZR DU TABLEAU PMATUUR
!     IMATUU            ADRESSE DANS ZR DU TABLEAU PMATUUR
!
! OUT ***          ***
! ----------------------------------------------------------------------

    integer(kind=8) :: n, nn, c, incm
    real(kind=8) :: posdef(3*nno)
    real(kind=8) :: covaini(3, 3), metrini(2, 2), jacini, cnvaini(3, 2), aini(2, 2)
    real(kind=8) :: covadef(3, 3), metrdef(2, 2), jacdef, cnvadef(3, 2), adef(2, 2)
    real(kind=8) :: sigpk2(2, 2), dsigpk2(2, 2, 2, 2), sighca(3), sigpk2temp(3)
    real(kind=8) :: ktgt(3*nno, 3*nno)
    real(kind=8) :: vecfie(3*nno)

! - CALCUL DES COORDONNEES COVARIANTES ET CONTRAVARIANTES DE LA SURFACE INITIALE
    call subaco(nno, dff, zr(igeom), covaini)
    call sumetr(covaini, metrini, jacini)
    call subacv(covaini, metrini, jacini, cnvaini, aini)

! - CALCUL DES COORDONNEES COVARIANTES ET CONTRAVARIANTES DE LA SURFACE DEFORMEE
!
    do n = 1, 3*nno
        posdef(n) = zr(igeom+n-1)+zr(ideplm+n-1)+zr(ideplp+n-1)
    end do

    call subaco(nno, dff, posdef, covadef)
    call sumetr(covadef, metrdef, jacdef)

    call subacv(covadef, metrdef, jacdef, cnvadef, adef)

! - ON APPELLE LA LDC HYPERELASTIQUE NEO-HOOKEENE
! - ON OBTIENT LES CONTRAINTES A L'ITERATION DE NEWTON (i-1) (SIGPK2: TENSEUR SYM)
!
    if (zk16(icompo) (1:16) .eq. 'ELAS_MEMBRANE_SV') then
        call mbhesv(imate, kpg, fami, aini, metrini, metrdef, sigpk2, dsigpk2)
    elseif (zk16(icompo) (1:16) .eq. 'ELAS_MEMBRANE_NH') then
        call mbhenh(imate, kpg, fami, aini, adef, jacini, jacdef, sigpk2, dsigpk2)
    else
        ASSERT(.false.)
    end if

! - SI LA NORME EUCLIDIENNE DE SIGPK2 EST NULLE, ON APPLIQUE DES PRECONTRAINTES

    if (sqrt(sigpk2(1, 1)**2+2*sigpk2(1, 2)**2+sigpk2(2, 2)**2) .lt. 1.0d-6) then
        sigpk2(1, 1) = sigpk2(1, 1)+preten
        sigpk2(2, 2) = sigpk2(2, 2)+preten
    end if

! - ON CALCUL LA MATRICE TANGENTE ELEMENTAIRE DUE AUX EFFORTS INTERNES
!
    if (lMatr) then
        call mbtgin(nno, kpg, dff, sigpk2, dsigpk2, ipoids, h, covadef, ktgt)
    end if

    if (lVect) then

! ---   ON EN DEDUIT LES CONTRAINTES INTEGREES (SUR L'EPAISSEUR) DE CAUCHY
! ---   SIGMA_CAUCHY = (SIGMA11, SIGMA22, SIGMA12)
        sigpk2temp(1) = sigpk2(1, 1)
        sigpk2temp(2) = sigpk2(2, 2)
        sigpk2temp(3) = sigpk2(1, 2)

        call mbpk2c(0, alpha, beta, h, covaini, jacini, jacdef, sigpk2temp, sighca)

! ---   CALCUL DU VECTEUR FORCE INTERNE ELEMENTAIRE
!
        call mbvfie(nno, kpg, dff, sigpk2, ipoids, h, covadef, vecfie)

    end if

! - RANGEMENT DES RESULTATS
!
    if (lVect) then
        do n = 1, 3*nno
            zr(ivectu+n-1) = zr(ivectu+n-1)+vecfie(n)*jacini
        end do
        do c = 1, ncomp
            zr(icontp+(kpg-1)*ncomp+c-1) = sighca(c)
        end do
    end if

    if (lMatr) then
        incm = 0
        do n = 1, 3*nno
            do nn = 1, n
                incm = incm+1
                zr(imatuu+incm-1) = zr(imatuu+incm-1)+ktgt(n, nn)*jacini
            end do
        end do
    end if

end subroutine
