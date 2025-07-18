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

subroutine te0342(option, nomte)
!
!
! --------------------------------------------------------------------------------------------------
!     CALCUL
!       - DU VECTEUR ELEMENTAIRE EFFORT GENERALISE,
!     POUR LES ELEMENTS DE POUTRE DE TIMOSHENKO AVEC GAUCHISSEMENT.
!     ------------------------------------------------------------------
! IN  OPTION : K16 : NOM DE L'OPTION A CALCULER
!        'SIEF_ELGA'
! IN  NOMTE  : K16 : NOM DU TYPE ELEMENT
!        'MECA_POU_D_TG': POUTRE DROITE DE TIMOSHENKO AVEC GAUCHISSEMENT
!
!
! --------------------------------------------------------------------------------------------------
!
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/jevech.h"
#include "asterfort/jsd1ff.h"
#include "asterfort/lonele.h"
#include "asterfort/matrot.h"
#include "asterfort/moytem.h"
#include "asterfort/poutre_modloc.h"
#include "asterfort/r8inir.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utpvgl.h"
!
    character(len=*) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: i, igau, j, jdepl, jeffo, k, lmater, iret
    integer(kind=8) :: lorien, nbpar, nbres, nc, nno, npg
    real(kind=8) :: a, alfay, alfaz, e, g, nu
    real(kind=8) :: phiy, phiz, valpar, xiy, xiz, xjg, xjx, xl, xl2
! --------------------------------------------------------------------------------------------------
    real(kind=8) :: b(7, 14)
    real(kind=8) :: pgl(14, 14), depl(14), depglo(14), epsgen(7), siggen(3, 7)
    character(len=8) :: nompar
! --------------------------------------------------------------------------------------------------
    parameter(nbres=2)
    integer(kind=8) :: codres(nbres)
    real(kind=8) :: valres(nbres)
    character(len=16) :: nomres(nbres)
    data nomres/'E', 'NU'/
! --------------------------------------------------------------------------------------------------
    integer(kind=8), parameter :: nb_cara = 7
    real(kind=8) :: vale_cara(nb_cara)
    character(len=8) :: noms_cara(nb_cara)
    data noms_cara/'A1', 'IY1', 'IZ1', 'AY1', 'AZ1', 'JX1', 'JG1'/
! --------------------------------------------------------------------------------------------------
!
    call r8inir(7*14, 0.0d0, b, 1)
    call r8inir(7*3, 0.0d0, siggen, 1)
!
    nbpar = 0
    nompar = '  '
    valpar = 0.d0
    valres(:) = 0.0d0
!
!   recuperation de la temperature :
    npg = 3
    call moytem('RIGI', npg, 1, '+', valpar, iret)
!
    nbpar = 1
    nompar = 'TEMP'
!
!   recuperation et interpolation des caracteristiques materiaux
    call jevech('PMATERC', 'L', lmater)
    call rcvalb('RIGI', npg, 1, '+', zi(lmater), ' ', 'ELAS', nbpar, nompar, [valpar], &
                nbres, nomres, valres, codres, 1)
!
    e = valres(1)
    nu = valres(2)
    g = e/(2.0d0*(1.0d0+nu))
!
!   recuperation des caracteristiques generales des sections :
    call poutre_modloc('CAGNPO', noms_cara, nb_cara, lvaleur=vale_cara)
!
    a = vale_cara(1)
    xiy = vale_cara(2)
    xiz = vale_cara(3)
    alfay = vale_cara(4)
    alfaz = vale_cara(5)
    xjx = vale_cara(6)
    xjg = vale_cara(7)
!
    nno = 2
    nc = 7
!
!   Calcul de la longueur de la poutre
    xl = lonele()
    xl2 = xl*xl
!
!   calcul des coefficients d'influence du cisaillement transverse
    phiy = e*xiz*12.0d0*alfay/(xl2*g*a)
    phiz = e*xiy*12.0d0*alfaz/(xl2*g*a)
!
!   recuperation des orientations alpha,beta,gamma  :
    call jevech('PCAORIE', 'L', lorien)
!
!   construction de la matrice de passage pgl du repere global au repere local
    call matrot(zr(lorien), pgl)
!
!   recuperation du champ de deplacement sur l'element
    call jevech('PDEPLAR', 'L', jdepl)
    do i = 1, 14
        depglo(i) = zr(jdepl+i-1)
    end do
!
!   passage des deplacements du repere global au repere local
    call utpvgl(nno, nc, pgl, depglo, depl)
!
!   boucle sur les points de gauss
    do igau = 1, 3
        call r8inir(7, 0.0d0, epsgen, 1)
!       calcul de la matrice (b) reliant les deformations generalisees aux deplacements
!           (DU/DX,GAMAXY,GAMAXZ,D(TETAX)/DX,D(TETAY)/DX,D(TETAZ/DX,D(GRX)/DX)
        call jsd1ff(igau, xl, phiy, phiz, b)
!       calcul des deformations generalisees au point d'integration courant
        do i = 1, 7
            do j = 1, 14
                epsgen(i) = epsgen(i)+b(i, j)*depl(j)
            end do
        end do
!       calcul des efforts generalises au point d'integration courant
        siggen(igau, 1) = e*a*epsgen(1)
        siggen(igau, 2) = alfay*g*a*epsgen(2)
        siggen(igau, 3) = alfaz*g*a*epsgen(3)
        siggen(igau, 4) = xjx*g*epsgen(4)
        siggen(igau, 5) = e*xiy*epsgen(5)
        siggen(igau, 6) = e*xiz*epsgen(6)
        siggen(igau, 7) = e*xjg*epsgen(7)
    end do
!
!   recuperation et affectation du vecteur des efforts generalises en sortie
    if (option .eq. 'SIEF_ELGA') then
        call jevech('PCONTRR', 'E', jeffo)
    else
!       option non programmee
        ASSERT(.false.)
    end if
!
    k = 0
    do igau = 1, 3
        do i = 1, 7
            k = k+1
            zr(jeffo+k-1) = siggen(igau, i)
        end do
    end do
!
end subroutine
