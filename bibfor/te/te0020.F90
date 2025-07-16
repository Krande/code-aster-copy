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

subroutine te0020(nomopt, nomte)
!
!
! --------------------------------------------------------------------------------------------------
!
!                   OPTION DE CALCUL 'CHAR_MECA_EPSI_R'
!
! --------------------------------------------------------------------------------------------------
!
!   IN
!       OPTION  : OPTION DE CALCUL 'CHAR_MECA_EPSI_R'
!       NOMTE   : NOM DU TYPE ELEMENT
!           POUTRES DROITE D'EULER
!               'MECA_POU_D_E'   : SECTION VARIABLE
!               'MECA_POU_D_EM'  : SECTION MULTIFIBRES
!           POUTRE DROITE DE TIMOSHENKO
!               'MECA_POU_D_T'   : SECTION VARIABLE
!               'MECA_POU_D_TG'  : AVEC GAUCHISSEMENT
!               'MECA_POU_D_TGM' : AVEC GAUCHISSEMENT, SECTION MULTIFIBRES
!
! --------------------------------------------------------------------------------------------------
!
    implicit none
    character(len=16) :: nomte, nomopt
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/angvxy.h"
#include "asterfort/fointe.h"
#include "asterfort/jevech.h"
#include "asterfort/matrot.h"
#include "asterfort/normev.h"
#include "asterfort/pmfinfo.h"
#include "asterfort/pmfitg.h"
#include "asterfort/pmfitx.h"
#include "asterfort/poutre_modloc.h"
#include "asterfort/provec.h"
#include "asterfort/rcvalb.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "asterfort/utpvlg.h"
#include "asterc/r8prem.h"
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: lmater, jacf, idefi, ivectu
    integer(kind=8) :: lorien, nno, nc, iabsc, nbpar
    real(kind=8) :: r8bid, e, xnu, g, carsec(6), fs(14)
    real(kind=8) :: a, xiy, xiz, alfay, alfaz, xjx, a2, xiy2, xiz2
    real(kind=8) :: epx, xky, xkz, vect_y(3), norm, vect_x(3)
    real(kind=8) :: pgl(3, 3), angl(3), dgamma, vect_n(3), xkn1, xkn2
!
    integer(kind=8) :: nbres
    parameter(nbres=4)
    integer(kind=8) :: codres(nbres)
    real(kind=8) :: valres(nbres)
    character(len=16) :: nomres(nbres)
!
    integer(kind=8) :: nbfibr, nbgrfi, tygrfi, nbcarm, nug(10), itemps, igeom, ier
    character(len=8) :: nompar(5)
    real(kind=8) :: valpar(5)
!
! --------------------------------------------------------------------------------------------------
    integer(kind=8), parameter :: nb_cara = 9
    real(kind=8) :: vale_cara(nb_cara)
    character(len=8) :: noms_cara(nb_cara)
    data noms_cara/'A1', 'IY1', 'IZ1', 'AY1', 'AZ1', 'JX1', 'A2', 'IY2', 'IZ2'/
! --------------------------------------------------------------------------------------------------
!
    ASSERT(nomopt(1:15) .eq. 'CHAR_MECA_EPSI_')
    nc = 6
!
    call jevech('PMATERC', 'L', lmater)
    call jevech('PVECTUR', 'E', ivectu)
!
    valres(:) = 0.0d+0
    nomres(1) = 'E'
    nomres(2) = 'NU'
    nomres(3) = 'ALPHA'
    nomres(4) = 'RHO'
    call rcvalb('FPG1', 1, 1, '+', zi(lmater), &
                ' ', 'ELAS', 0, ' ', [0.0d+0], &
                4, nomres, valres, codres, 0)
    if (codres(3) .ne. 0) valres(3) = 0.0d+0
    if (codres(4) .ne. 0) valres(4) = 0.0d+0
    e = valres(1)
    xnu = valres(2)
    g = e/(2.0d+0*(1.0d+0+xnu))
!
    call poutre_modloc('CAGNPO', noms_cara, nb_cara, lvaleur=vale_cara)
!   section initiale
    a = vale_cara(1)
    xiy = vale_cara(2)
    xiz = vale_cara(3)
    alfay = vale_cara(4)
    alfaz = vale_cara(5)
    xjx = vale_cara(6)
!   section finale
    a2 = vale_cara(7)
    xiy2 = vale_cara(8)
    xiz2 = vale_cara(9)
!
!   poutres homothétiques interdites car les forces doivent être
!   identiques aux deux noeuds
    if (abs(a-a2) .gt. 1d3*r8prem()) then
        call utmess('F', 'CHARGES_4')
    end if
    if (abs(xiy-xiy2) .gt. 1d3*r8prem()) then
        call utmess('F', 'CHARGES_4')
    end if
    if (abs(xiz-xiz2) .gt. 1d3*r8prem()) then
        call utmess('F', 'CHARGES_4')
    end if

! --------------------------------------------------------------------------------------------------
    if (nomte .eq. 'MECA_POU_D_TGM') then
!       Récupération des caractéristiques des fibres
        call pmfinfo(nbfibr, nbgrfi, tygrfi, nbcarm, nug, jacf=jacf)
        call pmfitg(tygrfi, nbfibr, nbcarm, zr(jacf), carsec)
        a = carsec(1)
        xiy = carsec(5)
        xiz = carsec(4)
    end if
!
    vect_n(:) = 0.d0
    if (nomopt(15:16) .eq. '_R') then
        call jevech('PEPSINR', 'L', idefi)
        epx = zr(idefi)
        xky = zr(idefi+1)
        xkz = zr(idefi+2)
        vect_n(1) = zr(idefi+3)
        vect_n(2) = zr(idefi+4)
        vect_n(3) = zr(idefi+5)
    else
        call jevech('PEPSINF', 'L', idefi)
        call jevech('PINSTR', 'L', itemps)
        call jevech('PGEOMER', 'L', igeom)
!       récupération de VECT_N
        call fointe('FM', zk8(idefi+3), 0, ' ', [0.d0], vect_n(1), ier)
        call fointe('FM', zk8(idefi+4), 0, ' ', [0.d0], vect_n(2), ier)
        call fointe('FM', zk8(idefi+5), 0, ' ', [0.d0], vect_n(3), ier)
!
        nompar(1) = 'X'
        nompar(2) = 'Y'
        nompar(3) = 'Z'
        nompar(4) = 'INST'
        valpar(4) = zr(itemps)
!       recuperation éventuelle de l'abscisse curviligne
        call tecach('ONO', 'PABSCUR', 'L', ier, iad=iabsc)
        nbpar = 4
        if (iabsc .ne. 0) then
            nbpar = nbpar+1
            nompar(nbpar) = 'ABSC'
            valpar(nbpar) = (zr(iabsc)+zr(iabsc+1))/2.d0
        end if
!       milieu de la poutre
        valpar(1) = (zr(igeom+(2-1)*3-1+1)+zr(igeom-1+1))/2.d0
        valpar(2) = (zr(igeom+(2-1)*3-1+2)+zr(igeom-1+2))/2.d0
        valpar(3) = (zr(igeom+(2-1)*3-1+3)+zr(igeom-1+3))/2.d0
        call fointe('FM', zk8(idefi), nbpar, nompar, valpar, &
                    epx, ier)
        call fointe('FM', zk8(idefi+1), nbpar, nompar, valpar, &
                    xky, ier)
        call fointe('FM', zk8(idefi+2), nbpar, nompar, valpar, &
                    xkz, ier)
    end if
!   Récupération des orientations alpha, beta, gamma
    call jevech('PCAORIE', 'L', lorien)
!   Matrice de rotation mgl
    call matrot(zr(lorien), pgl)
!   Traitement du repère imposé
    call normev(vect_n, norm)
    if (norm .gt. r8prem()) then
        vect_x(:) = pgl(1, :)
        call provec(vect_n, vect_x, vect_y)
        call normev(vect_y, norm)
        call angvxy(vect_x, vect_y, angl)
        dgamma = angl(3)-zr(lorien-1+3)
        xkn1 = xky
        xkn2 = xkz
        xky = cos(dgamma)*xkn1-sin(dgamma)*xkn2
        xkz = sin(dgamma)*xkn1+cos(dgamma)*xkn2
    end if
!
    fs(1) = e*a*epx
    fs(2) = 0.d0
    fs(3) = 0.d0
    fs(4) = 0.d0
    fs(5) = e*xiy*xky
    fs(6) = e*xiz*xkz
    if ((nomte .eq. 'MECA_POU_D_TG')) then
        fs(7) = 0.d0
        fs(8) = e*a*epx
        fs(9) = 0.d0
        fs(10) = 0.d0
        fs(11) = 0.d0
        fs(12) = e*xiy*xky
        fs(13) = e*xiz*xkz
        fs(14) = 0.d0
!
        nc = 7
    else if (nomte .eq. 'MECA_POU_D_TGM') then
!       Récupération des caractéristiques des fibres
        call pmfitx(zi(lmater), 1, carsec, r8bid)
        fs(1) = carsec(1)*epx
        fs(5) = carsec(5)*xky
        fs(6) = carsec(4)*xkz
        fs(7) = 0.d0
        fs(8) = carsec(1)*epx
        fs(9) = 0.d0
        fs(10) = 0.d0
        fs(11) = 0.d0
        fs(12) = carsec(5)*xky
        fs(13) = carsec(4)*xkz
        fs(14) = 0.d0
!
        nc = 7
    else if (nomte .eq. 'MECA_POU_D_EM') then
!       Récupération des caractéristiques des fibres
        call pmfitx(zi(lmater), 1, carsec, r8bid)
        fs(1) = carsec(1)*epx
        fs(5) = carsec(5)*xky
        fs(6) = carsec(4)*xkz
        fs(7) = carsec(1)*epx
        fs(8) = 0.d0
        fs(9) = 0.d0
        fs(10) = 0.d0
        fs(11) = carsec(5)*xky
        fs(12) = carsec(4)*xkz
    else
        fs(7) = e*a2*epx
        fs(8) = 0.d0
        fs(9) = 0.d0
        fs(10) = 0.d0
        fs(11) = e*xiy2*xky
        fs(12) = e*xiz2*xkz
    end if
    fs(1) = -fs(1)
    fs(2) = -fs(2)
    fs(3) = -fs(3)
    fs(4) = -fs(4)
    fs(5) = -fs(5)
    fs(6) = -fs(6)
!   Passage en repère global
    nno = 2
    call utpvlg(nno, nc, pgl, fs, zr(ivectu))
end subroutine
