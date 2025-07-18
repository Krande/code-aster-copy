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
subroutine dktrge(nomte, xyzl, pgl, rig)
    implicit none
#include "jeveux.h"
#include "asterc/r8dgrd.h"
#include "asterfort/assert.h"
#include "asterfort/coqrep.h"
#include "asterfort/cosiro.h"
#include "asterfort/dktbnl.h"
#include "asterfort/dxefro.h"
#include "asterfort/dxtloc.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/gtria3.h"
#include "asterfort/jevech.h"
#include "asterfort/prmama.h"
#include "asterfort/r8inir.h"
#include "asterfort/tecach.h"
!
    real(kind=8) :: xyzl(3, *), pgl(*), rig(*)
    character(len=16) :: nomte
!
! ======================================================================
!
!     matrice de rigidite geometrique de l'element de plaque dkt
!        option : rigi_meca_ge
!     ------------------------------------------------------------------
!     in  nomte  : nom du type element
!     in  xyzl   : coordonnees locales des trois noeuds
!     in  option : option rigi_meca_ge
!     in  pgl    : matrice de passage global/local
!     out rig    : matrice de rigidite geometrique
!     ------------------------------------------------------------------
!
    integer(kind=8) :: nbsig
    parameter(nbsig=6)
    integer(kind=8) :: nbcon
    parameter(nbcon=8)
!
    integer(kind=8) :: ndim, nno, nnoel, npg, ipoids, icoopg, ivf, idfdx, idfd2, jgano
    integer(kind=8) :: ipg, iret
    integer(kind=8) :: jsigm, nbsp, npgh, icou, idec, jtab(7)
    integer(kind=8) :: nbcou, jcoqu, j2, ier
    integer(kind=8) :: jgeom, i, j
!
    real(kind=8) :: poids, cara(25)
    real(kind=8) :: bnl(2, 9), bnli(9, 2)
    real(kind=8) :: flex(9, 9), memb(36), mefl(54)
    real(kind=8) :: flexi(9, 9)
    real(kind=8) :: ctor
    real(kind=8) :: effint(24), effgt(24), alpha, beta
    real(kind=8) :: t2ev(4), t2ve(4), c, s
    real(kind=8) :: epi, h, hb, hm, hh, cb, cm, ch
!
! --- contraintes
    real(kind=8) :: sixxb, siyyb, sixyb
    real(kind=8) :: sixxm, siyym, sixym
    real(kind=8) :: sixxh, siyyh, sixyh
    real(kind=8) :: nxx, nyy, nxy, normal(2, 2)
!
! deb ------------------------------------------------------------------
!
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnoel, npg=npg, &
                     jpoids=ipoids, jcoopg=icoopg, jvf=ivf, jdfde=idfdx, jdfd2=idfd2, &
                     jgano=jgano)
!
!     ------ mise a zero des matrices : flex et mefl -------------------
    call r8inir(81, 0.d0, flex, 1)
    call r8inir(36, 0.d0, memb, 1)
    call r8inir(54, 0.d0, mefl, 1)
!
    call jevech('PGEOMER', 'L', jgeom)
!
!     -- epaisseur :
!
    call jevech('PCACOQU', 'L', jcoqu)
    h = zr(jcoqu)
    alpha = zr(jcoqu+1)*r8dgrd()
    beta = zr(jcoqu+2)*r8dgrd()
    ctor = zr(jcoqu+3)
!
!     -- nombre de couches :
!
    if (nomte .eq. 'MEDKTR3') then
        call jevech('PNBSP_I', 'L', j2)
        nbcou = zi(j2)
    end if
!
!     -- contraintes dans les couches :
!     ----------------------------------
    call tecach('OOO', 'PCONTRR', 'L', iret, nval=7, &
                itab=jtab)
    jsigm = jtab(1)
    npg = jtab(3)
    nbsp = jtab(7)
    npgh = 3
!
    if (nomte .eq. 'MEDKTR3') then
        ASSERT(nbsp .eq. nbcou*npgh)
        ASSERT(jtab(2) .eq. nbsig*npg)
    end if
!
! ---     passage des contraintes dans le repere intrinseque :
!
    if (nomte .eq. 'MEDKTR3') then
        call cosiro(nomte, 'PCONTRR', 'L', 'UI', 'G', &
                    jsigm, 'S')
    else if (nomte .eq. 'MEDKTG3') then
        call tecach('OOO', 'PCONTRR', 'L', iret, nval=7, &
                    itab=jtab)
        jsigm = jtab(1)
        do i = 1, nbcon*npg
            effgt(i) = zr(jsigm-1+i)
        end do
        call coqrep(pgl, alpha, beta, t2ev, t2ve, &
                    c, s)
        call dxefro(npg, t2ev, effgt, effint)
    end if
!
!     ----- calcul des grandeurs geometriques sur le quadrangle --------
!
    call gtria3(xyzl, cara)
!
! --- calcul de la matrice bnl deformations non-lineaires de membrane
!
    call dktbnl(cara(9), bnl)
!
! ---- boucle sur les points d'integration :
!      ===================================
    do ipg = 1, npg
!
        nxx = 0.d0
        nyy = 0.d0
        nxy = 0.d0
!
        poids = zr(ipoids+ipg-1)*cara(7)
!
        if (nomte .eq. 'MEDKTR3') then
!       -- boucle sur les couches :
            hb = -h/2
            do icou = 1, nbcou
                idec = ((ipg-1)*nbcou+(icou-1))*npgh*nbsig
                epi = h/nbcou
                hm = hb+epi/2.d0
                hh = hm+epi/2.d0
!         -- sixxb, siyyb, ... : contraintes au bas de la couche
                sixxb = zr(jsigm-1+idec+1)
                siyyb = zr(jsigm-1+idec+2)
                sixyb = zr(jsigm-1+idec+4)
!         -- sixxm, siyym, ... : contraintes au milieu de la couche
                sixxm = zr(jsigm-1+idec+1+nbsig)
                siyym = zr(jsigm-1+idec+2+nbsig)
                sixym = zr(jsigm-1+idec+4+nbsig)
!         -- sixxh, siyyh, ... : contraintes en haut de la couche
                sixxh = zr(jsigm-1+idec+1+2*nbsig)
                siyyh = zr(jsigm-1+idec+2+2*nbsig)
                sixyh = zr(jsigm-1+idec+4+2*nbsig)
!         -- on integre dans l'epaisseur de chaque couche
!            avec une forrmule de newton-cotes a 3 points
!            les coefficients sont 1/6, 4/6 et 1/6
                cb = epi/6
                cm = 4.d0*epi/6
                ch = epi/6
!         -- nxx, nyy, nxy = somme de sixx, siyy, sixy :
                nxx = nxx+cb*sixxb+cm*sixxm+ch*sixxh
                nyy = nyy+cb*siyyb+cm*siyym+ch*siyyh
                nxy = nxy+cb*sixyb+cm*sixym+ch*sixyh
!         -- mise a jour de hb pour la couche suivante :
                hb = hb+epi
!
! --- fin de la boucle sur les couches
!
            end do
        else if (nomte .eq. 'MEDKTG3') then
            nxx = effint((ipg-1)*nbcon+1)
            nyy = effint((ipg-1)*nbcon+2)
            nxy = effint((ipg-1)*nbcon+3)
        end if
!
        normal(1, 1) = nxx*poids
        normal(2, 2) = nyy*poids
        normal(1, 2) = nxy*poids
        normal(2, 1) = normal(1, 2)
!
        call prmama(3, bnl, 2, 2, 9, &
                    normal, 2, 2, 2, bnli, &
                    9, 9, 2, ier)
        flexi = matmul(bnli, bnl)
!
        do i = 1, 9
            do j = 1, 9
                flex(i, j) = flex(i, j)+flexi(i, j)
            end do
        end do
!
! --- fin de la boucle sur le points d'integration
!
    end do
!
! --- stockage
!
    call dxtloc(flex, memb, mefl, ctor, rig)
!
end
