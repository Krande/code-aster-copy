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

subroutine dsqsie(option, fami, xyzl, pgl, depl, &
                  nbcou, cdl)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/dsqbfa.h"
#include "asterfort/dsqbfb.h"
#include "asterfort/dsqcis.h"
#include "asterfort/dsqdis.h"
#include "asterfort/dsqlxy.h"
#include "asterfort/dsxhft.h"
#include "asterfort/dsxhlt.h"
#include "asterfort/dxdmul.h"
#include "asterfort/dxhmft.h"
#include "asterfort/dxmate.h"
#include "asterfort/dxqbm.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/gquad4.h"
#include "asterfort/jevech.h"
#include "asterfort/jquad4.h"
    character(len=4) :: fami
    character(len=16) :: option
    real(kind=8) :: xyzl(3, *), pgl(3, *), depl(*), cdl(*)
    integer(kind=8) :: nbcou
!     RELATION ELAS_COQUE/ELAS_COQMU
!     CONTRAINTES DE L'ELEMENT DE PLAQUE DSQ (SIEF_ELGA)
!     ------------------------------------------------------------------
!     IN  XYZL   : COORDONNEES LOCALES DES QUATRE NOEUDS
!     IN  PGL    : MATRICE DE PASSAGE GLOBAL - LOCAL
!     IN  DEPL   : DEPLACEMENTS
!     OUT CDL    : CONTRAINTES AUX POINTS DE GAUSS DANS LE REPERE LOCAL
!                  LE CALCUL EST FAIT SUR UNE SEULE COUCHE (ELAS_COQUE)
!                  SUR 3 NIVEAUX : 3 PTS D INTEGRATION DANS L EPAISSEUR
!                  CORRESPONDANT AUX NIVEAUX INF, MOY, SUP
    integer(kind=8) :: nnomai
    parameter(nnomai=4)
    integer(kind=8) :: nddlme
    parameter(nddlme=2)
    integer(kind=8) :: nddlfl
    parameter(nddlfl=3)
!
    integer(kind=8) :: ndim, nno, nnos, npg, ipoids, icoopg, ivf, idfdx, idfd2, jgano
    integer(kind=8) :: jcaco, i, j, k, ie, icpg, ig, icou, iniv, multic
    real(kind=8) :: zic, epais, excen
    real(kind=8) :: depf(nddlfl*nnomai), depm(nddlme*nnomai)
    real(kind=8) :: vt(2), lambda(4)
    real(kind=8) :: df(3, 3), dm(3, 3), dmf(3, 3), dc(2, 2), dci(2, 2)
    real(kind=8) :: dmc(3, 2), dfc(3, 2)
    real(kind=8) :: h(3, 3), d1i(2, 2), d2i(2, 4)
    real(kind=8) :: bf(3, nddlfl*nnomai), bm(3, nddlme*nnomai)
    real(kind=8) :: sm(3), sf(3), hft2(2, 6), hlt2(4, 6)
    real(kind=8) :: eps(3), sig(3), cist(2), dcis(2)
    real(kind=8) :: qsi, eta, caraq4(25), t2iu(4), t2ui(4), t1ve(9)
    real(kind=8) :: an(4, 12), bc(2, 12), bcm(2, 8), hmft2(2, 6)
    real(kind=8) :: bfa(3, 4), bfb(3, 12), bfn(3, 12)
    real(kind=8) :: bca(2, 4), bcb(2, 12), bcn(2, 12)
    real(kind=8) :: jacob(5), hicou, zmin, zmax, quotient, a, b, c
    aster_logical :: coupmf, lcalct
!     ------------------------------------------------------------------
!
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jcoopg=icoopg, jvf=ivf, jdfde=idfdx, jdfd2=idfd2, &
                     jgano=jgano)
!
!     ----- RAPPEL DES MATRICES DE RIGIDITE DU MATERIAU EN FLEXION,
!           MEMBRANE ET CISAILLEMENT INVERSEES -------------------------
!     ----- CALCUL DES GRANDEURS GEOMETRIQUES SUR LE QUADRANGLE --------
    call gquad4(xyzl, caraq4)
!
!     ----- CARACTERISTIQUES DES MATERIAUX --------
    call dxmate(fami, df, dm, dmf, dc, &
                dci, dmc, dfc, nno, pgl, &
                multic, coupmf, t2iu, t2ui, t1ve)
!
!     -------- CALCUL DE LA MATRICE DE HOOKE EN MEMBRANE ---------------
    if (multic .eq. 0) then
        call jevech('PCACOQU', 'L', jcaco)
        epais = zr(jcaco)
        hicou = epais/nbcou
        excen = zr(jcaco-1+5)
        h = dm/epais
    end if
!
!     ----- COMPOSANTES DEPLACEMENT MEMBRANE ET FLEXION ----------------
    do j = 1, nnomai
        do i = 1, nddlme
            depm(i+2*(j-1)) = depl(i+6*(j-1))
        end do
        depf(1+3*(j-1)) = depl(1+2+6*(j-1))
        depf(2+3*(j-1)) = depl(3+2+6*(j-1))
        depf(3+3*(j-1)) = -depl(2+2+6*(j-1))
    end do
!     ---- CALCUL DE LA MATRICE AN -------------------------------------
    call dsqdis(xyzl, caraq4, df, dci, an)
!              ---------------------
!
!  BOUCLE SUR LES POINTS D INTEGRATION
!
    if (option .eq. 'EPSI_ELGA') then
        lcalct = .false.
    else
        lcalct = .true.
    end if
!   coefficients pour calcul de sixz et siyz
    if (multic .eq. 0) then
        zmin = excen-epais/2.d0
        zmax = excen+epais/2.d0
        quotient = 1.d0*zmax**3-3*zmax**2*zmin+3*zmax*zmin**2-1.d0*zmin**3
        a = -6.d0/quotient
        b = 6.d0*(zmin+zmax)/quotient
        c = -6.d0*zmax*zmin/quotient
    end if
!
    do ie = 1, npg
        qsi = zr(icoopg-1+ndim*(ie-1)+1)
        eta = zr(icoopg-1+ndim*(ie-1)+2)
!         ----- CALCUL DU JACOBIEN SUR LE QUADRANGLE -----------------
        call jquad4(xyzl, qsi, eta, jacob)
!         ----- CALCUL DE LA MATRICE BM AU POINT QSI ETA -------------
        call dxqbm(qsi, eta, jacob(2), bm)
!         ------ SM = BM.DEPM ----------------------------------------
        do i = 1, 3
            sm(i) = 0.d0
        end do
        do i = 1, 3
            do j = 1, nddlme*nnomai
                sm(i) = sm(i)+bm(i, j)*depm(j)
            end do
        end do
!         ----- CALCUL DE LA MATRICE BFB AU POINT QSI ETA -----------
        call dsqbfb(qsi, eta, jacob(2), bfb)
!         ----- CALCUL DE LA MATRICE BFA AU POINT QSI ETA -----------
        call dsqbfa(qsi, eta, jacob(2), caraq4, bfa)
!         ------ BF = BFB + BFA.AN ----------------------------------
        do i = 1, 3
            do j = 1, 12
                bfn(i, j) = 0.d0
                do k = 1, 4
                    bfn(i, j) = bfn(i, j)+bfa(i, k)*an(k, j)
                end do
                bf(i, j) = bfb(i, j)+bfn(i, j)
            end do
        end do
!         ------ SF = BF.DEPF ---------------------------------------
        do i = 1, 3
            sf(i) = 0.d0
        end do
        do i = 1, 3
            do j = 1, nddlfl*nnomai
                sf(i) = sf(i)+bf(i, j)*depf(j)
            end do
        end do
!
!  BOUCLE SUR LES COUCHES
!
        do icou = 1, nbcou
!
!  BOUCLE SUR LES POINTS D'INTEGRATION DANS L'EPAISSEUR DE LA COUCHE
!
            do ig = 1, 3
!
!           INDICE DANS LE CHAMP DE CONTRAINTES A ECRIRE
                icpg = 6*3*nbcou*(ie-1)+6*3*(icou-1)+6*(ig-1)
!
                if (multic .eq. 0) then
!             -- MONOCOUCHE
!             -- COTE DES POINTS D'INTEGRATION
!             --------------------------------
                    zic = excen-epais/2.d0+(icou-1)*hicou
                    if (ig .eq. 1) then
                        zic = zic
                    else if (ig .eq. 2) then
                        zic = zic+hicou/2.d0
                    else
                        zic = zic+hicou
                    end if
                    d1i(1, 1) = a*zic*zic+b*zic+c
                    d1i(2, 2) = d1i(1, 1)
                    d1i(1, 2) = 0.d0
                    d1i(2, 1) = 0.d0
                else
!             -- EN MULTICOUCHES
!             -- ON CALCULE TOUT D'UN COUP
                    iniv = ig-2
                    call dxdmul(lcalct, icou, iniv, t1ve, t2ui, &
                                h, d1i, d2i, zic, hicou)
                end if
!
                do i = 1, 3
                    eps(i) = sm(i)+zic*sf(i)
                    sig(i) = 0.d0
                end do
!
!
! ---       CALCUL DU PRODUIT HF.T2 :
!           ------------------------
                call dsxhft(df, jacob(2), hft2)
!
! ---       CALCUL DU PRODUIT HMF.T2 :
!           ------------------------
                call dxhmft(dmf, jacob(2), hmft2)
!
! ---       CALCUL DES MATRICES BCB, BCA ET BCM :
!           -----------------------------------
                call dsqcis(qsi, eta, caraq4, hmft2, hft2, &
                            bcm, bcb, bca)
!
!             ------ BC = BCB + BCA.AN ---------------------------------
                do i = 1, 2
                    do j = 1, 12
                        bcn(i, j) = 0.d0
                        do k = 1, 4
                            bcn(i, j) = bcn(i, j)+bca(i, k)*an(k, j)
                        end do
                        bc(i, j) = bcb(i, j)+bcn(i, j)
                    end do
                end do
!             ------ VT = BC.DEPF --------------------------------------
                vt(1) = 0.d0
                vt(2) = 0.d0
                do i = 1, 2
                    do j = 1, 12
                        vt(i) = vt(i)+bc(i, j)*depf(j)
                    end do
                end do
!
                if (option .eq. 'EPSI_ELGA') then
!           ------ DCIS = DCI.VT --------------------------------------
                    dcis(1) = dci(1, 1)*vt(1)+dci(1, 2)*vt(2)
                    dcis(2) = dci(2, 1)*vt(1)+dci(2, 2)*vt(2)
                    cdl(icpg+1) = eps(1)
                    cdl(icpg+2) = eps(2)
                    cdl(icpg+3) = 0.d0
!           --- PASSAGE DE LA DISTORSION A LA DEFORMATION DE CIS. ------
                    cdl(icpg+4) = eps(3)/2.d0
                    cdl(icpg+5) = dcis(1)/2.d0
                    cdl(icpg+6) = dcis(2)/2.d0
                else
!           SIEF_ELGA
                    do i = 1, 3
                        do j = 1, 3
                            sig(i) = sig(i)+h(i, j)*eps(j)
                        end do
                    end do
!
!             ------ CIST = D1I.VT ( + D2I.LAMBDA SI MULTICOUCHES ) ----
                    cist(1) = d1i(1, 1)*vt(1)+d1i(1, 2)*vt(2)
                    cist(2) = d1i(2, 1)*vt(1)+d1i(2, 2)*vt(2)
                    if (multic .gt. 0) then
!             ------- CALCUL DU PRODUIT HL.T2 ------------------------
                        call dsxhlt(df, jacob(2), hlt2)
                        call dsqlxy(qsi, eta, hlt2, an, depf, &
                                    caraq4(13), lambda)
                        do j = 1, 4
                            cist(1) = cist(1)+d2i(1, j)*lambda(j)
                            cist(2) = cist(2)+d2i(2, j)*lambda(j)
                        end do
                    end if
!
                    cdl(icpg+1) = sig(1)
                    cdl(icpg+2) = sig(2)
                    cdl(icpg+3) = 0.d0
                    cdl(icpg+4) = sig(3)
                    cdl(icpg+5) = cist(1)
                    cdl(icpg+6) = cist(2)
                end if
            end do
        end do
    end do
!
end subroutine
