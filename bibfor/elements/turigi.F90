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
subroutine turigi(nomte, nbrddl, k)
! aslint: disable=W1306
    implicit none
#include "jeveux.h"
#include "asterc/r8pi.h"
#include "asterfort/bcoudc.h"
#include "asterfort/bcoude.h"
#include "asterfort/carcou.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/kcoude.h"
#include "asterfort/klg.h"
#include "asterfort/klgcou.h"
#include "asterfort/mavec.h"
#include "asterfort/moytem.h"
#include "asterfort/poutre_modloc.h"
#include "asterfort/rcvala.h"
#include "asterfort/utmess.h"
    character(len=16) :: nomte
    integer(kind=8) :: nbrddl, nc
    integer(kind=8) :: ndim, nnos, idfdk, jdfd2, jgano
    real(kind=8) :: b(4, nbrddl), k(nbrddl, nbrddl)
! ......................................................................
!
!    - FONCTION REALISEE:  CALCUL DES MATRICES ELEMENTAIRES
!                          TUYAU
!                          OPTION : RIGI_MECA
!    - ARGUMENTS:
!        DONNEES:      B, MASS1,MASS,K : MATRICES
!                      DIMENSIONNEES PAR LE TE0582 APPELANT
!
! ......................................................................
!
!
!
!     VARIABLES LOCALES
!
    integer(kind=8) :: nbres, icoude, jnbspi, nbsecm, nbcoum, nspg
    parameter(nbres=9)
    character(len=8) :: nompar
    character(len=16) :: nomres(nbres)
    integer(kind=8) :: icodre(nbres)
    parameter(nbsecm=32, nbcoum=10)
    real(kind=8) :: poicou(2*nbcoum+1), poisec(2*nbsecm+1)
    real(kind=8) :: valres(nbres), valpar, theta
    real(kind=8) :: e, nu, h, a, l, r1
    real(kind=8) :: pi, deuxpi, fi, sinfi
    real(kind=8) :: beta, cisail, g, poids, r, rayon
    real(kind=8) :: c(4, 4), xpg(4)
    real(kind=8) :: pgl(3, 3), pgl4(3, 3)
    real(kind=8) :: pgl1(3, 3), pgl2(3, 3), pgl3(3, 3), omega
    integer(kind=8) :: nno, npg, nbcou, nbsec, m
    integer(kind=8) :: ipoids, ivf, iret
    integer(kind=8) :: imate, imatuu, igeom, nbpar
    integer(kind=8) :: igau, icou, isect, i, j, lorien, icoud2, mmt
    integer(kind=8) :: jcoopg
!
    integer(kind=8), parameter :: nb_cara1 = 2
    real(kind=8) :: vale_cara1(nb_cara1)
    character(len=8) :: noms_cara1(nb_cara1)
    data noms_cara1/'R1', 'EP1'/
!-----------------------------------------------------------------------
!
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jcoopg=jcoopg, jvf=ivf, jdfde=idfdk, jdfd2=jdfd2, &
                     jgano=jgano)
!
!     DIMENSION DE LA MATRICE STOCKEE SOUS FORME VECTEUR
    nc = nbrddl*(nbrddl+1)/2
!
    pi = r8pi()
    deuxpi = 2.d0*pi
!
    call jevech('PNBSP_I', 'L', jnbspi)
    nbcou = zi(jnbspi-1+1)
    nbsec = zi(jnbspi-1+2)
!
!     -- CALCUL DES POIDS DES COUCHES ET DES SECTEURS:
    poicou(1) = 1.d0/3.d0
    do i = 1, nbcou-1
        poicou(2*i) = 4.d0/3.d0
        poicou(2*i+1) = 2.d0/3.d0
    end do
    poicou(2*nbcou) = 4.d0/3.d0
    poicou(2*nbcou+1) = 1.d0/3.d0
    poisec(1) = 1.d0/3.d0
    do i = 1, nbsec-1
        poisec(2*i) = 4.d0/3.d0
        poisec(2*i+1) = 2.d0/3.d0
    end do
    poisec(2*nbsec) = 4.d0/3.d0
    poisec(2*nbsec+1) = 1.d0/3.d0
!
!
    call jevech('PGEOMER', 'L', igeom)
    call poutre_modloc('CAGEP1', noms_cara1, nb_cara1, lvaleur=vale_cara1)
    r1 = vale_cara1(1)
    h = vale_cara1(2)
    a = r1-h/2.d0
! A= RMOY, H = EPAISSEUR
!
!
    m = 3
    if (nomte .eq. 'MET6SEG3') m = 6
!
!
!
    do i = 1, npg
        xpg(i) = zr(jcoopg-1+i)
    end do
!     --- RECUPERATION DES ORIENTATIONS ---
    call jevech('PCAORIE', 'L', lorien)
    call carcou(zr(lorien), l, pgl, rayon, theta, &
                pgl1, pgl2, pgl3, pgl4, nno, &
                omega, icoud2)
    if (icoud2 .ge. 10) then
        icoude = icoud2-10
        mmt = 0
        if (h/a .gt. (0.25d0)) then
            call utmess('A', 'ELEMENTS4_54')
        end if
    else
        icoude = icoud2
        mmt = 1
    end if
    call jevech('PMATERC', 'L', imate)
    nomres(1) = 'E'
    nomres(2) = 'NU'
!
!       -- CALCUL DES TEMPERATURES INF, SUP ET MOY
!          (MOYENNE DES NNO NOEUDS) ET DES COEF. DES POLY. DE DEGRE 2 :
!          ------------------------------------------------------------
    nspg = (2*nbsec+1)*(2*nbcou+1)
    iret = 0
    call moytem('RIGI', npg, nspg, '+', valpar, &
                iret)
    if (iret .ne. 0) valpar = 0.d0
    nbpar = 1
    nompar = 'TEMP'
    call rcvala(zi(imate), ' ', 'ELAS', nbpar, nompar, &
                [valpar], 2, nomres, valres, icodre, &
                1)
    e = valres(1)
    nu = valres(2)
! DEFINITION DE LA MATRICE DE COMPORTEMENT C
!
    beta = 1.d0/(1.d0-nu**2)
    g = 1.d0/(2.d0*(1.d0+nu))
    cisail = 1.d0
!
! ON MULTIPLIERA PAR E PLUS LOIN
!
    c(1, 1) = beta
    c(1, 2) = nu*beta
    c(1, 3) = 0.d0
    c(1, 4) = 0.d0
!
    c(2, 1) = nu*beta
    c(2, 2) = beta
    c(2, 3) = 0.d0
    c(2, 4) = 0.d0
!
    c(3, 1) = 0.d0
    c(3, 2) = 0.d0
    c(3, 3) = g
    c(3, 4) = 0.d0
!
    c(4, 1) = 0.d0
    c(4, 2) = 0.d0
    c(4, 3) = 0.d0
    c(4, 4) = g*cisail
!
!
!     FIN DE LA CONSTRUCTION DE LA MATRICE DE COMPORTEMENT C
!
!  INITIALISATION DE LA MATRICE K
!
    k(:, :) = 0.d0
!
! BOUCLE SUR LES POINTS DE GAUSS
!
    do igau = 1, npg
!
! BOUCLE SUR LES POINTS DE SIMPSON DANS L'EPAISSEUR
!
        do icou = 1, 2*nbcou+1
            if (mmt .eq. 0) then
                r = a
            else
                r = a+(icou-1)*h/(2.d0*nbcou)-h/2.d0
            end if
!
! BOUCLE SUR LES POINTS DE SIMPSON SUR LA CIRCONFERENCE
!
            do isect = 1, 2*nbsec+1
                if (icoude .eq. 0) then
                    call bcoude(igau, icou, isect, l, h, &
                                a, m, nno, nbcou, nbsec, &
                                zr(ivf), zr(idfdk), zr(jdfd2), mmt, b)
                    poids = zr(ipoids-1+igau)*poicou(icou)*poisec(isect)*(l/2.d0)*h*deuxpi/(4.&
                            &d0*nbcou*nbsec)*r
!
                else if (icoude .eq. 1) then
                    fi = (isect-1)*deuxpi/(2.d0*nbsec)
!               FI = FI - OMEGA
                    sinfi = sin(fi)
                    l = theta*(rayon+r*sinfi)
                    call bcoudc(igau, icou, isect, h, a, &
                                m, omega, xpg, nno, nbcou, &
                                nbsec, zr(ivf), zr(idfdk), zr(jdfd2), rayon, &
                                theta, mmt, b)
                    poids = zr(ipoids-1+igau)*poicou(icou)*poisec(isect)*(l/2.d0)*h*deuxpi/(4.&
                            &d0*nbcou*nbsec)*r
!
                end if
!  LE DERNIER TERME C'EST R DU  R DFI DX DZETA
!
                call kcoude(nbrddl, poids, b, c, k)
            end do
        end do
    end do
!
! FIN DU CALCUL DE LA MATRICE DE RIGIDITE
!
! MULTIPLICATION PAR LE MODULE DE YOUNG E
!
    do i = 1, nbrddl
        do j = 1, nbrddl
            k(i, j) = e*k(i, j)
        end do
    end do
!
!  CHANGEMENT DE REPERE :
!  PASSAGE DU REPERE LOCAL AU REPERE GLOBAL
!
    if (icoude .eq. 0) then
        call klg(nno, nbrddl, pgl, k)
    else if (icoude .eq. 1) then
        call klgcou(nno, nbrddl, pgl1, pgl2, pgl3, &
                    pgl4, k)
    end if
!
! STOCKAGE DE LA MATRICE DE RIGIDITE
!
    call jevech('PMATUUR', 'E', imatuu)
    call mavec(k, nbrddl, zr(imatuu), nc)
!
end subroutine
