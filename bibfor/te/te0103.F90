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
subroutine te0103(option, nomte)
    implicit none
#include "jeveux.h"
#include "asterfort/cq3d2d.h"
#include "asterfort/dfdm1d.h"
#include "asterfort/dfdm2d.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/fointe.h"
#include "asterfort/jevech.h"
#include "asterfort/jevete.h"
#include "asterfort/r8inir.h"
!
    character(len=16) :: option, nomte
! ......................................................................
!    - FONCTION REALISEE:  CALCUL DES MATRICES ELEMENTAIRES
!                          OPTION : 'RIGI_THER_ECHA_F'
!                          CAS COQUE
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
! ......................................................................
    integer(kind=8) :: ndimax, nbvar, kp, kq, nbddl, mzr
    parameter(ndimax=27)
    parameter(nbvar=4)
    character(len=8) :: nompar(nbvar)
    real(kind=8) :: b(3, 3), theta, zero, rigith(ndimax, ndimax)
    real(kind=8) :: coor2d(18), valpar(nbvar), cour, cosa, sina
    real(kind=8) :: dfdx(9), dfdy(9), poids, x, y, z, xgau, ygau, zgau, pk, un
    real(kind=8) :: long, matp(3, 3), matn(3, 3), h, hmoins, hplus, poi1, poi2
    integer(kind=8) :: nno, nnos, npg2, gi, pi, gj, pj, k, imattt, ndim, jgano
    integer(kind=8) :: idfde, igeom, itemps, ipoids, ivf, i, j, icoefh, ier, nbv, ind
!
!
    if (nomte .ne. 'THCPSE3' .and. nomte .ne. 'THCASE3' .and. nomte .ne. 'THCOSE3' .and. &
        nomte .ne. 'THCOSE2 ') then
        call elrefe_info(fami='MASS', ndim=ndim, nno=nno, nnos=nnos, npg=npg2, &
                         jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
    else
        call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg2, &
                         jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
    end if
!
! --- INITIALISATIONS :
!     ---------------
    zero = 0.0d0
    un = 1.0d0
!
    do i = 1, ndimax
        do j = 1, ndimax
            rigith(i, j) = zero
        end do
    end do
!
! --- POUR LE CALCUL DU COEFFICIENT D'ECHANGE :
!     ---------------------------------------
    nompar(1) = 'X'
    nompar(2) = 'Y'
    nompar(3) = 'Z'
    nompar(4) = 'INST'
!
!
! --- RECUPERATION DES COORDONNEES DES NOEUDS DE L'ELEMENT :
!     ----------------------------------------------------
    call jevech('PGEOMER', 'L', igeom)
!
! --- RECUPERATION DU COEFFICIENT D'ECHANGE (QUI EST UNE FONCTION) :
!     ------------------------------------------------------------
    call jevech('PCOEFHF', 'L', icoefh)
!
! --- RECUPERATION DE L'INSTANT DU CALCUL ET
! --- DU PARAMETRE THETA DE LA METHODE 'THETA' UTILISEE
! --- POUR RESOUDRE L'EQUATION DIFFERENTIELLE EN TEMPS DE LA
! --- TEMPERATURE (EN STATIONNAIRE THETA =1 ) :
!     ---------------------------------------
    call jevech('PINSTR', 'L', itemps)
!
! --- ON RAPPELLE QUE NOMPAR(4) = 'INST' :
!     ----------------------------------
    valpar(4) = zr(itemps)
    theta = zr(itemps+2)
!
! --- DETERMINATION DE LA CONTRIBUTION A LA RIGIDITE THERMIQUE
! --- DES ECHANGES DE LA COQUE AVEC L'EXTERIEUR AU NIVEAU DES
! --- FEUILLETS INFERIEUR ET SUPERIEUR :
!     ================================
!
! --- CAS DES COQUES SURFACIQUES :
!     --------------------------
    if (nomte .ne. 'THCPSE3' .and. nomte .ne. 'THCASE3' .and. nomte .ne. 'THCOSE3' .and. &
        nomte .ne. 'THCOSE2') then
!
! --- DETERMINATION DES COORDONNEES COOR2D DES NOEUDS DE L'ELEMENT
! --- DANS LE REPERE DE L'ELEMENT :
!     ---------------------------
        call cq3d2d(nno, zr(igeom), un, zero, coor2d)
!
! --- BOUCLE SUR LES POINTS D'INTEGRATION :
!     -----------------------------------
        do kp = 1, npg2
            k = (kp-1)*nno
            call dfdm2d(nno, kp, ipoids, idfde, coor2d, &
                        poids, dfdx, dfdy)
!
            xgau = zero
            ygau = zero
            zgau = zero
!
! ---   COORDONNEES DES POINTS D'INTEGRATION :
!       ------------------------------------
            do i = 1, nno
                xgau = xgau+zr(igeom+3*(i-1))*zr(ivf+k+i-1)
                ygau = ygau+zr(igeom+3*(i-1)+1)*zr(ivf+k+i-1)
                zgau = zgau+zr(igeom+3*(i-1)+2)*zr(ivf+k+i-1)
            end do
!
            valpar(1) = xgau
            valpar(2) = ygau
            valpar(3) = zgau
!
! ---   CALCUL DU COEFFICIENT D'ECHANGE DU FEUILLET INFERIEUR DE LA
! ---   COQUE AVEC L'EXTERIEUR AU POINT D'INTEGRATION COURANT
! ---   ET A L'INSTANT COURANT :
!       ----------------------
            call fointe('FM', zk8(icoefh), nbvar, nompar, valpar, &
                        hmoins, ier)
!
! ---   CALCUL DU COEFFICIENT D'ECHANGE DU FEUILLET SUPERIEUR DE LA
! ---   COQUE AVEC L'EXTERIEUR AU POINT D'INTEGRATION COURANT
! ---   ET A L'INSTANT COURANT :
!       ----------------------
            call fointe('FM', zk8(icoefh+1), nbvar, nompar, valpar, &
                        hplus, ier)
!
! ---   CONTRIBUTION AU TENSEUR DE CONDUCTIVITE TRANSVERSE B DES
! ---   ECHANGES AVEC L'EXTERIEUR
! ---            (0 0  0 )
! ---       B =  (0 H- 0 )
! ---            (0 0  H+)
!       -------------------
            b(1, 1) = zero
            b(2, 1) = zero
            b(2, 2) = hmoins
            b(3, 1) = zero
            b(3, 2) = zero
            b(3, 3) = hplus
!
! ---   CALCUL DE LA RIGIDITE THERMIQUE DUE AU TERME D'ECHANGE B :
!       --------------------------------------------------------
            do gi = 1, nno
                do gj = 1, gi
                    do pi = 1, 3
                        do pj = 1, pi
                            pk = b(pi, pj)*zr(ivf+k+gi-1)*zr(ivf+k+gj-1)*poids*theta
!
! ---     AFFECTATION DES TERMES HORS DIAGONAUX DE LA TRIANGULAIRE
! ---     INFERIEURE DE LA SOUS-MATRICE :
!         -----------------------------
                            if ((pi .ne. pj) .and. (gi .ne. gj)) then
                                i = 3*(gi-1)+pj
                                j = 3*(gj-1)+pi
                                rigith(i, j) = rigith(i, j)+pk
                            end if
!
! ---     AFFECTATION DES TERMES DE LA TRIANGULAIRE SUPERIEURE
! ---     DE LA SOUS-MATRICE :
!         ------------------
                            i = 3*(gi-1)+pi
                            j = 3*(gj-1)+pj
                            rigith(i, j) = rigith(i, j)+pk
                        end do
                    end do
                end do
            end do
        end do
!
! --- CAS DES COQUES LINEIQUES (EN CONTRAINTES PLANES ET AXI) :
!     -------------------------------------------------------
    else if (nomte .eq. 'THCPSE3' .or. nomte .eq. 'THCASE3') &
        then
!
! --- BOUCLE SUR LES POINTS D'INTEGRATION :
!     -----------------------------------
        do kp = 1, npg2
            k = (kp-1)*nno
            call dfdm1d(nno, zr(ipoids+kp-1), zr(idfde+k), zr(igeom), dfdx, &
                        cour, poids, cosa, sina)
            xgau = zero
            ygau = zero
!
! ---   COORDONNEES DES POINTS D'INTEGRATION :
!       ------------------------------------
            do i = 1, nno
                xgau = xgau+zr(igeom+2*(i-1))*zr(ivf+k+i-1)
                ygau = ygau+zr(igeom+2*(i-1)+1)*zr(ivf+k+i-1)
            end do
!
            if (nomte .eq. 'THCASE3') poids = poids*xgau
!
            valpar(1) = xgau
            valpar(2) = ygau
!
            nbv = 3
!
! ---   CALCUL DU COEFFICIENT D'ECHANGE AVEC L'EXTERIEUR DU FEUILLET
! ---   INFERIEUR DE LA COQUE AU POINT D'INTEGRATION COURANT
! ---   ET A L'INSTANT COURANT :
!       ----------------------
            call fointe('FM', zk8(icoefh), nbv, nompar, valpar, &
                        hmoins, ier)
!
! ---   CALCUL DU COEFFICIENT D'ECHANGE AVEC L'EXTERIEUR DU FEUILLET
! ---   SUPERIEUR DE LA COQUE AU POINT D'INTEGRATION COURANT
! ---   ET A L'INSTANT COURANT :
!       ----------------------
            call fointe('FM', zk8(icoefh+1), nbv, nompar, valpar, &
                        hplus, ier)
!
! ---   CONTRIBUTION AU TENSEUR DE CONDUCTIVITE TRANSVERSE B DES
! ---   ECHANGES AVEC L'EXTERIEUR
! ---            (0 0  0 )
! ---       B =  (0 H- 0 )
! ---            (0 0  H+)
!       -------------------
            b(1, 1) = zero
            b(2, 1) = zero
            b(2, 2) = hmoins
            b(3, 1) = zero
            b(3, 2) = zero
            b(3, 3) = hplus
!
! ---   CALCUL DE LA RIGIDITE THERMIQUE DUE AU TERME D'ECHANGE B :
!       --------------------------------------------------------
            do gi = 1, nno
                do gj = 1, gi
                    do pi = 1, 3
                        do pj = 1, pi
                            pk = b(pi, pj)*zr(ivf+k+gi-1)*zr(ivf+k+gj-1)*poids*theta
!
! ---     AFFECTATION DES TERMES HORS DIAGONAUX DE LA TRIANGULAIRE
! ---     INFERIEURE DE LA SOUS-MATRICE :
!         -----------------------------
                            if ((pi .ne. pj) .and. (gi .ne. gj)) then
                                i = 3*(gi-1)+pj
                                j = 3*(gj-1)+pi
                                rigith(i, j) = rigith(i, j)+pk
                            end if
!
! ---     AFFECTATION DES TERMES DE LA TRIANGULAIRE SUPERIEURE
! ---     DE LA SOUS-MATRICE :
!         ------------------
                            i = 3*(gi-1)+pi
                            j = 3*(gj-1)+pj
                            rigith(i, j) = rigith(i, j)+pk
                        end do
                    end do
                end do
            end do
        end do
!
! --- CAS DES COQUES LINEIQUES (AUTRES QUE CONTRAINTES PLANES ET AXI) :
!     --------------------------------------------------------------
    else if (nomte .eq. 'THCOSE3' .or. nomte .eq. 'THCOSE2') &
        then
!
        call jevete('&INEL.'//nomte(1:8)//'.DEMR', ' ', mzr)
!CC     CALL JEVECH('PCACOQU','L',ICACOQ)
!
!CC     EP=ZR(ICACOQ)
!
        long = ( &
               zr(igeom+3)-zr(igeom))**2+(zr(igeom+3+1)-zr(igeom+1))**2+(zr(igeom+3+2)-zr(ig&
               &eom+2) &
               )**2
        long = sqrt(long)/2.d0
!       EP  =EP/2.D0
!
! ---   DETERMINATION DE LA 'PART' DE RIGIDITE THERMIQUE DU A L'ECHANGE
! ---   AVEC L'EXTERIEUR POUR LES COQUES LINEIQUES
! ---   ON RAPPELLE QUE LE TERME GENERIQUE DE CETTE MATRICE A POUR
! ---   EXPRESSION :
! ---   B(I,J) = SOMME_VOLUME(H*NI(X,Y,Z)*NJ(X,Y,Z).DX.DY.DZ)
! ---   SOIT
! ---   B(I,J) = SOMME_LONGUEUR (H*NI(X,Y)*NJ(X,Y).DX.DY)
! ---           *SOMME_EPAISSEUR(PK(Z)*PL(Z).DZ)
! ---   OU LES PK ET PL SONT LES FONCTIONS D'INTERPOLATION DANS
! ---   L'EPAISSEUR
! ---   PLUS EXACTEMENT P1(Z), P2(Z), P3(Z) SONT LES POLYNOMES
! ---   DE LAGRANGE (OU AUTRES) DE DEGRE 2 RELATIFS A L'INTERPOLATION
! ---   DE LA TEMPERATURE DANS L'EPAISSEUR TELS QUE
! ---   P1 EST RELATIF A LA TEMPERATURE MOYENNE
! ---   P2 EST RELATIF A LA TEMPERATURE INFERIEURE
! ---   P3 EST RELATIF A LA TEMPERATURE SUPERIEURE
! ---   (I.E. T(X,Y,Z) =    P1(Z)*TMOY(X,Y)
! ---                     + P2(Z)*TINF(X,Y)
! ---                     + P3(Z)*TSUP(X,Y)) :
!       ------------------------------------
        call r8inir(9, zero, matp, 1)
        call r8inir(9, zero, matn, 1)
!
! ---   DETERMINATION DE LA MATRICE MATP DONT LE TERME GENERIQUE
! ---   EST MATP(I,J) = SOMME_EPAISSEUR(PI(Z)*PJ(Z).DZ) :
!       -----------------------------------------------
        do kp = 1, npg2
            kq = (kp-1)*3
!
            poi1 = zr(mzr-1+12+kp)
!
            matp(1, 1) = matp(1, 1)+poi1*zr(mzr-1+kq+1)**2
            matp(1, 2) = matp(1, 2)+poi1*zr(mzr-1+kq+1)*zr(mzr-1+kq+2)
            matp(1, 3) = matp(1, 3)+poi1*zr(mzr-1+kq+1)*zr(mzr-1+kq+3)
            matp(2, 1) = matp(1, 2)
            matp(2, 2) = matp(2, 2)+poi1*zr(mzr-1+kq+2)**2
            matp(2, 3) = matp(2, 3)+poi1*zr(mzr-1+kq+2)*zr(mzr-1+kq+3)
            matp(3, 1) = matp(1, 3)
            matp(3, 2) = matp(2, 3)
            matp(3, 3) = matp(3, 3)+poi1*zr(mzr-1+kq+3)**2
        end do
!
        nompar(1) = 'X'
        nompar(2) = 'Y'
        nompar(3) = 'Z'
        nompar(4) = 'INST'
!
        do kp = 1, npg2
            k = (kp-1)*nno
!
            poi2 = zr(ipoids-1+kp)
!
            x = zero
            y = zero
            z = zero
            do i = 1, nno
                x = x+zr(igeom+3*(i-1))*zr(ivf+k+i-1)
                y = y+zr(igeom+3*(i-1)+1)*zr(ivf+k+i-1)
                z = z+zr(igeom+3*(i-1)+2)*zr(ivf+k+i-1)
            end do
!
            valpar(1) = x
            valpar(2) = y
            valpar(3) = z
            valpar(4) = zr(itemps)
            call fointe('FM', zk8(icoefh), 4, nompar, valpar, &
                        h, ier)
!
! ---   DETERMINATION DE LA MATRICE MATN DONT LE TERME GENERIQUE
! ---   EST MATN(I,J) = SOMME_LONGUEUR (H*NI(X,Y)*NJ(X,Y).DX.DY) :
!       --------------------------------------------------------
!
!      IMPORTANT: LAMB = CONV * EPAISSEUR
!
!         LAMB=LAMB*LONG*THETA*EP
            h = h*long*theta/2.d0
!
            matn(1, 1) = poi2*h*zr(ivf-1+k+1)**2
            matn(1, 2) = poi2*h*zr(ivf-1+k+1)*zr(ivf-1+k+2)
            matn(2, 1) = matn(1, 2)
            matn(2, 2) = poi2*h*zr(ivf-1+k+2)**2
!
            if (nomte .eq. 'THCOSE3') then
                matn(1, 3) = poi2*h*zr(ivf-1+k+1)*zr(ivf-1+k+3)
                matn(2, 3) = poi2*h*zr(ivf-1+k+2)*zr(ivf-1+k+3)
                matn(3, 1) = matn(1, 3)
                matn(3, 2) = matn(2, 3)
                matn(3, 3) = poi2*h*zr(ivf-1+k+3)**2
            end if
!
            rigith(1, 1) = rigith(1, 1)+matn(1, 1)*matp(1, 1)
            rigith(1, 2) = rigith(1, 2)+matn(1, 1)*matp(1, 2)
            rigith(1, 3) = rigith(1, 3)+matn(1, 1)*matp(1, 3)
            rigith(1, 4) = rigith(1, 4)+matn(1, 2)*matp(1, 1)
            rigith(1, 5) = rigith(1, 5)+matn(1, 2)*matp(1, 2)
            rigith(1, 6) = rigith(1, 6)+matn(1, 2)*matp(1, 3)
!
            rigith(2, 1) = rigith(1, 2)
            rigith(2, 2) = rigith(2, 2)+matn(1, 1)*matp(2, 2)
            rigith(2, 3) = rigith(2, 3)+matn(1, 1)*matp(2, 3)
            rigith(2, 4) = rigith(2, 4)+matn(1, 2)*matp(2, 1)
            rigith(2, 5) = rigith(2, 5)+matn(1, 2)*matp(2, 2)
            rigith(2, 6) = rigith(2, 6)+matn(1, 2)*matp(2, 3)
!
            rigith(3, 1) = rigith(1, 3)
            rigith(3, 2) = rigith(2, 3)
            rigith(3, 3) = rigith(3, 3)+matn(1, 1)*matp(3, 3)
            rigith(3, 4) = rigith(3, 4)+matn(1, 2)*matp(3, 1)
            rigith(3, 5) = rigith(3, 5)+matn(1, 2)*matp(3, 2)
            rigith(3, 6) = rigith(3, 6)+matn(1, 2)*matp(3, 3)
!
            rigith(4, 1) = rigith(1, 4)
            rigith(4, 2) = rigith(2, 4)
            rigith(4, 3) = rigith(3, 4)
            rigith(4, 4) = rigith(4, 4)+matn(2, 2)*matp(1, 1)
            rigith(4, 5) = rigith(4, 5)+matn(2, 2)*matp(1, 2)
            rigith(4, 6) = rigith(4, 6)+matn(2, 2)*matp(1, 3)
!
            rigith(5, 1) = rigith(1, 5)
            rigith(5, 2) = rigith(2, 5)
            rigith(5, 3) = rigith(3, 5)
            rigith(5, 4) = rigith(4, 5)
            rigith(5, 5) = rigith(5, 5)+matn(2, 2)*matp(2, 2)
            rigith(5, 6) = rigith(5, 6)+matn(2, 2)*matp(2, 3)
!
            rigith(6, 1) = rigith(1, 6)
            rigith(6, 2) = rigith(2, 6)
            rigith(6, 3) = rigith(3, 6)
            rigith(6, 4) = rigith(4, 6)
            rigith(6, 5) = rigith(5, 6)
            rigith(6, 6) = rigith(6, 6)+matn(2, 2)*matp(3, 3)
!
            if (nomte .eq. 'THCOSE3') then
!
                rigith(1, 7) = rigith(1, 7)+matn(1, 3)*matp(1, 1)
                rigith(1, 8) = rigith(1, 8)+matn(1, 3)*matp(1, 2)
                rigith(1, 9) = rigith(1, 9)+matn(1, 3)*matp(1, 3)
!
                rigith(2, 7) = rigith(2, 7)+matn(1, 3)*matp(2, 1)
                rigith(2, 8) = rigith(2, 8)+matn(1, 3)*matp(2, 2)
                rigith(2, 9) = rigith(2, 9)+matn(1, 3)*matp(2, 3)
!
                rigith(3, 7) = rigith(3, 7)+matn(1, 3)*matp(3, 1)
                rigith(3, 8) = rigith(3, 8)+matn(1, 3)*matp(3, 2)
                rigith(3, 9) = rigith(3, 9)+matn(1, 3)*matp(3, 3)
!
                rigith(4, 7) = rigith(4, 7)+matn(2, 3)*matp(1, 1)
                rigith(4, 8) = rigith(4, 8)+matn(2, 3)*matp(1, 2)
                rigith(4, 9) = rigith(4, 9)+matn(2, 3)*matp(1, 3)
!
                rigith(5, 7) = rigith(5, 7)+matn(2, 3)*matp(2, 1)
                rigith(5, 8) = rigith(5, 8)+matn(2, 3)*matp(2, 2)
                rigith(5, 9) = rigith(5, 9)+matn(2, 3)*matp(2, 3)
!
                rigith(6, 7) = rigith(6, 7)+matn(2, 3)*matp(3, 1)
                rigith(6, 8) = rigith(6, 8)+matn(2, 3)*matp(3, 2)
                rigith(6, 9) = rigith(6, 9)+matn(2, 3)*matp(3, 3)
!
                rigith(7, 1) = rigith(1, 7)
                rigith(7, 2) = rigith(2, 7)
                rigith(7, 3) = rigith(3, 7)
                rigith(7, 4) = rigith(4, 7)
                rigith(7, 5) = rigith(5, 7)
                rigith(7, 6) = rigith(6, 7)
                rigith(7, 7) = rigith(7, 7)+matn(3, 3)*matp(1, 1)
                rigith(7, 8) = rigith(7, 8)+matn(3, 3)*matp(1, 2)
                rigith(7, 9) = rigith(7, 9)+matn(3, 3)*matp(1, 3)
!
                rigith(8, 1) = rigith(1, 8)
                rigith(8, 2) = rigith(2, 8)
                rigith(8, 3) = rigith(3, 8)
                rigith(8, 4) = rigith(4, 8)
                rigith(8, 5) = rigith(5, 8)
                rigith(8, 6) = rigith(6, 8)
                rigith(8, 7) = rigith(7, 8)
                rigith(8, 8) = rigith(8, 8)+matn(3, 3)*matp(2, 2)
                rigith(8, 9) = rigith(8, 9)+matn(3, 3)*matp(2, 3)
!
                rigith(9, 1) = rigith(1, 9)
                rigith(9, 2) = rigith(2, 9)
                rigith(9, 3) = rigith(3, 9)
                rigith(9, 4) = rigith(4, 9)
                rigith(9, 5) = rigith(5, 9)
                rigith(9, 6) = rigith(6, 9)
                rigith(9, 7) = rigith(7, 9)
                rigith(9, 8) = rigith(8, 9)
                rigith(9, 9) = rigith(9, 9)+matn(3, 3)*matp(3, 3)
            end if
!
        end do
!
    end if
!
! --- RECUPERATION DE LA MATRICE DE RIGIDITE THERMIQUE EN SORTIE DU TE :
!     ----------------------------------------------------------------
    call jevech('PMATTTR', 'E', imattt)
!
! --- AFFECTATION DE LA MATRICE DE RIGIDITE THERMIQUE EN SORTIE DU TE :
!     ---------------------------------------------------------------
    nbddl = 3*nno
    ind = 0
    do i = 1, nbddl
        do j = 1, i
            ind = ind+1
            zr(imattt+ind-1) = rigith(i, j)
        end do
    end do
!
end subroutine
