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

subroutine te0402(option, nomte)
!
    implicit none
!
#include "jeveux.h"
!
#include "asterfort/antisy.h"
#include "asterfort/btdbma.h"
#include "asterfort/btkb.h"
#include "asterfort/btsig.h"
#include "asterfort/cosiro.h"
#include "asterfort/hsall.h"
#include "asterfort/jacbm1.h"
#include "asterfort/jevech.h"
#include "asterfort/jevete.h"
#include "asterfort/jm1dn2.h"
#include "asterfort/jm1dn3.h"
#include "asterfort/promat.h"
#include "asterfort/r8inir.h"
#include "asterfort/sigbar.h"
#include "asterfort/sigvte.h"
#include "asterfort/vectan.h"
#include "asterfort/vectgt.h"
    character(len=16) :: option, nomte
!
! ......................................................................
!     FONCTION :  CALCUL DE LA MATRICE DES CONTRAINTES INITIALES
!                 POUR LE FLAMBEMENT LINEAIRE
!
!                 COQUE_3D
!
!                 OPTION :  RIGI_GEOM
!
!    ARGUMENTS :
!    DONNEES   :       OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
! ......................................................................
!
!
!
!---- DECLARATIONS BIDONS
!
    real(kind=8) :: bid33(3, 3)
!
!
!
!
!---- DECLARATIONS LOCALES
!
    integer(kind=8) :: i, j
    integer(kind=8) :: in
    integer(kind=8) :: ii, jj
    integer(kind=8) :: irig
    integer(kind=8) :: kompt
    integer(kind=8) :: kpgs
!
!
    real(kind=8) :: sigmtd(5)
!
    real(kind=8) :: sigmt(3, 3)
!
    real(kind=8) :: sigma(3, 3)
!
    real(kind=8) :: barsig(9, 9)
!
    real(kind=8) :: vecni(3), antni(3, 3)
!
    real(kind=8) :: veczn(27)
    real(kind=8) :: antzi(3, 3)
!
    real(kind=8) :: rignc(3, 3)
    real(kind=8) :: vri(2601)
!
!
!
!---- DECLARATIONS STANDARDS
!
    integer(kind=8) :: igeom, icontr, imatuu
!
    integer(kind=8) :: lzi, lzr, jcara
!
    integer(kind=8) :: nb1, nb2
!
!
!
!
!---- DECLARATIONS PROPRES COQUE_3D
!
    integer(kind=8) :: inte, intsn
!
    real(kind=8) :: epais
!
    integer(kind=8) :: npge, npgsn, k1
!
    real(kind=8) :: vecta(9, 2, 3)
    real(kind=8) :: vectn(9, 3), vectpt(9, 2, 3)
!
    real(kind=8) :: vectg(2, 3), vectt(3, 3)
!
    real(kind=8) :: jm1(3, 3), detj
!
    real(kind=8) :: hstout(5, 9)
!
    real(kind=8) :: j1dn2(9, 51)
    real(kind=8) :: j1dn3(9, 27)
!
    real(kind=8) :: btild3(5, 27)
!
    parameter(npge=3)
    real(kind=8) :: epsval(npge), ksi3s2, poids(npge)

!
!       POIDS DES POINTS DE GAUSS DANS LA TRANCHE
!
    poids(1) = 0.33333333333333d0
    poids(2) = 1.33333333333333d0
    poids(3) = 0.33333333333333d0
!
!---- RECUPERATION DES POINTEURS ( L : LECTURE, E : ECRITURE )
!
!
!....... GEOMETRIE ( COORDONNEES DES NOEUDS )
!
    call jevech('PGEOMER', 'L', igeom)
!
!....... CONTRAINTES DE CAUCHY ( CONFONDUES AVEC PK2 )
!
!        -- PASSAGE DES CONTRAINTES DANS LE REPERE INTRINSEQUE :
    call cosiro(nomte, 'PCONTRR', 'L', 'UI', 'G', &
                icontr, 'S')
!
!....... MATRICE SYMETRISEE DE RIGIDITE GEOMETRIQUE
!
    call jevech('PMATUUR', 'E', imatuu)
!
!
!
!
!
!---- RECUPERATION DES OBJETS INITIALISES ( SAUF NPGSR )
!
!....... LES ENTIERS
!
    call jevete('&INEL.'//nomte(1:8)//'.DESI', ' ', lzi)
!
!------- NOMBRE DE NOEUDS ( NB1 : SERENDIP , NB2 : LAGRANGE )
!
    nb1 = zi(lzi-1+1)
    nb2 = zi(lzi-1+2)
!
!------- NBRE POINTS INTEGRATIONS ( NPGSR : REDUITE , NPGSN : NORMALE )
!
    npgsn = zi(lzi-1+4)
!
!....... LES REELS ( FONCTIONS DE FORMES, DERIVEES ET POIDS )
!
    call jevete('&INEL.'//nomte(1:8)//'.DESR', ' ', lzr)
!
!
!------ CARACTERISTIQUES DE COQUE
!
    call jevech('PCACOQU', 'L', jcara)
!
    epais = zr(jcara)
!
!       COORDONNEES DES POINTS DE GAUSS DANS LA TRANCHE
    epsval(1) = zr(lzr-1+1251)
    epsval(2) = zr(lzr-1+1252)
    epsval(3) = zr(lzr-1+1253)
!
!       POIDS DES POINTS DE GAUSS DANS LA TRANCHE
!
    poids(1) = 0.33333333333333d0
    poids(2) = 1.33333333333333d0
    poids(3) = 0.33333333333333d0
!
!
!
!
!---- VECTEURS DE BASE AUX NOEUDS
!
    call vectan(nb1, nb2, zr(igeom), zr(lzr), vecta, &
                vectn, vectpt)
!
!
!
!        CALCUL DE LA MATRICE DE RIGIDITE GEOMETRIQUE
!
!            RIG ( 6 * NB1 + 3 , 6 * NB1 + 3 )
!
!        DANS
!
!            VRI ( 6 * NB1 + 3 ) * ( 6 * NB1 + 3 ) )
!
!---- INITIALISATION DE LA MATRICE DE RIGIDITE GEOMETRIQUE
!
    call r8inir(51*51, 0.d0, vri, 1)
!
!---- INITIALISATION DE VECZN AVANT INTEGRATION
!
    call r8inir(27, 0.d0, veczn, 1)
!
!---- COMPTEUR DES POINTS D INTEGRATIONS ( EPAISSEUR * SURFACE )
!
    kpgs = 0
!
!---- BOUCLE SUR LES POINTS D INTEGRATION SUR L EPAISSEUR
!
    do inte = 1, npge
!
!------- COORDONNEE ISOPARAMETRIQUE SUR L EPAISSEUR  DIVISEE PAR DEUX
!
        ksi3s2 = epsval(inte)/2.d0
!
!------- BOUCLE SUR LES POINTS D INTEGRATION SUR LA SURFACE MOYENNE
!
        do intsn = 1, npgsn
!
            kpgs = kpgs+1
            k1 = 6*((intsn-1)*npge+inte-1)
!
!---------- VECTEUR 5 * 1 DES CONTRAINTES LOCALES
!
            sigmtd(1) = zr(icontr-1+k1+1)
            sigmtd(2) = zr(icontr-1+k1+2)
!
            sigmtd(3) = zr(icontr-1+k1+4)
!
            sigmtd(4) = zr(icontr-1+k1+5)
            sigmtd(5) = zr(icontr-1+k1+6)
!
!---------- TENSEUR 3 * 3 CONTRAINTES LOCALES TRIANGULAIRE SUPERIEURE
!
            call sigvte(sigmtd, sigmt)
!
!---------- MATRICE ROTATION GLOBAL --> LOCAL AUX POINTS D INTEGRATION
!
!                              ( T_1 )
!           VECTT ( 3 , 3 ) =  ( T_2 )  = ( LAMDA0 ) T
!                              ( N   )
!
            call vectgt(1, nb1, zr(igeom), ksi3s2, intsn, &
                        zr(lzr), epais, vectn, vectg, vectt)
!
!---------- ROTATION DU TENSEUR DES CONTRAINTES : LOCALES --> GLOBALES
!
!           SIGMA =  ( VECTT ) T * SIGMT * VECTT
!
            call btkb(3, 3, 3, sigmt, vectt, &
                      bid33, sigma)
!
!
!
!
!
!---------- POUR LE TERME NON CLASSIQUE
!
!---------- CALCUL DE    HSTOUT ( 5 , 9 ) = H ( 5 , 6 )  * S ( 6 , 9 )
!
            call hsall(vectt, hstout)
!
!---------- CALCUL DE LA MATRICE JACOBIENNE INVERSE       JM1 ( 3, 3 )
!
            call jacbm1(epais, vectg, vectt, bid33, jm1, &
                        detj)
!
!---------- CALCUL DE
!           J1DN3( 9 , 3 * NB2 )=JTILDM1( 9 , 9 )*DNDQSI3( 9 , 3 * NB2 )
!
            call jm1dn3(nb2, zr(lzr), epais, ksi3s2, intsn, &
                        jm1, j1dn3)
!
!---------- CALCUL DE
!           BTILD3 ( 5 , 27 ) = HSTOUT ( 5 , 9 ) * J1DN3 ( 9 , 3 * NB2 )
!
            call promat(hstout, 5, 5, 9, j1dn3, &
                        9, 9, 3*nb2, btild3)
!
!---------- VECZN ( 27 )  =     INTEGRALE  DE
!           ( BTILD3 ( 5 , 27 ) ) T * SIGMTD ( 5 ) *
!           POIDS SURFACE MOYENNE * DETJ * POIDS EPAISSEUR
!           VOIR ROUTINE INI080 , HSJ1F
!
            call btsig(3*nb2, 5, zr(lzr-1+127+intsn-1)*detj*poids(inte), btild3, &
                       sigmtd, veczn)
!---------- POUR LE TERME CLASSIQUE
!
!---------- BARSIG   ( 9 , 9 )
!
            call sigbar(sigma, barsig)
!
!---------- CALCUL DE
!           J1DN2 ( 9 , 6 * NB1 + 3 ) =
!           JTILDM1 ( 9 , 9 ) * DNDQSI2 ( 9 , 6 * NB1 + 3 )
!
!           INDN = 1 INTEGRATION NORMALE
!           INDC = 1 COMPLET
!
            call jm1dn2(1, 1, nb1, nb2, zr(lzr), &
                        epais, ksi3s2, intsn, vectn, jm1, &
                        j1dn2)
!
!---------- RIG  ( 6 * NB1 + 3 , 6 * NB1 + 3 )  = INTERALE
!           ( J1DN2 ( 9 , 6 * NB1 + 3 ) ) T * BARSIG ( 9 , 9 )
!           *                               J1DN2 ( 9 , 6 * NB1 + 3 ) *
!           POIDS SURFACE MOYENNE * DETJ * POIDS EPAISSEUR
!           VOIR ROUTINE INI080 , HSJ1F
!
            call btdbma(j1dn2, barsig, zr(lzr-1+127+intsn-1)*detj*poids(inte), 9, &
                        6*nb1+3, vri)
!
        end do
    end do
!
!---- PAS DE RIGIDITE DE ROTATION AUTOUR NORMALE
!
!
!
!---- RIGIDITE NON CLASSIQUE
!
!---- BOULE SUR TOUS LES NOEUDS
!
    do in = 1, nb2
!
!------- MATRICE ANTISYMETRIQUE    ANTZI ( 3 , 3 ) AU NOEUD
!
        call antisy(veczn((in-1)*3+1), 1.d0, antzi)
!
!------- NORMALE INITIALE ET SA MATRICE ANTISYM AU NOEUD
!
        do ii = 1, 3
            vecni(ii) = vectn(in, ii)
        end do
!
        call antisy(vecni, 1.d0, antni)
!
!------- RIGIDITE ROTATION NON CLASSIQUE RIGN ( 3 , 3 ) NON SYMETRIQUE
!
        call promat(antzi, 3, 3, 3, antni, &
                    3, 3, 3, rignc)
!
!------- RAJOUT DE LA PARTIE SYMETRIQUE DE RIGN ( 3 , 3 )
!
        if (in .le. nb1) then
!
!---------- NOEUDS DE SERENDIP
            do jj = 1, 3
                do ii = 1, 3
                    j = 6*(in-1)+jj+3
                    i = 6*(in-1)+ii+3
                    irig = (6*nb1+3)*(j-1)+i
                    vri(irig) = vri(irig)+(rignc(ii, jj) &
                                           +rignc(jj, ii))*0.5d0
                end do
            end do
        else
!
!---------- SUPERNOEUD
            do jj = 1, 3
                do ii = 1, 3
                    j = 6*nb1+jj
                    i = 6*nb1+ii
                    irig = (6*nb1+3)*(j-1)+i
                    vri(irig) = vri(irig)+(rignc(ii, jj) &
                                           +rignc(jj, ii))*0.5d0
                end do
            end do
        end if
    end do
!
!
!
!
!______________________________________________________________________
!
!---- STOCKAGE DE LA PARTIE TRIANGULAIRE SUPERIEURE  DANS
!
!                       ZR ( IMATUU )
!
!______________________________________________________________________
!     JEU D INDICES I J POUR LA PARTIE TRIANGULAIRE SUPERIEURE
!
!     ZR ( IMATUU ---> IMATUU + TAILLE  - 1 ) : TRIANGULAIRE SUP DE RIG
!
!     TAILLE = NDDLET * ( 1 + (NDDLET - 1)/2 ) : NDDLET = 6 * NB1 + 3
!
!     VOIR ROUTINE TRANLG
!
!
!
!
!---- COMPTEUR DE POSITION
!
    kompt = 0
!
    do j = 1, 6*nb1+3
        do i = 1, j
            kompt = kompt+1
            zr(imatuu-1+kompt) = +vri((6*nb1+3)*(j-1)+i)
        end do
    end do
!
end subroutine
