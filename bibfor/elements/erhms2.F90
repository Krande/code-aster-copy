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
subroutine erhms2(ino, nbs, theta, jac, &
                  nx, ny, sielnp, adsip, sielnm, &
                  nbcmp, typmav, tbref1, tbref2, ivois, &
                  tm2h1s)
!
! =====================================================================
!  ERREUR EN HYDRO-MECANIQUE - TERME DE SAUT - DIMENSION 2
!  **        *     *                    *                *
! =====================================================================
!  - FONCTION REALISEE :  CALCUL DE L'ERREUR DUE AUX TERMES DE SAUT
!                         DE LA MECANIQUE ET DE L'HYDRAULIQUE DANS
!                         L'INDICATEUR HM PERMANENT.
!
!     ARGUMENTS:
!     ----------
!
!      ENTREE :
!-------------
! IN   PERMAN : PERMANENT OU INSTATIONNAIRE ?
! IN   INO    : NUMERO DE L'ARETE
! IN   NBS    : NOMBRE DE SOMMETS SUR L'ARETE = NOMBRE D'ARETES
! IN   THETA  : PARAMETRE THETA DE DISCRETISATION TEMPORELLE
! IN   JAC    : VECTEUR DES JACOBIENS DE LA TRANSFORMATION AUX NOEUDS
! IN   NX     : VECTEUR DES ABSCISSES DES NORMALES AUX NOEUDS
! IN   NY     : VECTEUR DES ORDONNEES DES NORMALES AUX NOEUDS
! IN   SIELNP : CONTRAINTES AUX NOEUDS PAR ELEMENT A L'INSTANT ACTUEL
! IN   ADSIP  : ADRESSE DANS ZR DU TABLEAU DES CONTRAINTES DE PRESSION
!               DE LA MECANIQUE
! IN   SIELNM : CONTRAINTES AUX NOEUDS PAR ELEMENT A L'INSTANT PRECEDENT
! IN   NBCMP  : NOMBRE DE COMPOSANTES DU VECTEUR DES CONTRAINTES
!               GENERALISEES
! IN   TYPMAV : TYPE DE LA MAILLE VOISINE
! IN   TBREF1 : TABLEAU AUXILIAIRE CONTENANT DES ADRESSES/VALEURS
! IN   TBREF2 : TABLEAU AUXILIAIRE 2 CONTENANT DES ADRESSES/VALEURS
! IN   IVOIS  : ADRESSE DANS ZI DES VOISINS
!
!      SORTIE : SAUT AUX INTERFACES DES TERMES DIFFUSIFS
!               DES INDICATEURS HM
!-------------
! OUT TM2H1S : TABLEAU CONTENANT LES TERMES DE SAUT DES TERMES DIFFUSIFS
!              (2 POUR LA MECANIQUE, 1 POUR L'HYDRAULIQUE)
!  1 : MECANIQUE
!  2 : DERIVEE TEMPORELLE DE LA MECA
!  3 : HYDRAULIQUE
! ......................................................................
!
    implicit none
!
! DECLARATION PARAMETRES D'APPELS
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/indiis.h"
    integer(kind=8) :: ino, nbs
    integer(kind=8) :: ivois, nbcmp, adsip
    integer(kind=8) :: tbref1(12), tbref2(12)
    character(len=8) :: typmav
    real(kind=8) :: theta, jac(3), nx(3), ny(3)
    real(kind=8) :: sielnp(140), sielnm(140), tm2h1s(3)
!
!
!
!
! DECLARATION VARIABLES LOCALES
!
    real(kind=8) :: sixxp(3), siyyp(3), sixyp(3)
    real(kind=8) :: sixxm(3), siyym(3), sixym(3)
    real(kind=8) :: sixxpv(3), siyypv(3), sixypv(3)
    real(kind=8) :: sixxmv(3), siyymv(3), sixymv(3)
    real(kind=8) :: dsgxxp(3), dsgyyp(3), dsgxyp(3)
    real(kind=8) :: dsgxxm(3), dsgyym(3), dsgxym(3)
    real(kind=8) :: fluxp(3), fluyp(3), fluxm(3), fluym(3)
    real(kind=8) :: fluxpv(3), fluypv(3), fluxmv(3), fluymv(3)
    real(kind=8) :: dmxp(3), dmyp(3), dmxm(3), dmym(3)
    real(kind=8) :: dmxfin(3), dmyfin(3)
!
    real(kind=8) :: intme1(3), intme2(3), inthyd(3)
    real(kind=8) :: ta1
!
 integer(kind=8) :: iarepe, jceldp, jcelvp, jceldm, jcelvm, imav, igrel, iel, iavalp, iavalm, iconx1
    integer(kind=8) :: iconx2, admavp, admavm, adinov, adjnov, admnov
    integer(kind=8) :: jad, jadv, ncher
    integer(kind=8) :: nbnv, i, jno, mno, inov, jnov, mnov
    integer(kind=8) :: iaux
    character(len=2) :: formv, noeuv
! =====================================================================
! 1.  RECUPERATION SUR LA MAILLE COURANTE AUX NOEUDS INO ET JNO DE :
!     . CONTRAINTES EFFECTIVES (SIGMA MECANIQUE : SIXX, SIYY, SIXY)
!     . FLUX HYDRAULIQUE M
!
!              X1          X2          X3
!               O-----------O-----------O
!              INO         MNO         JNO
!
!         POINTS  1 --> INO PREMIER POINT DE L'ARETE COURANTE
!                 2 --> JNO DEUXIEME POINT  DE L'ARETE COURANTE
!                 3 --> MNO NOEUD MILIEU S'IL EXISTE
! =====================================================================
!
    ta1 = 1.d0-theta
!
    if (ino .eq. nbs) then
        jno = 1
    else
        jno = ino+1
    end if
!
    iaux = nbcmp*(ino-1)
    sixxp(1) = sielnp(iaux+1)
    siyyp(1) = sielnp(iaux+2)
    sixyp(1) = sielnp(iaux+4)
    fluxp(1) = sielnp(iaux+adsip+1+1+5)
    fluyp(1) = sielnp(iaux+adsip+1+2+5)
!
    iaux = nbcmp*(jno-1)
    sixxp(2) = sielnp(iaux+1)
    siyyp(2) = sielnp(iaux+2)
    sixyp(2) = sielnp(iaux+4)
    fluxp(2) = sielnp(iaux+adsip+1+1+5)
    fluyp(2) = sielnp(iaux+adsip+1+2+5)
!
! BIEN QUE LA PRESSION NE SOIT CALCULEE QU'AUX NOEUDS SOMMETS, LE
! GRADIENT L'EST EGALEMENT AU NOEUD MILIEU, VIA EPSTHM
!
    mno = nbs+ino
!
    iaux = nbcmp*(mno-1)
    sixxp(3) = sielnp(iaux+1)
    siyyp(3) = sielnp(iaux+2)
    sixyp(3) = sielnp(iaux+4)
!
!
! BIEN QUE LA PRESSION NE SOIT CALCULEE QU'AUX NOEUDS SOMMETS, LE
! GRADIENT L'EST EGALEMENT AU NOEUD MILIEU, VIA EPSTHM
!
    fluxp(3) = sielnp(iaux+adsip+1+1+5)
    fluyp(3) = sielnp(iaux+adsip+1+2+5)
!
    iaux = nbcmp*(ino-1)
    sixxm(1) = sielnm(iaux+1)
    siyym(1) = sielnm(iaux+2)
    sixym(1) = sielnm(iaux+4)
    fluxm(1) = sielnm(iaux+adsip+1+1+5)
    fluym(1) = sielnm(iaux+adsip+1+2+5)
!
    iaux = nbcmp*(jno-1)
    sixxm(2) = sielnm(iaux+1)
    siyym(2) = sielnm(iaux+2)
    sixym(2) = sielnm(iaux+4)
    fluxm(2) = sielnm(iaux+adsip+1+1+5)
    fluym(2) = sielnm(iaux+adsip+1+2+5)
!
! BIEN QUE LA PRESSION NE SOIT CALCULEE QU'AUX NOEUDS SOMMETS, LE
! GRADIENT L'EST EGALEMENT AU NOEUD MILIEU, VIA EPSTHM
!
    mno = nbs+ino
!
    iaux = nbcmp*(mno-1)
    sixxm(3) = sielnm(iaux+1)
    siyym(3) = sielnm(iaux+2)
    sixym(3) = sielnm(iaux+4)
    fluxm(3) = sielnm(iaux+adsip+1+1+5)
    fluym(3) = sielnm(iaux+adsip+1+2+5)

!
!
! =====================================================================
! 2. CARACTERISATION DE LA MAILLE VOISINE
! =====================================================================
!
! 2.1. ---- RECHERCHE DES ADRESSES POUR OBTENIR SIGMA SUR LE VOISIN ---
!
    iarepe = tbref1(1)
    jceldp = tbref1(2)
    jcelvp = tbref1(3)
    jceldm = tbref2(2)
    jcelvm = tbref2(3)
!
! ADRESSE DE LA MAILLE VOISINE PARTAGEANT L'ARETE INO
    imav = zi(ivois+ino)
!
! NUMERO DU LIGREL DE LA MAILLE VOISINE DE NUMERO GLOBAL IMAV
    igrel = zi(iarepe+2*(imav-1))
!
! INDICE DE LA MAILLE VOISINE DANS LE IGREL
    iel = zi(iarepe+2*(imav-1)+1)
!
! ADRESSE (DANS .CELV DE SIELNP) DU DEBUT DU GREL IGREL
    iavalp = jcelvp-1+zi(jceldp-1+zi(jceldp-1+4+igrel)+8)
! ADRESSE (DANS .CELV DE SIELNM) DU DEBUT DU GREL IGREL
    iavalm = jcelvm-1+zi(jceldm-1+zi(jceldm-1+4+igrel)+8)
!
! 2.2. ----- TESTS SUR LA MAILLE VOISINE ------------------------------
!
    formv = typmav(1:2)
    noeuv = typmav(5:5)
!
    if (formv .eq. 'TR') then
        if (noeuv .eq. '3') then
            nbnv = 3
        else
            nbnv = 6
        end if
    else if (formv .eq. 'QU') then
        if (noeuv .eq. '4') then
            nbnv = 4
        else if (noeuv .eq. '8') then
            nbnv = 8
        else
            nbnv = 9
        end if
    end if
!
! 2.3 ----- ADRESSE DANS LE VECTEUR DES CONTRAINTES DE LA PREMIERE
!           COMPOSANTE DU 1ER NOEUD DE L'ELEMENT IEL DE NUMERO IGREL
!
    admavp = iavalp-1+nbcmp*nbnv*(iel-1)
    admavm = iavalm-1+nbcmp*nbnv*(iel-1)
!
! ADRESSE DE LA COLLECTION CONNECTIVITE
!
    iconx1 = tbref1(11)
!
! ADRESSE DU POINTEUR DE LONGUEUR DE LA CONNECTIVITE
    iconx2 = tbref1(12)
!
    jad = iconx1-1+zi(iconx2+zi(ivois)-1)
    jadv = iconx1-1+zi(iconx2+zi(ivois+ino)-1)
!
!     NCHER : NUMERO GLOBAL DU NOEUD LOCAL INO DANS LA MAILLE COURANTE
!
    ncher = zi(jad-1+ino)
!
!     ON CHERCHE A QUEL NUMERO LOCAL (INOV) DANS LA MAILLE VOISINE
!     CORRESPOND LE NOEUD DE NUMERO GLOBAL NCHER
!
    inov = indiis(zi(jadv), ncher, 1, nbnv)
!
!     IDEM POUR JNO LOCAL DE LA MAILLE COURANTE
!
    ncher = zi(jad-1+jno)
    jnov = indiis(zi(jadv), ncher, 1, nbnv)
!
! ADRESSE DES NOEUDS INOV ET JNOV DE LA MAILLE VOISINE
! DE L'ELEMENT COURANT
!
    adinov = nbcmp*(inov-1)
    adjnov = nbcmp*(jnov-1)
!
! CONTRAINTES AUX NOEUDS SUR L'ELEMENT VOISIN, POUR LES 2 NOEUDS DE
! L'ARETE
!
    iaux = admavp+adinov
    sixxpv(1) = zr(iaux+1)
    siyypv(1) = zr(iaux+2)
    sixypv(1) = zr(iaux+4)
    fluxpv(1) = zr(iaux+adsip+1+1+5)
    fluypv(1) = zr(iaux+adsip+1+2+5)
!
    iaux = admavp+adjnov
    sixxpv(2) = zr(iaux+1)
    siyypv(2) = zr(iaux+2)
    sixypv(2) = zr(iaux+4)
    fluxpv(2) = zr(iaux+adsip+1+1+5)
    fluypv(2) = zr(iaux+adsip+1+2+5)
!
! 2.4. NOEUD MILIEU
!
    mno = nbs+ino
!
! MEME OPERATION QUE POUR LES EXTREMITES : ON PASSE PAR LE NUMERO GLOBAL
!
    ncher = zi(jad-1+mno)
    mnov = indiis(zi(jadv), ncher, 1, nbnv)
!
! ADRESSE DU NOEUD MNOV DE LA MAILLE VOISINE DE L'ELEMENT COURANT
!
    admnov = nbcmp*(mnov-1)
!
    iaux = admavp+admnov
    sixxpv(3) = zr(iaux+1)
    siyypv(3) = zr(iaux+2)
    sixypv(3) = zr(iaux+4)
    fluxpv(3) = zr(iaux+adsip+1+1+5)
    fluypv(3) = zr(iaux+adsip+1+2+5)
!
! CONTRAINTES MECANIQUES ET HYDRAULIQUES AUX NOEUDS
! SUR L'ELEMENT VOISIN, POUR LES 2 NOEUDS DE L'ARETE
!
    iaux = admavm+adinov
!
    sixxmv(1) = zr(iaux+1)
    siyymv(1) = zr(iaux+2)
    sixymv(1) = zr(iaux+4)
    fluxmv(1) = zr(iaux+adsip+1+1+5)
    fluymv(1) = zr(iaux+adsip+1+2+5)
!
    iaux = admavm+adjnov
!
    sixxmv(2) = zr(iaux+1)
    siyymv(2) = zr(iaux+2)
    sixymv(2) = zr(iaux+4)
    fluxmv(2) = zr(iaux+adsip+1+1+5)
    fluymv(2) = zr(iaux+adsip+1+2+5)
!
! NOEUD MILIEU
!
    mno = nbs+ino
!
! MEME OPERATION QUE POUR LES EXTREMITES : ON PASSE PAR LE NUMERO GLOBAL
!
    ncher = zi(jad-1+mno)
    mnov = indiis(zi(jadv), ncher, 1, nbnv)
!
! ADRESSE DU NOEUD MNOV DE LA MAILLE VOISINE DE L'ELEMENT COURANT
!
    admnov = nbcmp*(mnov-1)
!
    iaux = admavm+admnov
!
    sixxmv(3) = zr(iaux+1)
    siyymv(3) = zr(iaux+2)
    sixymv(3) = zr(iaux+4)
    fluxmv(3) = zr(iaux+adsip+1+1+5)
    fluymv(3) = zr(iaux+adsip+1+2+5)

!
! =====================================================================
! 3. --- CALCUL DES SAUTS DES TERMES DIFFUSIFS
! =====================================================================
!
! =====================================================================
! 3.A --- PARTIE MECANIQUE
! =====================================================================
!
! 3.1. SAUTS EN CHAQUE NOEUD DE L'ARETE
!
    do 10, i = 1, 3
!
        dsgxxp(i) = sixxp(i)-sixxpv(i)
        dsgyyp(i) = siyyp(i)-siyypv(i)
        dsgxyp(i) = sixyp(i)-sixypv(i)
!
        dmxp(i) = fluxp(i)-fluxpv(i)
        dmyp(i) = fluyp(i)-fluypv(i)
!
! POUR L'INTEGRATION NUMERIQUE, ON UTILISE :
!      . SI 3 NOEUDS SOMMETS, UNE FORMULATION DE SIMPSON
!      . SI 2 NOEUDS SOMMETS, UNE METHODE DES TRAPEZES
!
        intme1(i) = jac(i)*((dsgxxp(i)*nx(i)+dsgxyp(i)*ny(i))**2+(dsgxyp(i)*nx(i)+dsgyyp(i)*ny&
                    &(i))**2)
!
        dmxfin(i) = theta*dmxp(i)
        dmyfin(i) = theta*dmyp(i)
!
        dsgxxm(i) = sixxm(i)-sixxmv(i)
        dsgyym(i) = siyym(i)-siyymv(i)
        dsgxym(i) = sixym(i)-sixymv(i)
!
        dmxm(i) = fluxm(i)-fluxmv(i)
        dmym(i) = fluym(i)-fluymv(i)
!
        intme2(i) = jac(i)*(((dsgxxp(i)-dsgxxm(i))*nx(i)+(dsgxyp(i)-dsgxym(i))*ny(i))**2 &
                    &+((dsgxyp(i)-dsgxym(i))*nx(i)+(dsgyyp(i)-dsgyym(i))*ny(i))**2)
!
        dmxfin(i) = theta*dmxp(i)+ta1*dmxm(i)
        dmyfin(i) = theta*dmyp(i)+ta1*dmym(i)

!
! =====================================================================
! 3.B --- PARTIE HYDRAULIQUE
! =====================================================================
!
        inthyd(i) = jac(i)*((dmxfin(i)*nx(i)+dmyfin(i)*ny(i))**2)
!
10  end do
!
! =====================================================================
! 3.C --- ASSEMBLAGE DES DIFFERENTS TERMES
! =====================================================================
!
    tm2h1s(1) = tm2h1s(1)+((intme1(1)+4.d0*intme1(3)+intme1(2))/3.d0)
    tm2h1s(3) = tm2h1s(3)+((inthyd(1)+4.d0*inthyd(3)+inthyd(2))/3.d0)
    tm2h1s(2) = tm2h1s(2)+(intme2(1)+4.d0*intme2(3)+intme2(2))/3.d0

!
end subroutine
