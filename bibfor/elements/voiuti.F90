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
subroutine voiuti(numa, codvoi, nvoima, nscoma, iarepe, &
                  iaddvo, iadvoi, nbvois, livois, tyvois, &
                  nbnovo, nbsoco, lisoco)
    implicit none
!   IN  NUMA     : NUMERO DE MAILLE DU MAILLAGE
!       CODVOI   :  CODE DONNANT LES COMBINAISONS
!                   DES VOISINAGES POSSIBLES
!                   3D PAR FACE    : F3 : 1
!                   2D PAR FACE    : F2 : 2
!                   3D PAR ARRETE  : A3 : 3
!                   2D PAR ARRETE  : A2 : 4
!                   1D PAR ARRETE  : A1 : 5
!                   3D PAR SOMMET  : S3 : 6
!                   2D PAR SOMMET  : S2 : 7
!                   1D PAR SOMMET  : S1 : 8
!                   0D PAR SOMMET  : S0 : 9
!       NVOIMA   :  NOMBRE MAX DE VOISINS POSSIBLES
!       NSCOMA   :  NOMBRE MAX DE SOMMETS COMMUNS
!       IAREPE   :  ADRESSE JEVEUX D UN OPJET DE TYPE REPE
!       IADDVO   :  ADRESSE JEVEUX DU POINEUR DE VALEURS DANS VGE
!       IADVOI   :  ADRESSE JEVEUX DES VALEURS DE VGE
!  OUT
!       NBVOIS   :  NOMBRE DE VOISINS
!       LIVOIS   :  LISTE DE CES VOISINS (NUM. DE MAILLES AFFECTEES)
!       TYVOIS   :  TYPE(INTEGER) DES VOISINS
!       NBNOVO   :  NOMBRE DE NOEUDS DE CHACUN DE CES VOISINS
!       NBSOCO   :  POUR CHAQUE VOISINS NOMBRE DE SOMMETS PARTAGES
!       LISOCO   :  LISTE DE CES SOMMETS
!                  (IMA,IS,1) EN NUMEROTATION LOCALE MAILLE NUMA
!                  (IMA,IS,2) EN NUMEROTATION LOCALE MAILLE VOISINE
#include "jeveux.h"
#include "asterfort/lxlgut.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: numa, nvoima, nscoma, iarepe, iaddvo, iadvoi, nbvois
!     PARAMETER(NVOIMA=100,NSCOMA=4)
    integer(kind=8) :: livois(1:nvoima), tyvois(1:nvoima), nbnovo(1:nvoima)
    integer(kind=8) :: nbsoco(1:nvoima), lisoco(1:nvoima, 1:nscoma, 1:2)
    character(len=*) :: codvoi
    integer(kind=8) :: iv, is, iel, ideb, ifin, icode, jcode, numav, ielv, typev
    integer(kind=8) :: ntymax
    parameter(ntymax=3)
    integer(kind=8) :: lcod, ityvo, ntyvo, lityvo(1:ntymax)
    character(len=2) :: tybase(9)
!-----------FONCTIONS  D ACCES A VGE -----------------------------------
!     IADDVO : ADRESSE JEVEUX DU TABLEAU DE POINTEURS DANS LA SD EL_VOIS
!     IADVOI : ADRESSE JEVEUX DE LA SD EL_VOIS
!
!     DES DONNEES DES VOISINS DE LA MAILLE NUMA (0 SI MAILLE PAS ACTIVE)
#define zzadvo(numa) zi(iaddvo+numa-1)+iadvoi
!
!     NOMBBRE DE VOISINS DE NUMA EME MAILLE
#define zznbvo(numa) zi(zzadvo(numa)-1+1)
!
!     POUR LA MAILLE NUMA
!     POUR LE VOISIN IV
!     ADRESSE DES DONNEES
! FAUX ???      ZZADVE(NUMA,IV) = ZI(ZZADVO(NUMA)-1+1+IV)+IADVOI-1
#define zzadve(numa,iv) zi(zzadvo(numa)-1+1+iv)+zzadvo(numa)-1
!
!     POUR LA MAILLE NUMA
!     POUR LE VOISIN IV
!     TYPE DE VOISINAGE :
!        3D PAR FACE    : F3 : 1
!        2D PAR FACE    : F2 : 2
!        3D PAR ARRETE  : A3 : 3
!        2D PAR ARRETE  : A2 : 4
!        1D PAR ARRETE  : A1 : 5
!        3D PAR SOMMET  : S3 : 6
!        2D PAR SOMMET  : S2 : 7
!        1D PAR SOMMET  : S1 : 8
!        0D PAR SOMMET  : S0 : 9
#define zztyvo(numa,iv) zi(zzadve(numa,iv)-1+1)
!
!     POUR LA MAILLE NUMA
!     POUR LE VOISIN IV
!     NUMERO DE MAILLE
#define zzmavo(numa,iv) zi(zzadve(numa,iv)-1+2)
!
!     POUR LA MAILLE NUMA
!     POUR LE VOISIN IV
!     NOMBRE DE NOEUDS DE MAILLE
#define zznbno(numa,iv) zi(zzadve(numa,iv)-1+3)
!
!
!     POUR LA MAILLE NUMA
!     POUR LE VOISIN IV
!        NOMBRE DE SOMMETS COMMUNS
#define zznbsc(numa,iv) zi(zzadve(numa,iv)-1+4)
!
!
!     POUR LA MAILLE NUMA
!     POUR LE VOISIN IV
!     POUR LE SOMMET COMMUN IS
!     NUMERO LOCAL DANS NUMA
#define zzloc1(numa,iv,is) zi(zzadve(numa,iv)-1+4+1+2*(is-1))
!
!
!
!     POUR LA MAILLE NUMA
!     POUR LE VOISIN IV
!    POUR LE SOMMET COMMUN IS
!     NUMERO LOCAL DANS IV
#define zzloc2(numa,iv,is) zi(zzadve(numa,iv)-1+4+1+2*(is-1)+1)
!-----------FIN FONCTIONS  D ACCES A VGE -------------------------------
    data tybase/'F3', 'F2', 'A3', 'A2', 'A1', 'S3', 'S2', 'S1', 'S0'/
!
!  1 CETTE MAILLE EST ELLE AFFECTEE DANS LE MODELE ?
!
    iel = zi(iarepe-1+2*(numa-1)+2)
    if (iel .eq. 0) then
        nbvois = 0
        goto 80
!
    end if
!
!
!  1 RECHERCHE DES TYPES DES VOISINAGE ATTENDUS
!
    ntyvo = 0
    lcod = lxlgut(codvoi)
    if (lcod .gt. 2*ntymax) then
        call utmess('F', 'VOLUFINI_7', sk=codvoi, si=lcod)
    end if
    do icode = 1, lcod/2
        ideb = 2*(icode-1)+1
        ifin = 2*icode
        do jcode = 1, 9
            if (codvoi(ideb:ifin) .eq. tybase(jcode)) then
                ntyvo = ntyvo+1
                lityvo(ntyvo) = jcode
                goto 20
!
            end if
        end do
        call utmess('F', 'VOLUFINI_6', sk=codvoi(ideb:ifin))
20      continue
    end do
    if (ntyvo .eq. 0) then
        nbvois = 0
        goto 80
!
    end if
!
!      REMPLISSAGE DES TABLEAUX
!
    nbvois = 0
    do iv = 1, zznbvo(numa)
        numav = zzmavo(numa, iv)
        ielv = zi(iarepe-1+2*(numav-1)+2)
!
!  LE VOISIN EST IL UNE MAILLE AFFECTEE DANS LE MODELE
!
        if (ielv .eq. 0) then
            goto 70
!
        end if
        typev = zztyvo(numa, iv)
!
!  LE VOISIN EST IL D UN TYPE ATTENDU
!
        do ityvo = 1, ntyvo
            if (typev .eq. lityvo(ityvo)) then
                goto 50
!
            end if
        end do
        goto 70
!
50      continue
        nbvois = nbvois+1
        livois(nbvois) = numav
        tyvois(nbvois) = typev
        nbnovo(nbvois) = zznbno(numa, iv)
        nbsoco(nbvois) = zznbsc(numa, iv)
        do is = 1, zznbsc(numa, iv)
            lisoco(nbvois, is, 1) = zzloc1(numa, iv, is)
            lisoco(nbvois, is, 2) = zzloc2(numa, iv, is)
        end do
70      continue
    end do
80  continue
!
end subroutine
