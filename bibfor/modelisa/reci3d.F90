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

subroutine reci3d(lirela, mailla, nnoeca, noebe, nbcnx, &
                  cxma, itetra, xbar, immer)
    implicit none
!  DESCRIPTION : DETERMINATION DES RELATIONS CINEMATIQUES ENTRE LES DDLS
!  -----------   D'UN NOEUD DU CABLE ET LES DDLS DES NOEUDS VOISINS DE
!                LA STRUCTURE BETON
!                CAS OU LA STRUCTURE BETON EST MODELISEE PAR DES
!                ELEMENTS 3D
!                APPELANT : IMMECA
!
!  IN     : LIRELA : CHARACTER*19 , SCALAIRE
!                    NOM DE LA SD DE TYPE LISTE_DE_RELATIONS
!  IN     : MAILLA : CHARACTER*8 , SCALAIRE
!                    NOM DU CONCEPT MAILLAGE ASSOCIE A L'ETUDE
!  IN     : NNOECA : CHARACTER*8 , SCALAIRE
!                    NOM DU NOEUD DU CABLE
!  IN     : NOEBE  : INTEGER , SCALAIRE
!                    NUMERO DU NOEUD VOISIN DE LA STRUCTURE BETON LE
!                    PLUS PROCHE DU NOEUD DU CABLE
!  IN     : NBCNX  : INTEGER , SCALAIRE
!                    NOMBRE DE NOEUDS DE LA MAILLE VOISINE DE LA
!                    STRUCTURE BETON
!  IN     : CXMA   : INTEGER , VECTEUR DE DIMENSION AU PLUS NNOMAX
!                    CONTIENT LES NUMEROS DES NOEUDS DE LA MAILLE
!                    VOISINE DE LA STRUCTURE BETON
!                    (TABLE DE CONNECTIVITE)
!  IN     : ITETRA : INTEGER , SCALAIRE
!                    INDICATEUR DU SOUS-DOMAINE TETRAEDRE AUQUEL
!                    APPARTIENT LE NOEUD DU CABLE
!                    ITETRA = 1            SI IMMERSION DANS UNE
!                                          MAILLE TETRAEDRE
!                    ITETRA = 1 OU 2       SI IMMERSION DANS UNE
!                                          MAILLE PYRAMIDE
!                    ITETRA = 1 OU 2 OU 3  SI IMMERSION DANS UNE
!                                          MAILLE PENTAEDRE
!                    ITETRA = 1 OU 2 OU 3  SI IMMERSION DANS UNE
!                          OU 4 OU 5 OU 6  MAILLE HEXAEDRE
!  IN     : XBAR   : REAL*8 , VECTEUR DE DIMENSION 4
!                    COORDONNEES BARYCENTRIQUES DU NOEUD DU CABLE DANS
!                    LE SOUS-DOMAINE TETRAEDRE AUQUEL IL APPARTIENT
!  IN     : IMMER  : INTEGER , SCALAIRE
!                    INDICE D'IMMERSION
!                    IMMER =  0  LE NOEUD DU CABLE EST A L'INTERIEUR
!                                DE LA MAILLE
!                    IMMER = 100 + 10 * NUMERO DE FACE
!                                LE NOEUD DU CABLE EST SUR UNE FACE
!                                DE LA MAILLE
!                    IMMER = 100 + 10 * NUMERO DE FACE + NUMERO D'ARETE
!                                LE NOEUD DU CABLE EST SUR UNE ARETE
!                                DE LA MAILLE
!                    IMMER =  2  LE NOEUD DU CABLE COINCIDE AVEC UN DES
!                                NOEUDS DE LA MAILLE
!
!-------------------   DECLARATION DES VARIABLES   ---------------------
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/afrela.h"
#include "asterfort/ante3d.h"
#include "asterfort/elrfvf.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "asterfort/int_to_char8.h"
#include "asterfort/char8_to_int.h"
!
!
! ARGUMENTS
! ---------
    character(len=19) :: lirela
    character(len=8) :: mailla, nnoeca
    integer(kind=8) :: noebe, nbcnx, cxma(*), itetra, immer
    real(kind=8) :: xbar(*)
!
! VARIABLES LOCALES
! -----------------
    integer(kind=8) :: icnx, iterm, nbsom, nbterm
    integer(kind=8) :: nbtmax, nnomax, noeca
    real(kind=8) :: ksi1, ksi2, ksi3, zero
    complex(kind=8) :: cbid
    character(len=8) :: k8b
    aster_logical :: notlin
!
    real(kind=8) :: ffel3d, ff(27), x(3)
    real(kind=8), pointer :: coemur(:) => null()
    integer(kind=8), pointer :: dimens(:) => null()
    real(kind=8), pointer :: direct(:) => null()
    character(len=8), pointer :: nomddl(:) => null()
    character(len=8), pointer :: nomnoe(:) => null()
!
!-------------------   DEBUT DU CODE EXECUTABLE    ---------------------
!
    call jemarq()
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 1   CREATION DES OBJETS DE TRAVAIL - INITIALISATIONS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    nnomax = 27
    nbtmax = 1+nnomax
    AS_ALLOCATE(vr=coemur, size=nbtmax)
    AS_ALLOCATE(vk8=nomddl, size=nbtmax)
    AS_ALLOCATE(vk8=nomnoe, size=nbtmax)
    AS_ALLOCATE(vi=dimens, size=nbtmax)
    AS_ALLOCATE(vr=direct, size=3*nbtmax)
!
    notlin = (nbcnx .gt. 8)
    if (notlin) then
        if (nbcnx .eq. 10) then
            nbsom = 4
        else if (nbcnx .eq. 13) then
            nbsom = 5
        else if (nbcnx .eq. 15) then
            nbsom = 6
        else
            nbsom = 8
        end if
    else
        nbsom = nbcnx
    end if
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 2   DETERMINATION DE L'ANTECEDENT DU NOEUD DU CABLE DANS L'ELEMENT DE
!     REFERENCE ASSOCIE A L'ELEMENT REEL
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    if (immer .ne. 2) call ante3d(nbsom, itetra, xbar(1), ksi1, ksi2, &
                                  ksi3)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 3   DETERMINATION DES RELATIONS CINEMATIQUES
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    zero = 0.0d0
!
    nomnoe(1) = nnoeca
    nomddl(1) = 'DEPL'
    coemur(1) = 1.0d0
!
! 3.1.1 LE NOEUD DU CABLE COINCIDE TOPOLOGIQUEMENT
!       AVEC UN DES NOEUDS DE LA MAILLE
! ---
    noeca = char8_to_int(nnoeca)
    if (noeca .eq. noebe) goto 60
!
! 3.1.2 LE NOEUD DU CABLE COINCIDE GEOGRAPHIQUEMENT AVEC UN DES
!       NOEUDS DE LA MAILLE : PAS DE RELATIONS CINEMATIQUES
! ---
    if (immer .eq. 2) then
!
        nbterm = 2
        nomnoe(1+1) = int_to_char8(noebe)
        nomddl(1+1) = 'DEPL'
        coemur(1+1) = -1.0d0
!
!
! 3.3 LE NOEUD DU CABLE EST A L'INTERIEUR DE LA MAILLE
! ---
    else
!
        nbterm = 1+nbcnx
        do icnx = 1, nbcnx
            nomnoe(icnx+1) = int_to_char8(cxma(icnx))
            nomddl(icnx+1) = 'DEPL'
            x(1) = ksi1
            x(2) = ksi2
            x(3) = ksi3
            if (nbcnx .eq. 4) then
                call elrfvf('TE4', x, ff)
            else if (nbcnx .eq. 10) then
                call elrfvf('T10', x, ff)
            else if (nbcnx .eq. 5) then
                call elrfvf('PY5', x, ff)
            else if (nbcnx .eq. 13) then
                call elrfvf('P13', x, ff)
            else if (nbcnx .eq. 6) then
                call elrfvf('PE6', x, ff)
            else if (nbcnx .eq. 15) then
                call elrfvf('P15', x, ff)
            else if (nbcnx .eq. 8) then
                call elrfvf('HE8', x, ff)
            else if (nbcnx .eq. 20) then
                call elrfvf('H20', x, ff)
            else if (nbcnx .eq. 27) then
                call elrfvf('H27', x, ff)
            end if
            ffel3d = ff(icnx)
            coemur(icnx+1) = -ffel3d
!            ZR(JCMUR+ICNX) = -FFEL3D(NBCNX,ICNX,KSI1,KSI2,KSI3)
        end do
!
    end if
!
! 3.4 UNE RELATION PAR DDL DE TRANSLATION DU NOEUD DU CABLE
! ---
!.... LE VECTEUR ZI(JDIME) DOIT ETRE REINITIALISE AFIN DE PRENDRE
!.... EN COMPTE LES DIFFERENTS COEFFICIENTS PAR DIRECTION DEFINIS
!.... DANS LE VECTEUR ZR(JDIREC)
!
    do iterm = 1, nbterm
        dimens(iterm) = 3
    end do
!
!.... COEFFICIENTS PAR DIRECTIONS POUR LA PREMIERE RELATION (DDL DX)
!.... PUIS AFFECTATION
!
    do iterm = 1, nbterm
        direct(1+3*(iterm-1)) = 1.0d0
        direct(1+3*(iterm-1)+1) = 0.0d0
        direct(1+3*(iterm-1)+2) = 0.0d0
    end do
!
    call afrela(coemur, [cbid], nomddl, nomnoe, dimens, &
                direct, nbterm, zero, cbid, k8b, &
                'REEL', 'REEL', 0.d0, lirela)
!
!.... COEFFICIENTS PAR DIRECTIONS POUR LA DEUXIEME RELATION (DDL DY)
!.... PUIS AFFECTATION
!
    do iterm = 1, nbterm
        direct(1+3*(iterm-1)) = 0.0d0
        direct(1+3*(iterm-1)+1) = 1.0d0
        direct(1+3*(iterm-1)+2) = 0.0d0
    end do
!
    call afrela(coemur, [cbid], nomddl, nomnoe, dimens, &
                direct, nbterm, zero, cbid, k8b, &
                'REEL', 'REEL', 0.d0, lirela)
!
!.... COEFFICIENTS PAR DIRECTIONS POUR LA TROISIEME RELATION (DDL DZ)
!.... PUIS AFFECTATION
!
    do iterm = 1, nbterm
        direct(1+3*(iterm-1)) = 0.0d0
        direct(1+3*(iterm-1)+1) = 0.0d0
        direct(1+3*(iterm-1)+2) = 1.0d0
    end do
!
    call afrela(coemur, [cbid], nomddl, nomnoe, dimens, &
                direct, nbterm, zero, cbid, k8b, &
                'REEL', 'REEL', 0.d0, lirela)
!
60  continue
!
! --- MENAGE
    AS_DEALLOCATE(vr=coemur)
    AS_DEALLOCATE(vk8=nomddl)
    AS_DEALLOCATE(vk8=nomnoe)
    AS_DEALLOCATE(vi=dimens)
    AS_DEALLOCATE(vr=direct)
    call jedema()
!
! --- FIN DE RECI3D.
end subroutine
