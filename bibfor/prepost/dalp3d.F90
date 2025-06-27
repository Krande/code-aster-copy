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
subroutine dalp3d(nelem, nnoem, degre, nsommx, icnc, &
                  nelcom, numeli, xy, erreur, energi, &
                  volume, alpha, nalpha)
!
!*******************************************************************
!              BUT DE CETTE ROUTINE :                              *
! CALCULER LE DEGRE DE LA SINGULARITE PE EN CHAQUE NOEUD           *
! PUIS EN CHAQUE EF A PARTIR DE PE                                 *
! 1) COMMENT REPERE T-ON UN NOEUD SINGULIER ?                      *
!    L ERREUR LOCALE EST GRANDE EN CE NOEUD                        *
!    COMPAREE A L ERREUR DE REFERENCE SUR TOUTE LA STRUCTURE       *
!    DEROULEMENT DE LA ROUTINE :                                   *
!                                                                  *
!    CALCUL DE L ERREUR DE REFERENCE                               *
!    NOEU1 EST IL SINGULIER ?                                      *
!      CALCUL DE L ERREUR LOCALE SUR LES EFS CONNECTES A NOEU1     *
!      TEST1 : SI ERREUR (NOEU1) > ERREUR GLOBALE                  *
!        NOEU2 = NOEUDS CONNECTES A NOEU1                          *
!        CALCUL DE L ERREUR LOCALE SUR LES EFS CONNECTES A NOEU2   *
!        TEST2 : SI ERREUR (NOEU2) > ERREUR GLOBALE                *
!        LE NOEU1 EST SINGULIER                                    *
!    REMARQUE : CONTRAIREMENT AU 2D NOEU1 EST SINGULIER UNIQUEMENT *
!               SI L UN DE SES VOISINS L EST AUSSI                 *
!               => ON OUBLIE LES NOEUDS SINGULIERS ISOLES          *
! 2) COMMENT CALCULE T-ON PE=ALPHAN(NOEU1) ?                       *
!    SI NOEU1 REGULIER ALPHAN(NOEU1)=DEGRE D INTERPOLATION DE L EF *
!    SI NOEU1 SINGULIER CALCUL DE PE       *
!      ALPHAN(NOEU1)=MIN(ALPHAN(NOEU1),PE) ET                      *
!      ALPHAN(NOEU2)=MIN(ALPHAN(NOEU2),PE)                         *
! 3) COMMENT CALCULE T-ON PE ?                                     *
!    CONSTRUCTION DE LA COURBE ENERGIE EN FONCTION DU RAYON        *
!    EN 3D L ENERGIE EST CALCULEE SUR DES CYLINDRES D EXTREMITE    *
!    NOEU1 ET NOEU2 ET POUR DIFFERENTS RAYONS                      *
!    CALCUL DE PE PAR REFERENCE A L ENERGIE EN POINTE DE FISSURE   *
! 4) REMARQUE : EN 3D ON REGARDE UNIQUEMENT LES NOEUDS SOMMET BORD *
!*******************************************************************
!
! IN  NELEM                  : NOMBRE D ELEMENTS FINIS
! IN  NNOEM                  : NOMBRE DE NOEUDS
! IN  DEGRE                  : DEGRE DES EF (1 POUR P1 2 POUR P2)
! IN  NSOMMX                 : NBRE DE NOEUDS SOMMETS MAX PAR EF
! IN  ICNC(NSOMMX+2,NELEM)   : EF => NOEUDS SOMMETS CONNECTES A EF
!     1ERE VALEUR = NBRE DE NOEUDS SOMMETS CONNECTES A L EF N°X
!     2EME VALEUR = 1 SI EF UTILE 0 SINON
!     CONNECTIVITE  EF N°X=>N° DE NOEUDS SOMMETS CONNECTES A X
!     EN 2D EF UTILE = QUAD OU TRIA
!     EN 3D EF UTILE = TETRA OU HEXA
! IN  NELCOM                 : NBRE D EF SURFACIQUE MAX PAR NOEUD
! IN  NUMELI(NELCOM+2,NNOEM) : NOEUD =>EF SURFACIQUE CONNECTES A NOEUD
!     1ERE VALEUR = NBRE D EFS UTILES CONNECTES AU NOEUD N°X
!     2EME VALEUR = 0 NOEUD MILIEU OU NON CONNECTE A UN EF UTILE
!                   1 NOEUD SOMMET A L INTERIEUR + LIE A UN EF UTILE
!                   2 NOEUD SOMMET BORD + LIE A UN EF UTILE
!     CONNECTIVITE  NOEUD N
! IN  XY(3,NNOEM)            : COORDONNEES DES NOEUDS
! IN  ERREUR(NELEM)          : ERREUR SUR CHAQUE EF
! IN  ENERGI(NELEM)          : ENERGIE SUR CHAQUE EF
! IN  VOLUME(NELEM)          : VOLUME DE CHAQUE EF
! OUT ALPHA(NELEM)           : DEGRE DE LA SINGULARITE PAR ELEMENT
! OUT NALPHA                 : NOMBRE DE CPE PAR ELEMENT DIFFERENTS
!                              1 PAR DEFAUT SI PAS DE SINGULARITE
!
! aslint: disable=W1306
    implicit none
!
! DECLARATION GLOBALE
!
#include "asterfort/dfort3.h"
#include "asterfort/dzonfg.h"
    integer(kind=8) :: nelem, nnoem, degre, nsommx, nelcom
    integer(kind=8) :: icnc(nsommx+2, nelem), numeli(nelcom+2, nnoem)
    integer(kind=8) :: nalpha
    real(kind=8) :: xy(3, nnoem), erreur(nelem), energi(nelem), volume(nelem)
    real(kind=8) :: alpha(nelem)
!
! DECLARATION LOCALE
!
    integer(kind=8) :: i, inno, inel, nuef, noeu1, noeu2
    integer(kind=8) :: tbnozo(1000), nbnozo(3), tbelzo(1000), nbelzo(3)
    integer(kind=8) :: nbnoe
    real(kind=8) :: precmo, precre, prec, vol
    real(kind=8) :: dtyp, alphan(nnoem), pe
!
! 1 - DEGRE D INTERPOLATION
!
    if (degre .eq. 1) then
        dtyp = 1.d+0
    else
        dtyp = 2.d+0
    end if
!
! 2 - INITIALISATION DES ALPHA = DEGRE D INTERPOLATION
!
    do inno = 1, nnoem
        alphan(inno) = dtyp
    end do
    do inel = 1, nelem
        alpha(inel) = dtyp
    end do
!
! 3 - CALCUL DES PRECISIONS MOYENNE ET DE REFERENCE
!
    precmo = 0.d+0
    vol = 0.d+0
    do inel = 1, nelem
        precmo = precmo+erreur(inel)**2
        vol = vol+volume(inel)
    end do
    precmo = sqrt(precmo/vol)
!
    if (dtyp .eq. 1.0d0) then
        precre = 3.d0*precmo
    else if (dtyp .eq. 2.0d0) then
        precre = 2.d0*precmo
    end if
!
! 4 - BOUCLE SUR LES NOEUDS SOMMET BORD POUR DETECTER
!     SI LE NOEUD CONSIDERE EST SINGULIER OU PAS
!
    do inno = 1, nnoem
!
        if (numeli(2, inno) .ne. 2) goto 40
!
! 4.1 - ERREUR LOCALE SUR COUCHE 1 (=EF CONNECTES A INNO)
!
        prec = 0.0d0
        vol = 0.0d0
!
        do inel = 1, numeli(1, inno)
            nuef = numeli(2+inel, inno)
            prec = prec+erreur(nuef)**2
            vol = vol+volume(nuef)
        end do
        if (prec .ne. 0.d+0 .and. vol .ne. 0.d+0) then
            prec = sqrt(prec/vol)
        else
            prec = 0.d+0
        end if
!
! 4.2 - TEST1 : SI PREC > PRECREF ON VA CHERCHER LES NOEUDS NOEU2
!       CONNECTES A INNO
!
        if (prec .ge. precre) then
!
! 4.2.1 - RECHERCHE DES NOEUDS ET EFS COMPOSANTS LES COUCHES 1,2 ET 3
!
            call dzonfg(nsommx, icnc, nelcom, numeli, inno, &
                        tbelzo, nbelzo, tbnozo, nbnozo)
            noeu1 = inno
!
! 4.2.2 - TEST2 : POUR CHAQUE NOEU2 CONNECTE A NOEU1
!         ON DETECTE SI NOEU2 EST SINGULIER
!         SI OUI ON DETERMINE PE ET UNIQUEMENT SI NOEU2>NOEU1
!         POUR EVITER DE FAIRE DEUX FOIS LE CALCUL
!
            do i = 1, nbnozo(1)
                noeu2 = tbnozo(i)
                if (numeli(2, noeu2) .ne. 2) goto 60
                if (noeu2 .gt. noeu1) then
                    vol = 0.0d0
                    prec = 0.0d0
                    do inel = 1, numeli(1, noeu2)
                        nuef = numeli(2+inel, noeu2)
                        prec = prec+erreur(nuef)**2
                        vol = vol+volume(nuef)
                    end do
                    if (prec .ne. 0.d+0 .and. vol .ne. 0.d+0) then
                        prec = sqrt(prec/vol)
                    else
                        prec = 0.d+0
                    end if
!
                    if (prec .ge. precre) then
!
! 4.2.2.1 - CALCUL DE PE
!
                        nbnoe = nbnozo(1)+nbnozo(2)+nbnozo(3)
                        call dfort3(nsommx, icnc, noeu1, noeu2, tbelzo, &
                                    nbelzo(3), tbnozo, nbnoe, xy, volume, &
                                    energi, pe)
!
                        if (alphan(noeu1) .gt. pe) alphan(noeu1) = pe
                        if (alphan(noeu2) .gt. pe) alphan(noeu2) = pe
                    end if
                end if
60              continue
            end do
        end if
40      continue
    end do
!
! 5 - ON REMPLIT LE TABLEAU ALPHA = DEGRE DE LA SINGULARITE PAR EF
!
    nalpha = 1
    do inel = 1, nelem
        if (icnc(2, inel) .lt. 1) goto 110
        do inno = 1, icnc(1, inel)
            alpha(inel) = min(alpha(inel), alphan(icnc(inno+2, inel)))
        end do
        if (alpha(inel) .lt. dtyp) then
            nalpha = nalpha+1
        end if
110     continue
    end do
!
end subroutine
