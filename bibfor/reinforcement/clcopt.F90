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

subroutine clcopt(nb, fcttab, atab, ax, ay)
!_____________________________________________________________________
!
! CLCOPT
!
! RECHERCHE DU MINIMUM D'ARMATURES (METHODE DE CAPRA ET MAURY)
!
! I NB         NOMBRE DE FACETTES
! I FCTTAB     PARAMETRES TRIGONOMETRIQUES POUR LES (NB) FACETTES
! I ATAB       SECTIONS DES ARMATURES (DENSITES) CALCULEES PAR FACETTE
! O AX         SECTION DES ARMATURES EN X
! O AY         SECTION DES ARMATURES EN Y
!
!_____________________________________________________________________
!
!
    implicit none
!
    integer(kind=8) :: nb
    real(kind=8) :: fcttab(nb, 3)
    real(kind=8) :: atab(nb)
    real(kind=8) :: ax
    real(kind=8) :: ay
!
    real(kind=8) :: precis
    parameter(precis=1d-10)
!
!       DIMENSIONNEMENT DES TABLEAUX DE POINTS DU POLYGONE
!       BORDANT LE DOMAINE DE VALIDITE
    integer(kind=8) :: plgnnb
!
!         ******   POLYGONE BORDANT LE DOMAINE DE VALIDITE  ******
!           TABLEAUX DES POSITIONS (SURDIMENSIONNES)
    real(kind=8) :: xp(nb*2)
    real(kind=8) :: yp(nb*2)
!
!       ANGLES THETA
    integer(kind=8) :: ap(nb*2)
!
!       TABLEAU DE SUIVI DES INDICES
!       (CONTIENT A CHAQUE INDICE LE SUIVANT)
!          - DONNE LA SEQUENCE DES POINTS DU POLYGONE (CHAINE),
!            TERMINEE PAR 0
!          - FOURNIT LE NUMERO LIBRE SUIVANT QUAND CELUI-CI EST UTILISE
    integer(kind=8) :: inext(nb*2)
!
!       PREMIER INDICE LIBRE
!       (L'INDICE LIBRE SUIVANT IFREE <- INEXT(IFREE))
    integer(kind=8) :: ifree
!
!       INDICES DE SEGMENT COURANT
!         I...1 = INDICE DE DEBUT DU SEGMENT
!         I...2 = INDICE DE FIN DE SEGMENT
    integer(kind=8) :: icur1, icur2
!
!       INDICES DE NOUVEAU SEGMENT
!         I...1 = INDICE DE DEBUT DU SEGMENT
!         I...2 = INDICE DE FIN DE SEGMENT
    integer(kind=8) :: inew1, inew2
!
!       CARACTERISTIQUE DES SEGMENTS (CALCUL DE L'INTERSECTION)
    real(kind=8) :: phicur
    real(kind=8) :: phinew
!
!       RECHERCHE DU MIN
    real(kind=8) :: tmp0
    real(kind=8) :: tmp1
!
    integer(kind=8) :: i, ii, j
!
!       INITIALISATION DES SUIVANTS
    plgnnb = nb*2
    do 30 icur1 = 1, plgnnb
        inext(icur1) = icur1+1
        xp(icur1) = 1d+6
        yp(icur1) = 1d+6
30      continue
!
!       CREATION DE 2 PREMIERS SEGMENTS A INTERSECTER
        i = nb/2
!
        icur1 = 1
        tmp0 = atab(i+1)
        xp(icur1) = tmp0
        yp(icur1) = 1d+6
        ap(icur1) = 1
!
        icur1 = 2
        tmp0 = atab(i+1)
        xp(icur1) = tmp0
        tmp0 = atab(1)
        yp(icur1) = tmp0
        ap(icur1) = i+1
!
        icur1 = 3
        xp(icur1) = 1d+6
        tmp0 = atab(1)
        yp(icur1) = tmp0
        ap(icur1) = 1
!
!       FIN (PROVISOIRE) DE LA CHAINE DES INDICES DE POINTS DU POLYGONE
        inext(3) = 0
!       PREMIER INDICE LIBRE
        ifree = 4
!
        do 40 i = 2, nb
            if (i .ne. (nb/2)+1) then
                phicur = 1d0
                icur1 = 1
                inew2 = 0
                do 400 j = 1, plgnnb
                    icur2 = inext(icur1)
                    phinew = fcttab(i, 1)*xp(icur2)+fcttab(i, 2)*yp(icur2)-atab(i)
                    if (phinew .lt. -precis) then
                        if (phicur .lt. -precis) then
                            inext(icur1) = ifree
                            ifree = icur1
                        else
                            inew1 = ifree
                            ifree = inext(ifree)
                            ii = ap(icur2)
                            tmp0 = 1d0/(fcttab(i, 1)*fcttab(ii, 2)-fcttab(ii, 1)*fcttab(i, 2))
                            xp(inew1) = (atab(i)*fcttab(ii, 2)-atab(ii)*fcttab(i, 2))*tmp0
                            yp(inew1) = (atab(ii)*fcttab(i, 1)-atab(i)*fcttab(ii, 1))*tmp0
                            ap(inew1) = ii
                            inext(icur1) = inew1
                        end if
                    else
                        if (phicur .lt. -precis) then
                            inew2 = ifree
                            ifree = inext(ifree)
                            ii = ap(icur2)
                            tmp0 = 1d0/(fcttab(i, 1)*fcttab(ii, 2)-fcttab(ii, 1)*fcttab(i, 2))
                            xp(inew2) = (atab(i)*fcttab(ii, 2)-atab(ii)*fcttab(i, 2))*tmp0
                            yp(inew2) = (atab(ii)*fcttab(i, 1)-atab(i)*fcttab(ii, 1))*tmp0
                            ap(inew2) = i
                            inext(inew1) = inew2
                            inext(inew2) = icur2
                            inext(icur1) = ifree
                            ifree = icur1
                        end if
                    end if
                    if ((inext(icur2) .lt. 1) .or. (inew2 .ge. 1)) then
                        goto 401
                    end if
                    icur1 = icur2
                    phicur = phinew
400                 continue
401                 continue
                    end if
40                  continue
!
!       RECHERCHE DU MINIMUM DE XP(*)+YP(*)
                    icur1 = 1
                    tmp0 = 1d99
50                  continue
                    tmp1 = xp(icur1)+yp(icur1)
!         -- LORSQUE 2 FACETTES DONNENT QUASIMENT LE MEME XP+YP,
!            ON VEUT QUE L'ALGORITHME TROUVE TOUJOURS LA MEME FACETTE
!            QUELQUE SOIENT LES OPTIONS DE COMPILATION.
!            POUR CELA, ON FAVORISE LA 1ERE FACETTE RENCONTREE EN
!            NE CHANGEANT DE MINIMUM QUE SI LE GAIN EST APRECIABLE :
                    if (tmp1 .le. 0.9999d0*tmp0) then
                        tmp0 = tmp1
                        icur2 = icur1
                    end if
                    icur1 = inext(icur1)
                    if (icur1 .gt. 0) then
                        goto 50
                    end if
!
                    ax = xp(icur2)
                    ay = yp(icur2)
!
!
                    end subroutine
