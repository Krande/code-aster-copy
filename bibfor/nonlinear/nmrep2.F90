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
subroutine nmrep2(n, r, g, gu, rmin, &
                  rmax, rexm, rexp, posopt)
!
! person_in_charge: mickael.abbas at edf.fr
!
    implicit none
#include "asterf_types.h"
#include "asterc/r8gaem.h"
    integer(kind=8) :: n, posopt
    real(kind=8) :: r(*), g(*)
    real(kind=8) :: gu, rmin, rmax, rexm, rexp
!
! ----------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (ALGORITHME - PILOTAGE)
!
! CALCUL DU MINIMUM POUR LA RECHERCHE LINEAIRE PAR INTERPOLATION
! QUADRATIQUE AVEC BORNES
!
! ----------------------------------------------------------------------
!
!
! I/O N      : NOMBRE DE POINTS CONSIDERES  (N >= 2)
! I/O R      : ABSCISSES DES POINTS R(1) < ... R(I) ... < R(N)
! IN  G      : ORDONNEES DES POINTS
! IN  GU     : ORDONNEE A CONVERGENCE (GAIN)
! I/O RMIN   : ABSCISSE MINIMALE
! I/O RMAX   : ABSCISSE MAXIMALE
! I/O REXM   : INTERVALLE INTERDIT AUTOUR DE 0 EN NEGATIF
! I/O REXP   : INTERVALLE INTERDIT AUTOUR DE 0 EN POSITIF
! OUT POSOPT : POSITION DU NOUVEAU R :  R(POSOPT)
!
! ----------------------------------------------------------------------
!
    aster_logical :: gauche, droite
    integer(kind=8) :: i, c, j, pos
    real(kind=8) :: a, b, det, v, d, x, pente, dg, x1, y1, x2, y2, x3, y3
    real(kind=8) :: appuig, appuid, valg, vald, val, valopt, xopt
    real(kind=8) :: diff
    parameter(diff=1.d-8)
!
! ----------------------------------------------------------------------
!
    valopt = r8gaem()
!
!
    do i = 1, n
!
! --- PROPOSITION D'UN NOUVEAU POINT
!
        if (i .eq. 1) then
!
! ---    EXTRAPOLATION LINEAIRE POUR LE POINT LE PLUS A GAUCHE
!
            pente = (g(2)-g(1))/(r(2)-r(1))
            dg = gu/1.2d0-g(1)
            if (pente .le. 0.d0) goto 50
            if (abs(pente) .gt. (abs(dg)/r8gaem())) then
                x = r(1)+dg/pente
            else
                x = (r(1)+r(2))/2
            end if
        else if (i .eq. n) then
!
! ---    EXTRAPOLATION LINEAIRE POUR LE POINT LE PLUS A DROITE
!
            pente = (g(n)-g(n-1))/(r(n)-r(n-1))
            dg = gu/1.2d0-g(n-1)
            if (pente .ge. 0.d0) goto 50
            if (abs(pente) .gt. abs(dg)/r8gaem()) then
                x = r(n-1)+dg/pente
            else
                x = (r(n-1)+r(n))/2
            end if
        else if (i .gt. 1 .and. i .lt. n) then
!
! ---    INTERPOLATION QUADRATIQUE POUR UN POINT INTERMEDIAIRE
!
            j = i-1
            det = -(r(j)-r(i))*(r(i)-r(i+1))*(r(i+1)-r(j))
!
            if (abs(det) .lt. 1.d0/r8gaem()) then
                goto 50
            else
                a = ((g(j)-g(i))*(r(i)-r(i+1))-(g(i)-g(i+1))*(r(j)-r(i)))/det
            end if
!
            if (a .le. 0) then
                goto 50
            else
                b = ((g(i)-g(i+1))*(r(j)**2-r(i)**2)-(g(j)-g(i))*(r(i)**2-r(i+1)**2))/det
                x = -b/(2*a)
            end if
!
! ---    LE MINIMUM N'EST PAS DANS UN VOISINAGE DES 3 POINTS
!
            if (i-2 .ge. 1) then
                j = i-2
                if (x .le. r(j)) goto 50
            end if
            if (i+2 .le. n) then
                if (x .ge. r(i+2)) goto 50
            end if
        end if
!
!
! ----------------------------------------------------------------------
!                   PROJECTION DU POINT SUR LES BORNES
!                 ET EN DEHORS DES POINTS DEJA CALCULES
! ----------------------------------------------------------------------
!
! --- PROJECTION SUR L'INTERVALLE DE RECHERCHE
!
        if (x .lt. rmin) then
            x = rmin
        end if
        if (x .gt. rmax) then
            x = rmax
        end if
        if (x .lt. 0 .and. x .ge. rexm) then
            x = rexm
        end if
        if (x .ge. 0 .and. x .le. rexp) then
            x = rexp
        end if
!
!
!
!      X EST-IL CONFONDU AVEC UN POINT DEJA CALCULE
        c = 0
        do j = 1, n
            if (abs(r(j)-x) .le. diff) c = j
        end do
        if (c .eq. 0) goto 20
!      LES CHOIX VERS LA GAUCHE OU LA DROITE SONT-ILS LICITES
        gauche = r(c)-rmin .gt. diff
        droite = rmax-r(c) .gt. diff
!      POINTS D'APPUI A GAUCHE ET A DROITE
        if (gauche) then
            if (c .eq. 1) then
                appuig = rmin
                valg = g(1)+(appuig-r(1))/(r(2)-r(1))*(g(2)-g(1))
            else
                appuig = max(rmin, r(c-1))
                valg = g(c-1)+(appuig-r(c-1))/(r(c)-r(c-1))*(g(c)-g(c-1))
            end if
        end if
        if (droite) then
            if (c .eq. n) then
                appuid = rmax
                vald = g(n)+(appuid-r(n))/(r(n-1)-r(n))*(g(n-1)-g(n))
            else
                appuid = min(rmax, r(c+1))
                vald = g(c+1)+(appuid-r(c+1))/(r(c)-r(c+1))*(g(c)-g(c+1))
            end if
        end if
!      UNIQUEMENT LE CHOIX A GAUCHE
        if (.not. droite) x = (r(c)+appuig)/2
!      UNIQUEMENT LE CHOIX A DROITE
        if (.not. gauche) x = (r(c)+appuid)/2
!      LES DEUX CHOIX SONT LICITES
        if (gauche .and. droite) then
            if (valg .le. vald) then
                x = (r(c)+appuig)/2
            else
                x = (r(c)+appuid)/2
            end if
        end if
20      continue
!
!
! ----------------------------------------------------------------------
!             APPROXIMATION DE LA VALEUR DE LA FONCTION EN X
! ----------------------------------------------------------------------
!
!      RECHERCHE DE L'INTERVALLE DANS LEQUEL SE TROUVE X
        pos = n+1
        do j = 1, n
            if (x .le. r(j)) then
                pos = j
                goto 40
            end if
        end do
40      continue
!
!      SI DEUX POINTS : INTERPOLATION LINEAIRE
        if (n .eq. 2) then
            val = g(1)+(x-r(1))/(r(2)-r(1))*(g(2)-g(1))
!      EXTRAPOLATION A GAUCHE
        else if (pos .eq. 1) then
            val = g(1)+(x-r(1))/(r(2)-r(1))*(g(2)-g(1))
!      EXTRAPOLATION A DROITE
        else if (pos .eq. n+1) then
            val = g(n)+(x-r(n))/(r(n-1)-r(n))*(g(n-1)-g(n))
!
!      INTERPOLATION QUADRATIQUE
        else
            valg = r8gaem()
            vald = r8gaem()
!
!        INTERPOLATION ENTRE POS-1, POS ET POS+1
            if (pos+1 .le. n) then
                x1 = r(pos-1)
                x2 = r(pos)
                x3 = r(pos+1)
                y1 = g(pos-1)
                y2 = g(pos)
                y3 = g(pos+1)
                d = (x2-x1)*(x3-x2)*(x1-x3)
                a = ((y2-y1)*(x3-x2)-(y3-y2)*(x2-x1))/d
                b = ((x2**2-x1**2)*(y3-y2)-(x3**2-x2**2)*(y2-y1))/d
                v = y1-a*x1**2-b*x1
                valg = a*x**2+b*x+v
            end if
!        INTERPOLATION ENTRE POS-2, POS-1 ET POS
            if (pos-2 .ge. 1) then
                x1 = r(pos-2)
                x2 = r(pos-1)
                x3 = r(pos)
                y1 = g(pos-2)
                y2 = g(pos-1)
                y3 = g(pos)
                d = (x2-x1)*(x3-x2)*(x1-x3)
                a = ((y2-y1)*(x3-x2)-(y3-y2)*(x2-x1))/d
                b = ((x2**2-x1**2)*(y3-y2)-(x3**2-x2**2)*(y2-y1))/d
                v = y1-a*x1**2-b*x1
                vald = a*x**2+b*x+v
            end if
!
!        ON GARDE LA VALEUR DE L'INTERPOLATION LA PLUS BASSE
            val = min(valg, vald)
        end if
!
! ----------------------------------------------------------------------
!               ON GARDE LE MINIMUM DES X CONSTRUITS
! ----------------------------------------------------------------------
!
        if (val .lt. valopt) then
            valopt = val
            xopt = x
            posopt = pos
        end if
!
50      continue
    end do
!
!
! --- INSERTION DU MINIMUM
!
    do i = n, 1, -1
        if (xopt .gt. r(i)) then
            posopt = i+1
            goto 70
        end if
        r(i+1) = r(i)
        g(i+1) = g(i)
    end do
    posopt = 1
70  continue
!
    n = n+1
    r(posopt) = xopt
!
end subroutine
