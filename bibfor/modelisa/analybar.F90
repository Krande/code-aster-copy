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
subroutine analybar(x3d1, x3d2, x3d3, x3dp, xbar, &
                    excent, iproj, inoeu, icote)
    implicit none
!  DESCRIPTION :
!  -----------
!       APRES ECHEC DE TSTBAR SUR UNE FACE, ANALYSE A PARTIR D'UN
!       EXCENTREMENT DONNE SI LE POINT PROJETE EST SUFFISAMMENT PRES DE
!       LA FACE POUR POUVOIR TENTER UNE PROJECTION SUR COTE OU SUR NOEUD
!       DANS LE CADRE DE LA LIAISON CABLE/COQUE DE DEFI_CABLE_BP
!
!       LA VALEUR DE ALPHA_MAX EST FIXEE A 45°, CELA EXPRIME LE FAIT
!       QUE L'ON CONSIDERE QUE DEUX MAILLES ADJACENTES NE DOIVENT PAS
!       AVOIR UN ANGLE DE PLUS DE 45° ENTRE ELLES (LE PLAT ETANT 0°)
!
!       OUT : IPROJ = -1 si noeud trop loin de la maille
!                   = 20 si projection sur côté possible
!                   = 30 si projection sur noeud possible
!       OUT : INOEU = 1, 2 ou 3 (numéro du noeud)
!       OUT : ICOTE = 1, 2 ou 3 (numéro du côté)
!-------------------   DECLARATION DES VARIABLES   ---------------------
!
!
#include "asterc/r8pi.h"
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/calc_h_tria.h"
#include "asterfort/tstbar.h"
! ARGUMENTS
! ---------
    integer(kind=8) :: iproj, inoeu, icote
    real(kind=8) :: xbar(3), excent, x3d1(3), x3d2(3), x3d3(3), x3dp(3)
!
! VARIABLES LOCALES
! -----------------
    real(kind=8) :: h, d, alpha, alpha_max, r8bid3(3), xbar2(2)
    integer(kind=8) :: j, ino, jno, kno, lno_neg(2), iproj2
!
!
!-------------------   DEBUT DU CODE EXECUTABLE    ---------------------
!
!
    alpha_max = r8pi()/4.d0
    if (excent .lt. r8prem()) then
        iproj = -1
        goto 77
    end if
!
    j = 0
    do ino = 1, 3
        if (xbar(ino) .lt. 0.d0) then
!           calcul de la hauteur
            call calc_h_tria(ino, x3d1, x3d2, x3d3, h)
!           distance au côté opposé au noeud
            d = -h*xbar(ino)
!           calcul de l'angle
            alpha = atan(d/excent)
            if (alpha .le. alpha_max) then
                j = j+1
                lno_neg(j) = ino
            else
!           pas de projection segment ou noeud possible
                iproj = -1
                goto 77
            end if
        end if
    end do
    if (j .eq. 2) then
!           projection possible sur le noeud pas dans la liste lno_neg
        iproj = 30
        do ino = 1, 3
            do jno = 1, 2
                if (ino .eq. lno_neg(jno)) goto 20
            end do
            exit
20          continue
        end do
        inoeu = ino
    else if (j .eq. 1) then
!           projection possible sur le côté opposé au noeud ou sur un
!           des deux autres noeuds
        ino = lno_neg(1)
        jno = ino+1
        if (jno .gt. 3) jno = 1
        kno = jno+1
        if (kno .gt. 3) kno = 1
        if (jno .eq. 1) then
            call tstbar(2, x3d1, x3d2, r8bid3, r8bid3, &
                        x3dp, xbar2, iproj2)
        else if (jno .eq. 2) then
            call tstbar(2, x3d2, x3d3, r8bid3, r8bid3, &
                        x3dp, xbar2, iproj2)
        else
            call tstbar(2, x3d3, x3d1, r8bid3, r8bid3, &
                        x3dp, xbar2, iproj2)
        end if
        if (iproj2 .eq. 0) then
!               projection possible sur un côté
            iproj = 20
            icote = jno
        else
            iproj = 30
            if (xbar2(1) .le. 0.01d0) then
                inoeu = kno
            else if (xbar2(2) .le. 0.01d0) then
                inoeu = jno
            else
                ASSERT(.false.)
            end if
        end if
    end if
77  continue
!
end subroutine
