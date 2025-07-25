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
subroutine norlin(typma, l, knumai, coor, dfonc, &
                  in, prec, a, b, c)
!
    implicit none
#include "asterfort/assert.h"
#include "asterfort/utmess.h"
    character(len=3) :: typma
    character(len=8) :: knumai
    real(kind=8) :: coor(3, *), dfonc(*), a, b, c, prec
    integer(kind=8) :: in, l
!
!     OBJECTIF : CALCULER LA NORMALE OU LA TANGENTE A UNE MAILLE
!                QUADRATIQUE EN LINEAIRE ( SEULS LES NOEUDS SOMMETS
!                IMPORTENT DONC ) SUR UN NOEUD
!
!     ARGUMENT D ENTREE
!     TYPMA  : TYPE DE LA MAILLE
!     L      : VAUT 2 SI ON CHERCHE LA TANGENTE ET 1 POUR LA NORMALE
!     KNUMAI : NUMERO DE LA MAILLE ACTUELLE ( UTILE SEULEMENT POUR LES
!              MESSAGES D ALARME/D ERREUR )
!     COOR   : COORDONNEES DES NOEUDS DE LA MAILLE
!     DFONC  : VALEURS DES DERIVES DES FONCTIONS DE FORMES AUX NOEUDS
!     IN     : NUMERO DU NOEUD DANS LA CONNECTIVITE DE LA MAILLE
!     PREC   : PRECISION SUR LA NORME DE LA NORMALE AU DESSOUS DE
!              LAQUELLE ON CONSIDERE QUE LA MAILLE EST DEGENEREE
!     A      : VALEUR SELON X DU VECTEUR NORMAL
!     B      : VALEUR SELON Y DU VECTEUR NORMAL
!     C      : VALEUR SELON Z DU VECTEUR NORMAL
!
    integer(kind=8) :: ifonc, nn
    real(kind=8) :: eks1x, eks1y, eks1z, eet1x, eet1y, eet1z, norme
    real(kind=8) :: eks2x, eks2y, eks2z, eet2x, eet2y, eet2z
!
!
!     ON CALCULE LE NOMBRE DE NOEUDS SOMMETS
    if (typma .eq. 'QU8') then
        nn = 4
    else if (typma .eq. 'TR6') then
        nn = 3
    else if (typma .eq. 'SE3') then
        nn = 2
    else
        ASSERT(.false.)
    end if
!     CAS 3D
    if (typma .ne. 'SE3') then
!     ON EST SUR UN NOEUD SOMMET : ON CALCULE SA NORMALE
        if (in .le. nn) then
            eks1x = 0.d0
            eks1y = 0.d0
            eks1z = 0.d0
            eet1x = 0.d0
            eet1y = 0.d0
            eet1z = 0.d0
            do ifonc = 1, nn
                eks1x = eks1x+coor(1, ifonc)*dfonc((in-1)*nn*2+ifonc)
                eks1y = eks1y+coor(2, ifonc)*dfonc((in-1)*nn*2+ifonc)
                eks1z = eks1z+coor(3, ifonc)*dfonc((in-1)*nn*2+ifonc)
!
                eet1x = eet1x+coor(1, ifonc)*dfonc((in-1)*nn*2+nn+ifonc)
                eet1y = eet1y+coor(2, ifonc)*dfonc((in-1)*nn*2+nn+ifonc)
                eet1z = eet1z+coor(3, ifonc)*dfonc((in-1)*nn*2+nn+ifonc)
            end do
!         CALCUL DU VECTEUR NORMAL ET NORMALISATION
            a = eks1y*eet1z-eks1z*eet1y
            b = eks1z*eet1x-eks1x*eet1z
            c = eks1x*eet1y-eks1y*eet1x
            norme = sqrt(a*a+b*b+c*c)
            if (norme .gt. prec) then
                a = a/norme
                b = b/norme
                c = c/norme
            else
                call utmess('F', 'CHARGES2_26', sk=knumai)
            end if
!     ON EST SUR UN NOEUD MILIEUX : ON CALCULE SA NORMALE EN FAISANT
!     LA MOYENNE DES DEUX NOEUDS DE SON ARETE
        else
            eks1x = 0.d0
            eks1y = 0.d0
            eks1z = 0.d0
            eet1x = 0.d0
            eet1y = 0.d0
            eet1z = 0.d0
            eks2x = 0.d0
            eks2y = 0.d0
            eks2z = 0.d0
            eet2x = 0.d0
            eet2y = 0.d0
            eet2z = 0.d0
            do ifonc = 1, nn
                eks1x = eks1x+coor(1, ifonc)*dfonc((in-1-nn)*nn*2+ifonc)
                eks1y = eks1y+coor(2, ifonc)*dfonc((in-1-nn)*nn*2+ifonc)
                eks1z = eks1z+coor(3, ifonc)*dfonc((in-1-nn)*nn*2+ifonc)
                eks2x = eks2x+coor(1, ifonc)*dfonc((in-1-nn)*nn*2+ifonc)
                eks2y = eks2y+coor(2, ifonc)*dfonc((in-1-nn)*nn*2+ifonc)
                eks2z = eks2z+coor(3, ifonc)*dfonc((in-1-nn)*nn*2+ifonc)
!
                eet1x = eet1x+coor(1, ifonc)*dfonc((in-nn)*nn*2+nn+ifonc)
                eet1y = eet1y+coor(2, ifonc)*dfonc((in-nn)*nn*2+nn+ifonc)
                eet1z = eet1z+coor(3, ifonc)*dfonc((in-nn)*nn*2+nn+ifonc)
                eet2x = eet2x+coor(1, ifonc)*dfonc((in-nn)*nn*2+nn+ifonc)
                eet2y = eet2y+coor(2, ifonc)*dfonc((in-nn)*nn*2+nn+ifonc)
                eet2z = eet2z+coor(3, ifonc)*dfonc((in-nn)*nn*2+nn+ifonc)
            end do
            a = (eks1y*eet1z-eks1z*eet1y+eks2y*eet2z-eks2z*eet2y)/2.d0
            b = (eks1z*eet1x-eks1x*eet1z+eks2z*eet2x-eks2x*eet2z)/2.d0
            c = (eks1x*eet1y-eks1y*eet1x+eks2x*eet2y-eks2y*eet2x)/2.d0
            norme = sqrt(a*a+b*b+c*c)
            if (norme .gt. prec) then
                a = a/norme
                b = b/norme
                c = c/norme
            else
                call utmess('F', 'CHARGES2_26', sk=knumai)
            end if
        end if
!     ON EST EN 2D
    else
!     ON EST SUR UN NOEUD SOMMET : ON CALCULE SA NORMALE
        if (in .le. nn) then
            eks1x = 0.d0
            eks1y = 0.d0
            do ifonc = 1, nn
                eks1x = eks1x+coor(1, ifonc)*dfonc((in-1)*nn+ifonc)
                eks1y = eks1y+coor(2, ifonc)*dfonc((in-1)*nn+ifonc)
            end do
!         ON S INTERESSE AU VECTEUR TANGENT
            if (l .eq. 2) then
                norme = sqrt(eks1x**2+eks1y**2)
                if (norme .gt. prec) then
                    a = eks1x/norme
                    b = eks1y/norme
                else
                    call utmess('F', 'CHARGES2_23', sk=knumai)
                end if
!         ON S INTERESSE AU VECTEUR NORMAL
            else if (l .eq. 1) then
                norme = sqrt(eks1x**2+eks1y**2)
                if (norme .gt. prec) then
                    a = eks1y/norme
                    b = -eks1x/norme
                else
                    call utmess('F', 'CHARGES2_24', sk=knumai)
                end if
            end if
!     ON EST SUR UN NOEUD MILIEUX : ON CALCULE SA NORMALE EN FAISANT
!     LA MOYENNE DES DEUX NOEUDS DE SON ARETE
        else
            eks1x = 0.d0
            eks1y = 0.d0
            eks2x = 0.d0
            eks2y = 0.d0
            do ifonc = 1, nn
                eks1x = eks1x+coor(1, ifonc)*dfonc((in-1-nn)*nn+ifonc)
                eks1y = eks1y+coor(2, ifonc)*dfonc((in-1-nn)*nn+ifonc)
                eks2x = eks2x+coor(1, ifonc)*dfonc((in-nn)*nn+ifonc)
                eks2y = eks2y+coor(2, ifonc)*dfonc((in-nn)*nn+ifonc)
            end do
            eks1x = (eks1x+eks2x)/2
            eks1y = (eks1y+eks2y)/2
!         ON S INTERESSE AU VECTEUR TANGENT
            if (l .eq. 2) then
                norme = sqrt(eks1x**2+eks1y**2)
                if (norme .gt. prec) then
                    a = eks1x/norme
                    b = eks1y/norme
                else
                    call utmess('F', 'CHARGES2_23', sk=knumai)
                end if
!         ON S INTERESSE AU VECTEUR NORMAL
            else if (l .eq. 1) then
                norme = sqrt(eks1x**2+eks1y**2)
                if (norme .gt. prec) then
                    a = eks1y/norme
                    b = -eks1x/norme
                else
                    call utmess('F', 'CHARGES2_24', sk=knumai)
                end if
            end if
        end if
    end if
!
end subroutine
