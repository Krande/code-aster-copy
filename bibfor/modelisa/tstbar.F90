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
subroutine tstbar(nbsom, x3d1, x3d2, x3d3, x3d4, &
                  x3dp, xbar, itest)
    implicit none
!  DESCRIPTION : TEST D'APPARTENANCE DU POINT X3DP(3) AU DOMAINE
!  -----------   GEOMETRIQUE DEFINI
!                  - PAR LE SEGMENT  DE SOMMETS  X3D1(3), X3D2(3)
!                    SI NBSOM = 2
!                  - PAR LE TRIANGLE DE SOMMETS  X3D1(3), X3D2(3),
!                    X3D3(3) SI NBSOM = 3
!                  - PAR LE TETRAEDRE DE SOMMETS X3D1(3), X3D2(3),
!                    X3D3(3), X3D4(3) SI NBSOM = 4
!
!                CALCUL DES COORDONNEES BARYCENTRIQUES PUIS TEST
!
!                APPELANT : IMMEHX, IMMEPN, IMMEPY, IMMETT,
!                           PROJSG, PROJTQ
!
!  IN     : NBSOM  : INTEGER , SCALAIRE
!                    NOMBRE DE SOMMETS A CONSIDERER
!                    NBSOM=2 : TEST D'APPARTENANCE AU SEGMENT   1-2
!                    NBSOM=3 : TEST D'APPARTENANCE AU TRIANGLE  1-2-3
!                    NBSOM=4 : TEST D'APPARTENANCE AU TETRAEDRE 1-2-3-4
!  IN     : X3D1   : REAL*8 , VECTEUR DE DIMENSION 3
!                    COORDONNEES DU PREMIER SOMMET
!  IN     : X3D2   : REAL*8 , VECTEUR DE DIMENSION 3
!                    COORDONNEES DU DEUXIEME SOMMET
!  IN     : X3D3   : REAL*8 , VECTEUR DE DIMENSION 3
!                    COORDONNEES DU TROISIEME SOMMET
!  IN     : X3D4   : REAL*8 , VECTEUR DE DIMENSION 3
!                    COORDONNEES DU QUATRIEME SOMMET
!  IN     : X3DP   : REAL*8 , VECTEUR DE DIMENSION 3
!                    COORDONNEES DU POINT CONSIDERE
!
!  OUT    : XBAR   : REAL*8 , VECTEUR DE DIMENSION 2 OU 3 OU 4
!                    COORDONNEES BARYCENTRIQUES DU POINT CONSIDERE
!                    POUR UN SEGMENT   : BARYCENTRE DES SOMMETS 1-2
!                    POUR UN TRIANGLE  : BARYCENTRE DES SOMMETS 1-2-3
!                    POUR UN TETRAEDRE : BARYCENTRE DES SOMMETS 1-2-3-4
!
!  OUT    : ITEST  : INTEGER , SCALAIRE
!                    INDICATEUR DE RESULTAT DU TEST
!
!                    ITEST =  -1  LE POINT CONSIDERE N'APPARTIENT PAS
!                                 AU DOMAINE
!
!                    ITEST =   0  LE POINT CONSIDERE EST A L'INTERIEUR
!                                 DU DOMAINE
!
!                    ITEST =   2  LE POINT CONSIDERE COINCIDE AVEC
!                                 UN DES NOEUDS SOMMETS
!
!                    POUR UN DOMAINE TRIANGLE : SI LE POINT CONSIDERE
!                    SE TROUVE SUR UNE ARETE
!                    ITEST =  10 + NUMERO D'ARETE
!
!                    POUR UN DOMAINE TETRAEDRE : SI LE POINT CONSIDERE
!                    SE TROUVE SUR UNE FACE
!                    ITEST = 100 + 10 * NUMERO DE FACE
!                    SI LE POINT CONSIDERE SE TROUVE SUR UNE ARETE
!                    ITEST = 100 + 10 * NUMERO DE FACE + NUMERO D'ARETE
!
!-------------------   DECLARATION DES VARIABLES   ---------------------
!
! ARGUMENTS
! ---------
#include "asterc/r8prem.h"
#include "asterfort/mtcrog.h"
#include "asterfort/rslsvd.h"
#include "blas/dcopy.h"
    integer(kind=8) :: nbsom, itest
    real(kind=8) :: x3d1(*), x3d2(*), x3d3(*), x3d4(*), x3dp(*), xbar(*)
!
! VARIABLES LOCALES
! -----------------
    integer(kind=8) :: ierr
    real(kind=8) :: a(4, 4), b(4), eps, epsg, s(3), u(4, 4), v(4, 3), work(4)
    blas_int :: b_incx, b_incy, b_n
!
!
!-------------------   DEBUT DU CODE EXECUTABLE    ---------------------
!
    itest = -1
!
    eps = 1.0d+02*r8prem()
    epsg = 1.0d-5
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 1   CALCUL DES COORDONNEES BARYCENTRIQUES
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! 1.1 CONSTRUCTION DE LA MATRICE A ET DU SECOND MEMBRE
! --- DU SYSTEME A RESOUDRE
!
    a(1, 1) = x3d1(1)-x3dp(1)
    a(2, 1) = x3d1(2)-x3dp(2)
    a(3, 1) = x3d1(3)-x3dp(3)
    a(4, 1) = 1.0d0
    a(1, 2) = x3d2(1)-x3dp(1)
    a(2, 2) = x3d2(2)-x3dp(2)
    a(3, 2) = x3d2(3)-x3dp(3)
    a(4, 2) = 1.0d0
    if (nbsom .gt. 2) then
        a(1, 3) = x3d3(1)-x3dp(1)
        a(2, 3) = x3d3(2)-x3dp(2)
        a(3, 3) = x3d3(3)-x3dp(3)
        a(4, 3) = 1.0d0
        if (nbsom .eq. 4) then
            a(1, 4) = x3d4(1)-x3dp(1)
            a(2, 4) = x3d4(2)-x3dp(2)
            a(3, 4) = x3d4(3)-x3dp(3)
            a(4, 4) = 1.0d0
        end if
    end if
!
    b(1) = 0.0d0
    b(2) = 0.0d0
    b(3) = 0.0d0
    b(4) = 1.0d0
!
! 1.2 RESOLUTION DU SYSTEME
! ---
    if (nbsom .lt. 4) then
!
! 1.2.1  SI NBSOM = 2 OU 3 : SYSTEME SUR-CONTRAINT
! .....  RESOLUTION AU SENS DES MOINDRES CARRES PAR DECOMPOSITION AUX
!        VALEURS SINGULIERES DE LA MATRICE A
!
        call rslsvd(4, 4, nbsom, a(1, 1), s(1), &
                    u(1, 1), v(1, 1), 1, b(1), eps, &
                    ierr, work(1))
        if (ierr .ne. 0) goto 9999
        b_n = to_blas_int(nbsom)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, b(1), b_incx, xbar(1), b_incy)
!
    else
!
! 1.2.2  SI NBSOM = 4 : SYSTEME EQUI-CONTRAINT
! .....  RESOLUTION EXACTE
!
        call mtcrog(a(1, 1), b(1), 4, 4, 1, &
                    xbar(1), work(1), ierr)
        if (ierr .ne. 0) goto 9999
!
    end if
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 2   DISCUSSION DE L'APPARTENANCE DU POINT CONSIDERE AU DOMAINE :
!     LES COORDONNEES BARYCENTRIQUES DOIVENT ETRE POSITIVES
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    if ((xbar(1) .lt. 0.0d0) .and. (dble(abs(xbar(1))) .gt. epsg)) goto 9999
    if ((xbar(2) .lt. 0.0d0) .and. (dble(abs(xbar(2))) .gt. epsg)) goto 9999
    if (nbsom .gt. 2) then
        if ((xbar(3) .lt. 0.0d0) .and. (dble(abs(xbar(3))) .gt. epsg)) goto 9999
        if (nbsom .eq. 4) then
            if ((xbar(4) .lt. 0.0d0) .and. (dble(abs(xbar(4))) .gt. epsg)) goto 9999
        end if
    end if
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 3   DISCUSSION PRECISE DE LA SITUATION DU POINT CONSIDERE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    itest = 0
!
! 3.1 DOMAINE SEGMENT
! ---
    if (nbsom .eq. 2) then
!
        if ((dble(abs(xbar(1))) .lt. epsg) .or. (dble(abs(xbar(2))) .lt. epsg)) itest = &
            2
!
! 3.2 DOMAINE TRIANGLE
! ---
    else if (nbsom .eq. 3) then
!
! 3.2.1  POINT SUR L'ARETE ( SOMMET 1 , SOMMET 2 )
! .....
        if (dble(abs(xbar(3))) .lt. epsg) then
            itest = 11
!.......... COINCIDENCE AVEC LE SOMMET 1 OU LE SOMMET 2
            if ((dble(abs(xbar(2))) .lt. epsg) .or. (dble(abs(xbar(1))) .lt. epsg)) itest = &
                2
!
! 3.2.2  POINT SUR L'ARETE ( SOMMET 2 , SOMMET 3 )
! .....
        else if (dble(abs(xbar(1))) .lt. epsg) then
            itest = 12
!.......... COINCIDENCE AVEC LE SOMMET 3
            if (dble(abs(xbar(2))) .lt. epsg) itest = 2
!
! 3.2.3  POINT SUR L'ARETE ( SOMMET 3 , SOMMET 1 )
! .....
        else if (dble(abs(xbar(2))) .lt. epsg) then
            itest = 13
        end if
!
! 3.3 DOMAINE TETRAEDRE
! ---
    else
!
! 3.3.1  POINT SUR LA FACE ( SOMMET 1 , SOMMET 2 , SOMMET 3 )
! .....
        if (dble(abs(xbar(4))) .lt. epsg) then
            itest = 110
!.......... POINT SUR L'ARETE ( SOMMET 1 , SOMMET 2 )
            if (dble(abs(xbar(3))) .lt. epsg) then
                itest = 111
!............. COINCIDENCE AVEC LE SOMMET 1 OU LE SOMMET 2
                if ((dble(abs(xbar(2))) .lt. epsg) .or. (dble(abs(xbar(1))) .lt. epsg)) &
                    itest = 2
!.......... POINT SUR L'ARETE ( SOMMET 2 , SOMMET 3 )
            else if (dble(abs(xbar(1))) .lt. epsg) then
                itest = 112
!............. COINCIDENCE AVEC LE SOMMET 3
                if (dble(abs(xbar(2))) .lt. epsg) itest = 2
!.......... POINT SUR L'ARETE ( SOMMET 3 , SOMMET 1 )
            else if (dble(abs(xbar(2))) .lt. epsg) then
                itest = 113
            end if
!
! 3.3.2  POINT SUR LA FACE ( SOMMET 1 , SOMMET 2 , SOMMET 4 )
! .....
        else if (dble(abs(xbar(3))) .lt. epsg) then
            itest = 120
!.......... POINT SUR L'ARETE ( SOMMET 2 , SOMMET 4 )
            if (dble(abs(xbar(1))) .lt. epsg) then
                itest = 122
!............. COINCIDENCE AVEC LE SOMMET 4
                if (dble(abs(xbar(2))) .lt. epsg) itest = 2
!.......... POINT SUR L'ARETE ( SOMMET 4 , SOMMET 1 )
            else if (dble(abs(xbar(2))) .lt. epsg) then
                itest = 123
            end if
!
! 3.3.3  POINT SUR LA FACE ( SOMMET 2 , SOMMET 3 , SOMMET 4 )
! .....
        else if (dble(abs(xbar(1))) .lt. epsg) then
            itest = 130
!.......... POINT SUR L'ARETE ( SOMMET 3 , SOMMET 4 )
            if (dble(abs(xbar(2))) .lt. epsg) itest = 132
!
! 3.3.4  POINT SUR LA FACE ( SOMMET 3 , SOMMET 1 , SOMMET 4 )
! .....
        else if (dble(abs(xbar(2))) .lt. epsg) then
            itest = 140
        end if
!
    end if
!
9999 continue
!
! --- FIN DE TSTBAR.
end subroutine
