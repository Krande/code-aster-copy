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
subroutine projtq(nbcnx, xyzma, icnx, x3dp, excent, &
                  itria, inoeu, icote, xbar, iproj)
    implicit none
!  DESCRIPTION : TEST D'APPARTENANCE DU POINT PROJETE X3DP(3)
!  -----------   AU DOMAINE GEOMETRIQUE DEFINI PAR UNE MAILLE
!                TRIANGLE OU QUADRANGLE
!
!                POUR CELA, ON CALCULE LES COORDONNEES BARYCENTRIQUES
!                DU POINT PROJETE
!                  - DANS LE TRIANGLE 1-2-3 POUR UNE MAILLE TRIANGLE
!                  - DANS LE TRIANGLE 1-2-3 PUIS SI NECESSAIRE DANS
!                    LE TRIANGLE 3-4-1 POUR UNE MAILLE QUADRANGLE
!                CETTE METHODE EST EXACTE POUR DES MAILLES A BORDS
!                DROITS SUPPORTANT DES ELEMENTS LINEAIRES, APPROCHEE
!                POUR DES MAILLES A BORDS COURBES SUPPORTANT DES
!                ELEMENTS QUADRATIQUES OU AUTRES
!
!                APPELANT : PROJKM
!
!  IN     : NBCNX  : INTEGER , SCALAIRE
!                    NOMBRE DE NOEUDS DE LA MAILLE
!  IN     : XYZMA  : REAL*8 , TABLEAU DE DIMENSIONS (3,NNOMAX)
!                    CONTIENT LES COORDONNEES DES NOEUDS DE LA MAILLE
!  IN     : ICNX   : INTEGER , SCALAIRE
!                    RANG DU NOEUD BETON LE PLUS PROCHE DU POINT PROJETE
!                    DANS LA TABLE DE CONNECTIVITE DE LA MAILLE
!  IN     : X3DP   : REAL*8 , VECTEUR DE DIMENSION 3
!                    COORDONNEES DU POINT PROJETE
!  IN     : EXCENT : REAL*8 (POSITIF), SCALAIRE
!                    EXCENTREMENT (VAL ABS) DU POINT PAR RAPPORT A LA MAILLE
!  OUT    : ITRIA  : INTEGER , SCALAIRE
!                    SI PROJECTION REUSSIE : INDICATEUR DU SOUS-DOMAINE
!                    AUQUEL APPARTIENT LE POINT PROJETE X3DP(3) :
!                    ITRIA = 1 : TRIANGLE 1-2-3
!                    ITRIA = 2 : TRIANGLE 3-4-1
!  OUT    : INOEU  : INTEGER , SCALAIRE
!                    SI PROJECTION ECHOUE MAIS POSSIBLE PROJECTION
!                    SUR UN NOEUD, NUMERO DU NOEUD (1, 2, 3 ou 4)
!  OUT    : ICOTE  : INTEGER , SCALAIRE
!                    SI PROJECTION ECHOUE MAIS POSSIBLE PROJECTION
!                    SUR UN COTE, NUMERO DU COTE (1, 2, 3 ou 4)
!  OUT    : XBAR   : REAL*8 , VECTEUR DE DIMENSION 3
!                    SI PROJECTION REUSSIE : COORDONNEES BARYCENTRIQUES
!                    DU POINT PROJETE (BARYCENTRE DES SOMMETS DU
!                    TRIANGLE 1-2-3 OU 3-4-1)
!  OUT    : IPROJ  : INTEGER , SCALAIRE
!                    INDICE DE PROJECTION
!                    IPROJ = -1  PROJECTION NON REUSSIE
!                    IPROJ =  0  LE POINT PROJETE EST A L'INTERIEUR
!                                DE LA MAILLE
!                    IPROJ =  1X LE POINT PROJETE EST SUR UNE FRONTIERE
!                                DE LA MAILLE
!                    IPROJ =  2  LE POINT PROJETE COINCIDE AVEC UN DES
!                                NOEUDS DE LA MAILLE
!
!-------------------   DECLARATION DES VARIABLES   ---------------------
!
! ARGUMENTS
! ---------
#include "asterf_types.h"
#include "asterc/matfpe.h"
#include "asterc/r8prem.h"
#include "asterfort/analybar.h"
#include "asterfort/assert.h"
#include "asterfort/tstbar.h"
#include "blas/dnrm2.h"
    integer(kind=8) :: nbcnx, icnx, itria, inoeu, icote, iproj
    real(kind=8) :: xyzma(3, *), x3dp(*), excent, xbar(*)
!
! VARIABLES LOCALES
! -----------------
    integer(kind=8) :: ino, nbsom
    real(kind=8) :: d, dx, dy, dz, epsg, nrm2, r8bid3(3)
    aster_logical :: notlin
    blas_int :: b_incx, b_n
!
!
!-------------------   DEBUT DU CODE EXECUTABLE    ---------------------
!
    call matfpe(-1)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 1   INITIALISATIONS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    epsg = 1.0d+08*r8prem()
    r8bid3(1) = 0.d0
    r8bid3(2) = 0.d0
    r8bid3(3) = 0.d0
!
    notlin = (nbcnx .gt. 4)
    if ((nbcnx .eq. 3) .or. (nbcnx .eq. 6)) then
        nbsom = 3
    else
        nbsom = 4
    end if
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 2   TEST POUR UNE MAILLE TRIANGLE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    if (nbsom .eq. 3) then
!
        itria = 1
!
!....... TEST D'APPARTENANCE AU TRIANGLE 1-2-3, PAR DETERMINATION DES
!....... COORDONNEES BARYCENTRIQUES
!
        call tstbar(3, xyzma(1, 1), xyzma(1, 2), xyzma(1, 3), r8bid3, &
                    x3dp, xbar, iproj)
        if (iproj .lt. 0) then
            call analybar(xyzma(1, 1), xyzma(1, 2), xyzma(1, 3), x3dp, xbar, &
                          excent, iproj, inoeu, icote)
        end if
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 3   TESTS POUR UNE MAILLE QUADRANGLE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! 3.1 APPARTENANCE PLUS PROBABLE AU TRIANGLE 3-4-1
! ---
    else if ((icnx .eq. 4) .or. (icnx .eq. 7) .or. (icnx .eq. 8)) then
!
        itria = 2
!
!....... TEST D'APPARTENANCE AU TRIANGLE 3-4-1, PAR DETERMINATION DES
!....... COORDONNEES BARYCENTRIQUES
!
        call tstbar(3, xyzma(1, 3), xyzma(1, 4), xyzma(1, 1), r8bid3, &
                    x3dp(1), xbar(1), iproj)
!
!....... REAJUSTEMENT DE IPROJ SI APPARTENANCE A UN BORD
!
        if (iproj .eq. 11) then
            iproj = 13
        else if (iproj .eq. 12) then
            iproj = 14
        else if (iproj .eq. 13) then
            iproj = 0
        else if (iproj .lt. 0) then
            call analybar(xyzma(1, 3), xyzma(1, 4), xyzma(1, 1), x3dp, xbar, &
                          excent, iproj, inoeu, icote)
            if (iproj .lt. 0) then
                if (.not. ( &
                    xbar(2) .lt. 0.d0 .and. xbar(1) .ge. 0.d0 .and. xbar(3) .ge. 0.d0)) then
                    goto 999
                end if
            else if (iproj .eq. 20 .and. icote .eq. 3) then
                iproj = -1
            else if (iproj .eq. 20) then
                if (icote .eq. 1) then
                    icote = 3
                else if (icote .eq. 2) then
                    icote = 4
                else
!                   on ne doit pas passer par la
                    ASSERT(.false.)
                end if
            else if (iproj .eq. 30) then
                if (inoeu .eq. 1) then
                    inoeu = 3
                else if (inoeu .eq. 2) then
                    inoeu = 4
                else if (inoeu .eq. 3) then
                    inoeu = 1
                else
                    ASSERT(.false.)
                end if
            else
                ASSERT(.false.)
            end if
        end if
!
!....... EN CAS D'ECHEC TEST D'APPARTENANCE AU TRIANGLE 1-2-3,
!....... PAR DETERMINATION DES COORDONNEES BARYCENTRIQUES
!
        if (iproj .lt. 0) then
            itria = 1
            call tstbar(3, xyzma(1, 1), xyzma(1, 2), xyzma(1, 3), r8bid3, &
                        x3dp(1), xbar(1), iproj)
!.......... REAJUSTEMENT DE IPROJ SI PROJECTION SUR LE TROISIEME COTE
!.......... DU TRIANGLE 1-2-3
            if (iproj .eq. 13) then
!           on ne doit normalement pas passer par la !
                iproj = 0
            else if (iproj .lt. 0) then
                call analybar(xyzma(1, 1), xyzma(1, 2), xyzma(1, 3), x3dp, xbar, &
                              excent, iproj, inoeu, icote)
                if (iproj .lt. 0) then
                    goto 999
                else if (iproj .eq. 20) then
                    ASSERT(icote .eq. 1 .or. icote .eq. 2)
                else if (iproj .eq. 30) then
                    ASSERT(inoeu .eq. 2)
                else
                    ASSERT(.false.)
                end if
            end if
        end if
!
! 3.2 APPARTENANCE PLUS PROBABLE AU TRIANGLE 1-2-3
! ---
    else
!
        itria = 1
!
!....... TEST D'APPARTENANCE AU TRIANGLE 1-2-3, PAR DETERMINATION DES
!....... COORDONNEES BARYCENTRIQUES
!
        call tstbar(3, xyzma(1, 1), xyzma(1, 2), xyzma(1, 3), r8bid3, &
                    x3dp(1), xbar(1), iproj)
!
!....... REAJUSTEMENT DE IPROJ SI PROJECTION SUR LE TROISIEME COTE
!....... DU TRIANGLE 1-2-3
!
        if (iproj .eq. 13) then
            iproj = 0
        else if (iproj .lt. 0) then
            call analybar(xyzma(1, 1), xyzma(1, 2), xyzma(1, 3), x3dp, xbar, &
                          excent, iproj, inoeu, icote)
            if (iproj .lt. 0) then
                if (.not. ( &
                    xbar(2) .lt. 0.d0 .and. xbar(1) .ge. 0.d0 .and. xbar(3) .ge. 0.d0)) then
                    goto 999
                end if
            else if (iproj .eq. 20 .and. icote .eq. 3) then
                iproj = -1
            end if
        end if
!
!....... EN CAS D'ECHEC TEST D'APPARTENANCE AU TRIANGLE 3-4-1,
!....... PAR DETERMINATION DES COORDONNEES BARYCENTRIQUES
!
        if (iproj .lt. 0) then
            itria = 2
            call tstbar(3, xyzma(1, 3), xyzma(1, 4), xyzma(1, 1), r8bid3, &
                        x3dp(1), xbar(1), iproj)
!.......... REAJUSTEMENT DE IPROJ SI APPARTENANCE A UN BORD
            if (iproj .eq. 11) then
                iproj = 13
            else if (iproj .eq. 12) then
                iproj = 14
            else if (iproj .eq. 13) then
                iproj = 0
            else if (iproj .lt. 0) then
                call analybar(xyzma(1, 3), xyzma(1, 4), xyzma(1, 1), x3dp, xbar, &
                              excent, iproj, inoeu, icote)
                if (iproj .lt. 0) then
                    goto 999
                else if (iproj .eq. 20) then
                    if (icote .eq. 1) then
                        icote = 3
                    else if (icote .eq. 2) then
                        icote = 4
                    else
                        ASSERT(.false.)
                    end if
                else if (iproj .eq. 30) then
                    ASSERT(inoeu .eq. 2)
                    inoeu = 4
                else
                    ASSERT(.false.)
                end if
            end if
        end if
!
    end if
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 4   TESTS COMPLEMENTAIRES POUR LES MAILLES A BORDS COURBES SUPPORTANT
!     DES ELEMENTS NON LINEAIRES
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    if (notlin) then
!
! 4.1    TEST DE COINCIDENCE AVEC UN NOEUD MILIEU
! ---
        if (iproj .gt. 10) then
            ino = iproj-10+nbsom
            b_n = to_blas_int(3)
            b_incx = to_blas_int(1)
            nrm2 = dnrm2(b_n, xyzma(1, ino), b_incx)
            if (nrm2 .eq. 0.0d0) nrm2 = 1.0d0
            dx = xyzma(1, ino)-x3dp(1)
            dy = xyzma(2, ino)-x3dp(2)
            dz = xyzma(3, ino)-x3dp(3)
            d = dble(sqrt(dx*dx+dy*dy+dz*dz))
            if (d/nrm2 .lt. epsg) iproj = 2
        end if
!
! 4.2    TEST DE COINCIDENCE AVEC LE NOEUD CENTRAL POUR UNE MAILLE QUAD9
! ---
        if ((iproj .eq. 0) .and. (nbcnx .eq. 9)) then
            b_n = to_blas_int(3)
            b_incx = to_blas_int(1)
            nrm2 = dnrm2(b_n, xyzma(1, 9), b_incx)
            if (nrm2 .eq. 0.0d0) nrm2 = 1.0d0
            dx = xyzma(1, 9)-x3dp(1)
            dy = xyzma(2, 9)-x3dp(2)
            dz = xyzma(3, 9)-x3dp(3)
            d = dble(sqrt(dx*dx+dy*dy+dz*dz))
            if (d/nrm2 .lt. epsg) iproj = 2
        end if
!
    end if
!
! --- FIN DE PROJTQ.
999 continue
    call matfpe(1)
!
end subroutine
