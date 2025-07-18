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
function mefin3(nbz, nbgrp, imod, icyl, jmod, &
                jcyl, z, f1, f2, f3, &
                g, h)
    implicit none
!
    integer(kind=8) :: nbz, nbgrp, imod, icyl, jmod, jcyl
    real(kind=8) :: z(*), f1(nbz*nbgrp, *), f2(nbz*nbgrp, *), f3(*), g(*), h(*)
!     CALCUL DE L'INTEGRALE SUR (0,L) DE F3(Z)*F1(Z)*F2''(Z)
!     OU F1 EST LA DEFORMEE DU MODE (IMOD) SUR LE CYLINDRE
!     (ICYL) ET F2 CELLE DU MODE (JMOD) SUR LE CYLINDRE (JCYL)
!     OPERATEUR APPELANT : OP0144 , FLUST3, MEFIST, MEFMAT
! ----------------------------------------------------------------------
!     OPTION DE CALCUL   : CALC_FLUI_STRU , CALCUL DES PARAMETRES DE
!     COUPLAGE FLUIDE-STRUCTURE POUR UNE CONFIGURATION DE TYPE "FAISCEAU
!     DE TUBES SOUS ECOULEMENT AXIAL"
! ----------------------------------------------------------------------
! IN  : NBZ    : NOMBRE DE NOEUDS DE LA DISCRETISATION AXIALE
! IN  : NBGRP  : NOMBRE DE GROUPES D EQUIVALENCE
! IN  : IMOD   : NUMERO DU MODE POUR LA FONCTION F1
! IN  : ICYL   : INDICE DU CYLINDRE POUR LA FONCTION F1
! IN  : JMOD   : NUMERO DU MODE POUR LA FONCTION F2
! IN  : JCYL   : INDICE DU GROUPE DE CYLINDRE POUR LA FONCTION F2
! IN  : Z      : COORDONNEES 'Z' DANS LE REPERE AXIAL DES
!                POINTS DISCRETISES, IDENTIQUES POUR TOUS LES CYLINDRES
! IN  : F1     : PREMIERE FONCTION
! IN  : F2     : DEUXIEME FONCTION
! IN  : F3     : TROISIEME FONCTION
! --  : G      : TABLEAU DE TRAVAIL
! --  : H      : TABLEAU DE TRAVAIL
! OUT : MEFIN3 : INTEGRALE CALCULEE
! ----------------------------------------------------------------------
    real(kind=8) :: mefin3
! ----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
    integer(kind=8) :: n, nbz1, nbz2
!-----------------------------------------------------------------------
    nbz1 = nbz*(icyl-1)
    nbz2 = nbz*(jcyl-1)
!
!     ------------------------------------------------
!     CALCUL DE F"(J) -> G
!     F'(J) -> H : MINIMISATION QUADRATIQUE DES RESTES
!     DES DEVELOPPEMENTS DE TAYLOR DE F(J,N)
!     A GAUCHE ET A DROITE
!     H' -> G : MINIMISATION QUADRATIQUE DES RESTES
!     DES DEVELOPPEMENTS DE TAYLOR DE H(N)
!     A GAUCHE ET A DROITE
!     ------------------------------------------------
!
    h(1) = (f2(nbz2+2, jmod)-f2(nbz2+1, jmod))/(z(2)-z(1))
!
    do n = 2, nbz-1
        h(n) = ( &
               ( &
              f2(n+nbz2+1, jmod)-f2(n+nbz2, jmod))*(z(n+1)-z(n))+(f2(n+nbz2-1, jmod)-f2(n+nbz2, jmo&
               &d))*(z(n-1)-z(n)))/((z(n+1)-z(n))*(z(n+1)-z(n))+(z(n-1)-z(n))*(z(n-1)-z(n) &
               ) &
               )
    end do
!
    h(nbz) = (f2(nbz*jcyl, jmod)-f2(nbz*jcyl-1, jmod))/(z(nbz)-z(nbz-1))
!
    g(1) = (h(2)-h(1))/(z(2)-z(1))
!
    do n = 2, nbz-1
        g(n) = ( &
               ( &
               h(n+1)-h(n))*(z(n+1)-z(n))+(h(n-1)-h(n))*(z(n-1)-z(n)))/((z(n+1)-z(n))*(z(n+1)-&
               &z(n))+(z(n-1)-z(n))*(z(n-1)-z(n) &
               ) &
               )
    end do
!
    g(nbz) = (h(nbz)-h(nbz-1))/(z(nbz)-z(nbz-1))
!
    mefin3 = 0.d0
!
    do n = 1, nbz-1
        mefin3 = mefin3+0.5d0*( &
                 z(n+1)-z(n))*(f3(n+1)*f1(n+nbz1+1, imod)*g(n+1)+f3(n)*f1(n+nbz1, imod)*g(n))
    end do
!
end function
