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

subroutine wp2ayc(lmatra, lmasse, lamor, sigma, lbloq, &
                  yh, yb, zh, zb, u1, &
                  u2, u3, n, solveu)
    implicit none
#include "jeveux.h"
#include "asterfort/mcmult.h"
#include "asterfort/resoud.h"
!
    complex(kind=8) :: sigma, u1(*), u2(*), u3(*), yh(*), yb(*), zh(*), zb(*)
    integer(kind=8) :: lmatra, lmasse, lamor, n, lbloq(*)
    character(len=19) :: solveu
!                   T               T
!     CALCUL (ZH,ZB)   = A * (YH,YB)
!
!     OU A EST L' OPERATEUR (COMPLEXE) DONT ON CHERCHE UNE
!     APPROXIMATION DES ELEMENTS PROPRES
!     ------------------------------------------------------------------
! IN  LMATRA : I : FACTORISEE LDLT (DANS C) DE LA MATRICE DYNAMIQUE
! IN  LMASSE : I : MATRICE DE MASSE
! IN  LAMOR  : I : MATRICE D' AMORTISSEMENT
! IN  SIGMA  : C : VALEUR DU SHIFT
! IN  YH     : C : PARTIE SUPERIEUR DE Y
! IN  YB     : C : PARTIE INFERIEURE DE Y
! IN  N      : I : DIMENSION DE MATRICE
! IN  LBLOQ  : I : TYPE DES DDL (LBOLOQ(I) = 0 <=> DDL(I) = BLOQUE)
! OUT ZH     : C : PARTIE SUPERIEURE DU RESULTAT
! OUT ZB     : C : PARTIE INFERIEURE DU RESULTAT
! VAR U1     : C : VECTEUR DE TRAVAIL, EN SORTIE VAUT AMOR *YH
! VAR U2     : C : VECTEUR DE TRAVAIL, EN SORTIE VAUT MASSE*YB
! VAR U3     : C : VECTEUR DE TRAVAIL, EN SORTIE VAUT MASSE*YH
! VAR V      : C : VECTEUR DE TRAVAIL
! IN  SOLVEU : K19 : SD SOLVEUR POUR PARAMETRER LE SOLVEUR LINEAIRE
!     ------------------------------------------------------------------
!
!
    integer(kind=8) :: i
    character(len=1) :: kbid
    character(len=19) :: k19bid, matass, chcine, criter
    integer(kind=8) :: iret
!     ------------------------------------------------------------------
!
! INIT. OBJETS ASTER
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    matass = zk24(zi(lmatra+1))
    chcine = ' '
    criter = ' '
    k19bid = ' '
!
    call mcmult('ZERO', lamor, yh, u1, 1, &
                .false._1)
    call mcmult('ZERO', lmasse, yb, u2, 1, &
                .false._1)
    call mcmult('ZERO', lmasse, yh, u3, 1, &
                .false._1)
!-RM-DEB
!     LA BOUCLE 5 REALISE LE PRODUIT PAR MASSE*INV(MASSE_REG)*MASSR
!     OR CETTE MATRICE EST EGALE A MASSE
!---> VOIR CE QUI SE PASSE QUAND LA BOUCLE EST SUPPRIMEE
    do i = 1, n, 1
        u3(i) = u3(i)*lbloq(i)
        u2(i) = u2(i)*lbloq(i)
    end do
!-RM-FIN
    do i = 1, n, 1
        u1(i) = u1(i)+sigma*u3(i)+u2(i)
    end do
    call resoud(matass, k19bid, solveu, chcine, 1, &
                k19bid, k19bid, kbid, [0.d0], u1, &
                criter, .false._1, 0, iret)
    do i = 1, n, 1
        zh(i) = -u1(i)
        zb(i) = (yh(i)-sigma*u1(i))*lbloq(i)
    end do
end subroutine
