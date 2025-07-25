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

subroutine wp2ayl(appr, lmatra, lmasse, lamor, sigma, &
                  lbloq, yh, yb, zh, zb, &
                  u1, u2, u3, u4, v, &
                  n, solveu)
    implicit none
#include "jeveux.h"
#include "asterfort/mrmult.h"
#include "asterfort/resoud.h"
!
    character(len=1) :: appr
    complex(kind=8) :: v(*), sigma
    real(kind=8) :: u1(*), u2(*), u3(*), u4(*), yh(*), yb(*), zh(*), zb(*)
    integer(kind=8) :: lmatra, lmasse, lamor, n, lbloq(*)
    character(len=19) :: solveu
!                   T               T
!     CALCUL (ZH,ZB)   = A * (YH,YB)
!
!     OU A EST L' OPERATEUR (REEL) DONT ON CHERCHE UNE APPROXIMATION
!     DES ELEMENTS PROPRES
!     ------------------------------------------------------------------
! IN  APPR   : K : INDICATEUR D' APPROCHE ('R'/'I') POUR A
! IN  LMATRA : I : FACTORISEE LDLT (DANS C) DE LA MATRICE DYNAMIQUE
! IN  LMASSE : I : MATRICE DE MASSE
! IN  LAMOR  : I : MATRICE D' AMORTISSEMENT
! IN  SIGMA  : C : VALEUR DU SHIFT
! IN  YH     : R : PARTIE SUPERIEUR DE Y
! IN  YB     : R : PARTIE INFERIEURE DE Y
! IN  N      : I : DIMENSION DE MATRICE
! IN  LBLOQ  : I : TYPE DES DDL (LBOLOQ(I) = 0 <=> DDL(I) = BLOQUE)
! OUT ZH     : R : PARTIE SUPERIEURE DU RESULTAT
! OUT ZB     : R : PARTIE INFERIEURE DU RESULTAT
! VAR U1     : R : VECTEUR DE TRAVAIL, EN SORTIE VAUT AMOR *YH
! VAR U2     : R : VECTEUR DE TRAVAIL, EN SORTIE VAUT MASSE*YB
! VAR U3     : R : VECTEUR DE TRAVAIL, EN SORTIE VAUT MASSE*YH
! VAR U4     : R : VECTEUR DE TRAVAIL
! VAR V      : R : VECTEUR DE TRAVAIL
! IN  SOLVEU : K19: SD SOLVEUR POUR PARAMETRER LE SOLVEUR LINEAIRE
!     ------------------------------------------------------------------
!
!
    real(kind=8) :: zero, sr
    integer(kind=8) :: i
    complex(kind=8) :: cbid
    character(len=1) :: kbid
    character(len=19) :: k19bid, matass, chcine, criter
    integer(kind=8) :: iret
!     ------------------------------------------------------------------
! INIT. OBJETS ASTER
!-----------------------------------------------------------------------
    real(kind=8) :: si
    cbid = dcmplx(0.d0, 0.d0)
!-----------------------------------------------------------------------
    matass = zk24(zi(lmatra+1))
    chcine = ' '
    criter = ' '
    k19bid = ' '
    zero = 0.0d0
    sr = dble(sigma)
    si = dimag(sigma)
!
    call mrmult('ZERO', lamor, yh, u1, 1, &
                .false._1)
    call mrmult('ZERO', lmasse, yb, u2, 1, &
                .false._1)
    call mrmult('ZERO', lmasse, yh, u3, 1, &
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
    if (si .ne. zero) then
        do i = 1, n, 1
            v(i) = dcmplx(u1(i))+sigma*dcmplx(u3(i))+dcmplx(u2(i))
        end do
        call resoud(matass, k19bid, solveu, chcine, 1, &
                    k19bid, k19bid, kbid, [0.d0], v, &
                    criter, .false._1, 0, iret)
        if (appr .eq. 'R') then
            do i = 1, n, 1
                zh(i) = -dble(v(i))
                zb(i) = (yh(i)-dble(sigma*v(i)))*lbloq(i)
            end do
        else
            do i = 1, n, 1
                zh(i) = -dimag(v(i))
                zb(i) = -dimag(sigma*v(i))*lbloq(i)
            end do
        end if
    else
        do i = 1, n, 1
            u4(i) = u1(i)+sr*u3(i)+u2(i)
        end do
        call resoud(matass, k19bid, solveu, chcine, 1, &
                    k19bid, k19bid, kbid, u4, [cbid], &
                    criter, .false._1, 0, iret)
        do i = 1, n, 1
            zh(i) = -u4(i)
            zb(i) = (yh(i)-sr*u4(i))*lbloq(i)
        end do
    end if
end subroutine
