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
subroutine wp2biy(lm, lc, lk, s2, dsr, &
                  isi, yh, yb, zh, zb, &
                  lbloq, u1, u2, u3, u4, &
                  n)
    implicit none
#include "asterfort/mrmult.h"
    real(kind=8) :: u1(*), u2(*), u3(*), u4(*), yh(*), yb(*), zh(*), zb(*)
    real(kind=8) :: s2, dsr, isi
    integer(kind=8) :: lm, lc, lk, n, lbloq(*)
!                    T               T
!     CALCUL (ZH  ZB)   = B * (YH YB)
!     OU B EST L' OPERATEUR (REEL) DU PSEUDO PRODUIT SCALAIRE POUR
!     L' APPROCHE EN PARTIE IMAGINAIRE
!     ------------------------------------------------------------------
! IN  LM   : I : MATRICE DE MASSE
! IN  LC   : I : MATRICE D' AMORTISSEMENT
! IN  LK   : I : MATRICE DE RAIDEUR
! IN  DSR  : C : VALEUR DE 2*RE(SHIFT)
! IN  S2   : C : VALEUR DU CARRE DU MODULE DU SHIFT
! IN  ISI  : C : VALEUR DE 1/IM(SHIFT)
! IN  YH   : R : PARTIE SUPERIEUR DE Y
! IN  YB   : R : PARTIE INFERIEURE DE Y
! IN  N    : I : DIMENSION DES MATRICES
! IN  LBLOQ  : I : TYPE DES DDL (LBOLOQ(I) = 0 <=> DDL(I) = BLOQUE)
! OUT ZH   : R : PARTIE SUPERIEURE DU RESULTAT
! OUT ZB   : R : PARTIE INFERIEURE DU RESULTAT
! VAR U1   : R : VECTEUR DE TRAVAIL
! VAR U2   : R : VECTEUR DE TRAVAIL
! VAR U3   : R : VECTEUR DE TRAVAIL
! VAR U4   : R : VECTEUR DE TRAVAIL
!     ------------------------------------------------------------------
    integer(kind=8) :: i
    real(kind=8) :: zero
!     ------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    zero = 0.0d0
!
    call mrmult('ZERO', lk, yh, u1, 1, &
                .false._1)
    call mrmult('ZERO', lc, yh, u2, 1, &
                .false._1)
    call mrmult('ZERO', lm, yb, u3, 1, &
                .false._1)
    call mrmult('ZERO', lm, yh, u4, 1, &
                .false._1)
!
    if (dsr .ne. zero) then
!        --- PARTIE REELLE DU DECALLAGE NON NULLE ---
        do i = 1, n, 1
            zh(i) = -dsr*u1(i)-s2*(u2(i)+u3(i))
            zb(i) = dsr*u3(i)+(-s2*u4(i)+u1(i))*lbloq(i)
        end do
    else
!        --- PARTIE REELLE DU DECALLAGE NULLE ---
        do i = 1, n, 1
            zh(i) = -s2*(u2(i)+u3(i))
            zb(i) = (-s2*u4(i)+u1(i))*lbloq(i)
        end do
    end if
!
    call mrmult('CUMU', lk, yb, zh, 1, &
                .false._1)
    call mrmult('CUMU', lc, yb, zb, 1, &
                .false._1)
!
    do i = 1, n, 1
        zh(i) = isi*zh(i)
        zb(i) = isi*zb(i)
    end do
!
end subroutine
