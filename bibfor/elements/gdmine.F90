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
subroutine gdmine(kp, nno, pjacob, en, grani, &
                  alfnmk, delnmk, pas, rot0, rotm, &
                  rotkm1, rotk, rmkm1, rmk, omgkm, &
                  ompgkm, omgk, ompgk, rigi)
!
! FONCTION: POUR UN ELEMENT DE POUTRE EN GRAND DEPLACEMENT, CALCULE LA
!           CONTRIBUTION DU POINT DE GAUSS KP A LA PARTIE ROTATOIRE DE
!           LA MATRICE D'INERTIE.  CETTE MATRICE EST DIRECTEMENT AJOUTEE
!           A LA MATRICE DE RIGIDITE.
!
!     IN  : KP        : NUMERO DU POINT DE GAUSS
!           NNO       : NOMBRE DE NOEUDS
!           PJACOB    : POIDS * JACOBIEN
!           EN        : FONCTIONS DE FORME
!           GRANI     : DIAGONALE DU TENSEUR D'INERTIE EN AXES LOCAUX
!                       POUR LES 3 1ERES COMPOSANTES, RHO*A POUR LA 4EME
!           ALFNMK    : ALPHA DE NEWMARK
!           DELNMK    : DELTA DE NEWMARK
!           PAS       : PAS DE TEMPS
!           ROT0      : MATRICE DE ROTATION DES AXES PRINCIPAUX D'INERT.
!                       AU POINT DE GAUSS DANS LA POSITION DE REFERENCE,
!                       PAR RAPPORT AUX AXES GENERAUX
!           ROTM      : MATRICE DE ROTATION A L'INSTANT N
!           ROTKM1    : MATRICE DE ROTATION A L'ITER. I   DE L'INST. N+1
!           ROTK      : MATRICE DE ROTATION A L'ITER. I+1 DE L'INST. N+1
!           RMKM1     : INCREMENT DE ROTATION ENTRE L'INSTANT N  ET
!                       L'ITER. I   DE L'INSTANT N+1
!           RMK       : INCREMENT DE ROTATION ENTRE L'INSTANT N  ET
!                       L'ITER. I+1 DE L'INSTANT N+1
!           OMGKM     : VITESSE      ANGULAIRE, A L'ITERATION I
!           OMPGKM    : ACCELERATION ANGULAIRE, A L'ITERATION I
!           OMGK      : VITESSE      ANGULAIRE, A L'ITERATION I+1
!           OMPGK     : ACCELERATION ANGULAIRE, A L'ITERATION I+1
!
!     OUT : RIGI      : MATRICE DE RIGIDITE (CUMUL DES CONTRIBUTIONS DES
!                       POINTS DE GAUSS)
! ------------------------------------------------------------------
    implicit none
#include "asterfort/antisy.h"
#include "asterfort/cumuma.h"
#include "asterfort/gdt.h"
#include "asterfort/promat.h"
    real(kind=8) :: inert(6, 6), iro(3, 3), iroomt(3, 3)
    real(kind=8) :: en(3, 2), grani(4), rot0(3, 3), rotm(3, 3), rotkm1(3, 3)
    real(kind=8) :: rotk(3, 3), rmkm1(3), rmk(3), omgkm(3), ompgkm(3), omgk(3)
    real(kind=8) :: ompgk(3), rigi(18, 18), rotabs(3, 3), rotmt(3, 3)
    real(kind=8) :: rotkmt(3, 3), omtiro(3, 3), amat1(3, 3), amat2(3, 3)
    real(kind=8) :: amat3(3, 3), amat4(3, 3), amat5(3, 3), amat6(3, 3)
    real(kind=8) :: amat7(3, 3), vect1(3), vect2(3)
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, j, kp, nno
    real(kind=8) :: alfnmk, coef, delnmk, pas, pjacob, un, zero
!
!-----------------------------------------------------------------------
    zero = 0.d0
    un = 1.d0
!
    call promat(rotk, 3, 3, 3, rot0, &
                3, 3, 3, rotabs)
!
    do j = 1, 3
        do i = 1, 3
            amat2(i, j) = grani(i)*rotabs(j, i)
        end do
    end do
    call promat(rotabs, 3, 3, 3, amat2, &
                3, 3, 3, iro)
    call promat(iro, 3, 3, 3, omgk, &
                3, 3, 1, vect1)
    call antisy(vect1, un, iroomt)
!
    call antisy(omgk, un, amat2)
    call promat(amat2, 3, 3, 3, iro, &
                3, 3, 3, omtiro)
    do j = 1, 3
        do i = 1, 3
            amat5(i, j) = -iroomt(i, j)+omtiro(i, j)
            amat1(i, j) = (iro(i, j)+delnmk*pas*amat5(i, j))/alfnmk/pas/pas
        end do
    end do
!
    do j = 1, 3
        do i = 1, 3
            rotmt(i, j) = rotm(j, i)
            rotkmt(i, j) = rotkm1(j, i)
        end do
        vect1(j) = rmk(j)-rmkm1(j)
    end do
!
    call promat(rotmt, 3, 3, 3, vect1, &
                3, 3, 1, vect2)
    call promat(rotk, 3, 3, 3, vect2, &
                3, 3, 1, vect1)
    coef = -un
    call antisy(vect1, coef, amat2)
!
    call gdt(rmk, amat3)
    call promat(rotmt, 3, 3, 3, amat3, &
                3, 3, 3, amat4)
    call promat(rotk, 3, 3, 3, amat4, &
                3, 3, 3, amat3)
!
    do j = 1, 3
        do i = 1, 3
            amat2(i, j) = amat2(i, j)+amat3(i, j)
        end do
    end do
!
    call promat(amat1, 3, 3, 3, amat2, &
                3, 3, 3, amat3)
!C
    call promat(rotkmt, 3, 3, 3, ompgkm, &
                3, 3, 1, vect1)
    call promat(rotk, 3, 3, 3, vect1, &
                3, 3, 1, vect2)
    call antisy(vect2, coef, amat2)
    call promat(iro, 3, 3, 3, amat2, &
                3, 3, 3, amat1)
!C
    call promat(rotkmt, 3, 3, 3, omgkm, &
                3, 3, 1, vect1)
    call promat(rotk, 3, 3, 3, vect1, &
                3, 3, 1, vect2)
    call antisy(vect2, coef, amat4)
    call promat(amat5, 3, 3, 3, amat4, &
                3, 3, 3, amat2)
!CC
    call promat(iro, 3, 3, 3, ompgk, &
                3, 3, 1, vect1)
    call antisy(vect1, un, amat4)
    call antisy(ompgk, un, amat5)
    call promat(iro, 3, 3, 3, amat5, &
                3, 3, 3, amat6)
    do j = 1, 3
        do i = 1, 3
            amat4(i, j) = -amat4(i, j)+amat6(i, j)
        end do
    end do
    call antisy(omgk, un, amat5)
    call promat(amat5, 3, 3, 3, iroomt, &
                3, 3, 3, amat6)
    call promat(omtiro, 3, 3, 3, amat5, &
                3, 3, 3, amat7)
    do j = 1, 3
        do i = 1, 3
            amat5(i, j) = -amat6(i, j)+amat7(i, j)
        end do
    end do
!CC
    do j = 1, 6
        do i = 1, 6
            inert(i, j) = zero
        end do
    end do
!
!* ON CALCULE CE QU'APPORTE LA MATRICE DE MASSE A LA PARTIE ROTATOIRE,
!* POUR LE RETRANCHER ICI.
    do j = 1, 3
        do i = 1, 3
            amat6(i, j) = grani(i)*rot0(j, i)
        end do
    end do
    call promat(rot0, 3, 3, 3, amat6, &
                3, 3, 3, amat7)
!
    coef = un/alfnmk/pas/pas
    do j = 1, 3
        do i = 1, 3
            inert(3+i, 3+j) = amat3(i, j)+amat1(i, j)+amat2(i, j)+amat4(i, j)+amat5(i, j)-coe&
                             &f*amat7(i, j)
!C    INERT(3+I,3+J) = COEF * ( AMAT3(I,J) + AMAT1(I,J) + AMAT2(I,J) +
!C   &                          AMAT4(I,J) + AMAT5(I,J) - AMAT7(I,J) )
!    &                          AMAT4(I,J) + AMAT5(I,J)              )
        end do
    end do
    do j = 1, nno
        do i = 1, nno
            coef = pjacob*en(i, kp)*en(j, kp)
            call cumuma(i, j, inert, coef, rigi)
        end do
    end do
end subroutine
