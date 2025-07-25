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

subroutine hujori(sens, nmat, reorie, angl_naut, vec, &
                  mat)
    implicit none
!
!      DPASSA  -- CALCUL DE LA MATRICE DE PASSAGE DU REPERE
!                 D'ORTHOTROPIE AU REPERE GLOBAL POUR LE
!                 TENSEUR D'ELASTICITE
!                 CETTE MATRICE EST CONSTRUITE EN PARTANT
!                 DE LA CONSIDERATION QUE L'ENERGIE DE DEFORMATION
!                 EXPRIMEE DANS LE REPERE GLOBAL EST EGALE A
!                 L'ENERGIE DE DEFORMATION EXPRIMEE DANS LE REPERE
!                 D'ORTHOTROPIE
!
!   ARGUMENT        E/S  TYPE         ROLE
!    ANGL_NAUT(3)   IN     R        ANGLES NAUTIQUES DEFINISSANT LE REPERE
!                                   D'ORTHOTROPIE
!    IREP           OUT    I        = 0
!                                     SI LE CHANGEMENT DE REPERE EST
!                                     TRIVIAL (I.E. PASSAG = IDENTITE)
!                                   = 1 SINON
!    PASSAG(6,6)    OUT    R        MATRICE DE PASSAGE DU REPERE
!                                   D'ORTHOTROPIE AU REPERE GLOBAL
!                                   POUR LE TENSEUR D'ELASTICITE
!
!========================= DEBUT DES DECLARATIONS ====================
#include "asterf_types.h"
#include "asterfort/matrot.h"
    integer(kind=8) :: i, j, nmat
    real(kind=8) :: zero, deux, dsqr, isqr
    real(kind=8) :: angl_naut(3), p(3, 3), passag(6, 6), passal(6, 6)
    real(kind=8) :: vec1(6), vec(6), mat(6, 6), mat1(6, 6), work(6, 6)
    character(len=5) :: sens
    aster_logical :: reorie
!
    data zero/0.d0/
    data deux/2.d0/
    data dsqr/1.41421356237d0/
    data isqr/.707106781187d0/
!
    if (.not. reorie) goto 9999
!
    do i = 1, 3
        do j = 1, 3
            p(i, j) = zero
        end do
    end do
!
    do i = 1, 6
        vec1(i) = zero
    end do
!
!
! ----   CONSTRUCTION DE LA MATRICE DE PASSAGE (POUR DES VECTEURS)
! ----   DU REPERE D'ORTHOTROPIE AU REPERE GLOBAL
!        ----------------------------------------
    call matrot(angl_naut, p)
!
!
! calcul de PASSAGT * SIG *PASSAG et PASSAGT * DEPS *PASSAG
    if (nmat .eq. 1) then
!
        if (sens .eq. 'LOCAL') then
!
            passal(1, 1) = p(1, 1)*p(1, 1)
            passal(1, 2) = p(1, 2)*p(1, 2)
            passal(1, 3) = p(1, 3)*p(1, 3)
            passal(1, 4) = deux*isqr*p(1, 1)*p(1, 2)
            passal(1, 5) = deux*isqr*p(1, 1)*p(1, 3)
            passal(1, 6) = deux*isqr*p(1, 2)*p(1, 3)
!
            passal(2, 1) = p(2, 1)*p(2, 1)
            passal(2, 2) = p(2, 2)*p(2, 2)
            passal(2, 3) = p(2, 3)*p(2, 3)
            passal(2, 4) = deux*isqr*p(2, 1)*p(2, 2)
            passal(2, 5) = deux*isqr*p(2, 1)*p(2, 3)
            passal(2, 6) = deux*isqr*p(2, 2)*p(2, 3)
!
            passal(3, 1) = p(3, 1)*p(3, 1)
            passal(3, 2) = p(3, 2)*p(3, 2)
            passal(3, 3) = p(3, 3)*p(3, 3)
            passal(3, 4) = deux*isqr*p(3, 1)*p(3, 2)
            passal(3, 5) = deux*isqr*p(3, 1)*p(3, 3)
            passal(3, 6) = deux*isqr*p(3, 2)*p(3, 3)
!
            passal(4, 1) = dsqr*p(1, 1)*p(2, 1)
            passal(4, 2) = dsqr*p(1, 2)*p(2, 2)
            passal(4, 3) = dsqr*p(1, 3)*p(2, 3)
            passal(4, 4) = p(1, 1)*p(2, 2)+p(1, 2)*p(2, 1)
            passal(4, 5) = p(1, 1)*p(2, 3)+p(1, 3)*p(2, 1)
            passal(4, 6) = p(1, 2)*p(2, 3)+p(1, 3)*p(2, 2)
!
            passal(5, 1) = dsqr*p(1, 1)*p(3, 1)
            passal(5, 2) = dsqr*p(1, 2)*p(3, 2)
            passal(5, 3) = dsqr*p(1, 3)*p(3, 3)
            passal(5, 4) = p(1, 1)*p(3, 2)+p(1, 2)*p(3, 1)
            passal(5, 5) = p(1, 1)*p(3, 3)+p(1, 3)*p(3, 1)
            passal(5, 6) = p(1, 2)*p(3, 3)+p(1, 3)*p(3, 2)
!
            passal(6, 1) = dsqr*p(2, 1)*p(3, 1)
            passal(6, 2) = dsqr*p(2, 2)*p(3, 2)
            passal(6, 3) = dsqr*p(2, 3)*p(3, 3)
            passal(6, 4) = p(2, 1)*p(3, 2)+p(2, 2)*p(3, 1)
            passal(6, 5) = p(2, 1)*p(3, 3)+p(2, 3)*p(3, 1)
            passal(6, 6) = p(2, 2)*p(3, 3)+p(2, 3)*p(3, 2)
!
            do i = 1, 6
                do j = 1, 6
                    vec1(i) = vec1(i)+passal(i, j)*vec(j)
                end do
            end do
!
        else if (sens .eq. 'GLOBA') then
!
            passag(1, 1) = p(1, 1)*p(1, 1)
            passag(1, 2) = p(2, 1)*p(2, 1)
            passag(1, 3) = p(3, 1)*p(3, 1)
            passag(1, 4) = deux*isqr*p(1, 1)*p(2, 1)
            passag(1, 5) = deux*isqr*p(1, 1)*p(3, 1)
            passag(1, 6) = deux*isqr*p(2, 1)*p(3, 1)
!
            passag(2, 1) = p(1, 2)*p(1, 2)
            passag(2, 2) = p(2, 2)*p(2, 2)
            passag(2, 3) = p(3, 2)*p(3, 2)
            passag(2, 4) = deux*isqr*p(1, 2)*p(2, 2)
            passag(2, 5) = deux*isqr*p(1, 2)*p(3, 2)
            passag(2, 6) = deux*isqr*p(2, 2)*p(3, 2)
!
            passag(3, 1) = p(1, 3)*p(1, 3)
            passag(3, 2) = p(2, 3)*p(2, 3)
            passag(3, 3) = p(3, 3)*p(3, 3)
            passag(3, 4) = deux*isqr*p(1, 3)*p(2, 3)
            passag(3, 5) = deux*isqr*p(1, 3)*p(3, 3)
            passag(3, 6) = deux*isqr*p(2, 3)*p(3, 3)
!
            passag(4, 1) = dsqr*p(1, 1)*p(1, 2)
            passag(4, 2) = dsqr*p(2, 1)*p(2, 2)
            passag(4, 3) = dsqr*p(3, 1)*p(3, 2)
            passag(4, 4) = p(1, 1)*p(2, 2)+p(2, 1)*p(1, 2)
            passag(4, 5) = p(1, 1)*p(3, 2)+p(3, 1)*p(1, 2)
            passag(4, 6) = p(2, 1)*p(3, 2)+p(3, 1)*p(2, 2)
!
            passag(5, 1) = dsqr*p(1, 1)*p(1, 3)
            passag(5, 2) = dsqr*p(2, 1)*p(2, 3)
            passag(5, 3) = dsqr*p(3, 1)*p(3, 3)
            passag(5, 4) = p(1, 1)*p(2, 3)+p(2, 1)*p(1, 3)
            passag(5, 5) = p(1, 1)*p(3, 3)+p(3, 1)*p(1, 3)
            passag(5, 6) = p(2, 1)*p(3, 3)+p(3, 1)*p(2, 3)
!
            passag(6, 1) = dsqr*p(1, 2)*p(1, 3)
            passag(6, 2) = dsqr*p(2, 2)*p(2, 3)
            passag(6, 3) = dsqr*p(3, 2)*p(3, 3)
            passag(6, 4) = p(1, 2)*p(2, 3)+p(2, 2)*p(1, 3)
            passag(6, 5) = p(1, 2)*p(3, 3)+p(3, 2)*p(1, 3)
            passag(6, 6) = p(2, 2)*p(3, 3)+p(3, 2)*p(2, 3)
!
            do i = 1, 6
                do j = 1, 6
                    vec1(i) = vec1(i)+passag(i, j)*vec(j)
                end do
            end do
!
        end if
!
        do i = 1, 6
            vec(i) = vec1(i)
        end do
!
!
! calcul de PASSAG * DSDE *PASSAL et PASSAG * DEPS *PASSAL
    else if (nmat .eq. 2) then
!
        passal(1, 1) = p(1, 1)*p(1, 1)
        passal(1, 2) = p(1, 2)*p(1, 2)
        passal(1, 3) = p(1, 3)*p(1, 3)
        passal(1, 4) = deux*isqr*p(1, 1)*p(1, 2)
        passal(1, 5) = deux*isqr*p(1, 1)*p(1, 3)
        passal(1, 6) = deux*isqr*p(1, 2)*p(1, 3)
!
        passal(2, 1) = p(2, 1)*p(2, 1)
        passal(2, 2) = p(2, 2)*p(2, 2)
        passal(2, 3) = p(2, 3)*p(2, 3)
        passal(2, 4) = deux*isqr*p(2, 1)*p(2, 2)
        passal(2, 5) = deux*isqr*p(2, 1)*p(2, 3)
        passal(2, 6) = deux*isqr*p(2, 2)*p(2, 3)
!
        passal(3, 1) = p(3, 1)*p(3, 1)
        passal(3, 2) = p(3, 2)*p(3, 2)
        passal(3, 3) = p(3, 3)*p(3, 3)
        passal(3, 4) = deux*isqr*p(3, 1)*p(3, 2)
        passal(3, 5) = deux*isqr*p(3, 1)*p(3, 3)
        passal(3, 6) = deux*isqr*p(3, 2)*p(3, 3)
!
        passal(4, 1) = dsqr*p(1, 1)*p(2, 1)
        passal(4, 2) = dsqr*p(1, 2)*p(2, 2)
        passal(4, 3) = dsqr*p(1, 3)*p(2, 3)
        passal(4, 4) = p(1, 1)*p(2, 2)+p(1, 2)*p(2, 1)
        passal(4, 5) = p(1, 1)*p(2, 3)+p(1, 3)*p(2, 1)
        passal(4, 6) = p(1, 2)*p(2, 3)+p(1, 3)*p(2, 2)
!
        passal(5, 1) = dsqr*p(1, 1)*p(3, 1)
        passal(5, 2) = dsqr*p(1, 2)*p(3, 2)
        passal(5, 3) = dsqr*p(1, 3)*p(3, 3)
        passal(5, 4) = p(1, 1)*p(3, 2)+p(1, 2)*p(3, 1)
        passal(5, 5) = p(1, 1)*p(3, 3)+p(1, 3)*p(3, 1)
        passal(5, 6) = p(1, 2)*p(3, 3)+p(1, 3)*p(3, 2)
!
        passal(6, 1) = dsqr*p(2, 1)*p(3, 1)
        passal(6, 2) = dsqr*p(2, 2)*p(3, 2)
        passal(6, 3) = dsqr*p(2, 3)*p(3, 3)
        passal(6, 4) = p(2, 1)*p(3, 2)+p(2, 2)*p(3, 1)
        passal(6, 5) = p(2, 1)*p(3, 3)+p(2, 3)*p(3, 1)
        passal(6, 6) = p(2, 2)*p(3, 3)+p(2, 3)*p(3, 2)
!
        passag(1, 1) = p(1, 1)*p(1, 1)
        passag(1, 2) = p(2, 1)*p(2, 1)
        passag(1, 3) = p(3, 1)*p(3, 1)
        passag(1, 4) = deux*isqr*p(1, 1)*p(2, 1)
        passag(1, 5) = deux*isqr*p(1, 1)*p(3, 1)
        passag(1, 6) = deux*isqr*p(2, 1)*p(3, 1)
!
        passag(2, 1) = p(1, 2)*p(1, 2)
        passag(2, 2) = p(2, 2)*p(2, 2)
        passag(2, 3) = p(3, 2)*p(3, 2)
        passag(2, 4) = deux*isqr*p(1, 2)*p(2, 2)
        passag(2, 5) = deux*isqr*p(1, 2)*p(3, 2)
        passag(2, 6) = deux*isqr*p(2, 2)*p(3, 2)
!
        passag(3, 1) = p(1, 3)*p(1, 3)
        passag(3, 2) = p(2, 3)*p(2, 3)
        passag(3, 3) = p(3, 3)*p(3, 3)
        passag(3, 4) = deux*isqr*p(1, 3)*p(2, 3)
        passag(3, 5) = deux*isqr*p(1, 3)*p(3, 3)
        passag(3, 6) = deux*isqr*p(2, 3)*p(3, 3)
!
        passag(4, 1) = dsqr*p(1, 1)*p(1, 2)
        passag(4, 2) = dsqr*p(2, 1)*p(2, 2)
        passag(4, 3) = dsqr*p(3, 1)*p(3, 2)
        passag(4, 4) = p(1, 1)*p(2, 2)+p(2, 1)*p(1, 2)
        passag(4, 5) = p(1, 1)*p(3, 2)+p(3, 1)*p(1, 2)
        passag(4, 6) = p(2, 1)*p(3, 2)+p(3, 1)*p(2, 2)
!
        passag(5, 1) = dsqr*p(1, 1)*p(1, 3)
        passag(5, 2) = dsqr*p(2, 1)*p(2, 3)
        passag(5, 3) = dsqr*p(3, 1)*p(3, 3)
        passag(5, 4) = p(1, 1)*p(2, 3)+p(2, 1)*p(1, 3)
        passag(5, 5) = p(1, 1)*p(3, 3)+p(3, 1)*p(1, 3)
        passag(5, 6) = p(2, 1)*p(3, 3)+p(3, 1)*p(2, 3)
!
        passag(6, 1) = dsqr*p(1, 2)*p(1, 3)
        passag(6, 2) = dsqr*p(2, 2)*p(2, 3)
        passag(6, 3) = dsqr*p(3, 2)*p(3, 3)
        passag(6, 4) = p(1, 2)*p(2, 3)+p(2, 2)*p(1, 3)
        passag(6, 5) = p(1, 2)*p(3, 3)+p(3, 2)*p(1, 3)
        passag(6, 6) = p(2, 2)*p(3, 3)+p(3, 2)*p(2, 3)
!
        if (sens .eq. 'LOCAL') then
!
            work = matmul(mat, passag)
            mat1 = matmul(passal, work)
!
        else if (sens .eq. 'GLOBA') then
!
            work = matmul(mat, passal)
            mat1 = matmul(passag, work)
!
        end if
!
        do j = 1, 6
            do i = 1, 6
                mat(i, j) = mat1(i, j)
            end do
        end do
!
    end if
!
9999 continue
end subroutine
