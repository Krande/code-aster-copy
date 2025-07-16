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
subroutine dil2gr(imate, ndim, dimdef, &
                  defgep, sigp, dsde2g)
! --- BUT : CALCUL DE LA LOI DE COMPORTEMENT POUR LA PARTIE --
! ---       SECOND GRADIENT --------------------------------------------
! ======================================================================
    implicit none
#include "asterfort/rcvalb.h"
    integer(kind=8) :: imate, ndim, dimdef
    real(kind=8) :: sigp(ndim), dsde2g(ndim, ndim), defgep(dimdef)
! ======================================================================
! --- VARIABLES LOCALES ------------------------------------------------
! ======================================================================
    integer(kind=8) :: i, j
    real(kind=8) :: val(5)
    integer(kind=8) :: icodre(5), kpg, spt
    character(len=8) :: ncra(5), fami, poum
! ======================================================================
! --- DEFINITION DES DONNEES INITIALES ---------------------------------
! ======================================================================
    data ncra/'A1', 'A2', 'A3', 'A4', 'A5'/
    fami = 'FPG1'
    kpg = 1
    spt = 1
    poum = '+'

    do i = 1, ndim
        do j = 1, ndim
            dsde2g(j, i) = 0.0d0
        end do
    end do

    call rcvalb(fami, kpg, spt, poum, imate, &
                ' ', 'SECOND_ELAS', 0, ' ', [0.0d0], &
                1, ncra(1), val(1), icodre(1), 1)

    do i = 1, ndim
        dsde2g(i, i) = (1+ndim)*val(1)
    end do
!
    do i = 1, ndim
        sigp(i) = 0.0d0
        do j = 1, ndim
            sigp(i) = sigp(i)+dsde2g(i, j)*defgep(j)
        end do
    end do

! ======================================================================
end subroutine
