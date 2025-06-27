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
subroutine eiinit(nomte, iu, il)
!
!
    implicit none
#include "asterfort/assert.h"
    character(len=16) :: nomte
    integer(kind=8) :: iu(:, :), il(:, :)
! ----------------------------------------------------------------------
!            DECALAGE D'INDICE POUR LES ELEMENTS D'INTERFACE
! ----------------------------------------------------------------------
! IN  NOMTE  NOM DE L'ELEMENT FINI
! OUT IU     DECALAGE D'INDICE POUR ACCEDER AUX DDL DE DEPLACEMENT
! OUT IL     DECALAGE D'INDICE POUR ACCEDER AUX DDL DE LAGRANGE
! ----------------------------------------------------------------------
    integer(kind=8) :: n
    integer(kind=8) :: uh20(16), lh20(4)
    integer(kind=8) :: up15(12), lp15(3)
    integer(kind=8) :: uq8(6), lq8(2)
! ----------------------------------------------------------------------
    data uh20/1, 2, 3, 4, 9, 10, 11, 12, 5, 6, 7, 8, 17, 18, 19, 20/
    data lh20/13, 14, 15, 16/
    data up15/1, 2, 3, 7, 8, 9, 4, 5, 6, 13, 14, 15/
    data lp15/10, 11, 12/
    data uq8/1, 2, 5, 4, 3, 7/
    data lq8/8, 6/
! ----------------------------------------------------------------------
!
    if ((nomte .eq. 'MEEI_HEXA20') .or. (nomte .eq. 'MEEI_HEXS20')) then
        do n = 1, 16
            iu(1, n) = 1+(uh20(n)-1)*3
            iu(2, n) = 2+(uh20(n)-1)*3
            iu(3, n) = 3+(uh20(n)-1)*3
        end do
!
        do n = 1, 4
            il(1, n) = 1+(lh20(n)-1)*3
            il(2, n) = 2+(lh20(n)-1)*3
            il(3, n) = 3+(lh20(n)-1)*3
        end do
!
!
    else if ((nomte .eq. 'MEEI_PENTA15') .or. (nomte .eq. 'MEEI_PENTS15')) &
        then
        do n = 1, 12
            iu(1, n) = 1+(up15(n)-1)*3
            iu(2, n) = 2+(up15(n)-1)*3
            iu(3, n) = 3+(up15(n)-1)*3
        end do
!
        do n = 1, 3
            il(1, n) = 1+(lp15(n)-1)*3
            il(2, n) = 2+(lp15(n)-1)*3
            il(3, n) = 3+(lp15(n)-1)*3
        end do
!
!
    else if ((nomte .eq. 'EIPLQU8') .or. (nomte .eq. 'EIPLQS8') .or. ( &
             nomte .eq. 'EIAXQU8') .or. (nomte .eq. 'EIAXQS8')) then
!
        do n = 1, 6
            iu(1, n) = 1+(uq8(n)-1)*2
            iu(2, n) = 2+(uq8(n)-1)*2
        end do
!
        do n = 1, 2
            il(1, n) = 1+(lq8(n)-1)*2
            il(2, n) = 2+(lq8(n)-1)*2
        end do
!
    else
        ASSERT(.false.)
    end if
!
end subroutine
