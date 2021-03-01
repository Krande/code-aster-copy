! --------------------------------------------------------------------
! Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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

subroutine nume_ddl_matr(nume, list_matr, nb_matr)
!
implicit none
!
#include "asterfort/promor.h"
#include "asterfort/numddl.h"
#include "asterfort/jeveuo.h"
!
!
    character(len=*), intent(in) :: nume
    character(len=*), intent(in) :: list_matr
    integer, intent(in) :: nb_matr
!
! --------------------------------------------------------------------------------------------------
!
! Factor
!
! Numbering - Create NUME_EQUA objects with matrix
!
! --------------------------------------------------------------------------------------------------
!
! In  nume_ddl       : name of nume_ddl object
! In  list_matr      : list of elementary matrixes
! In  nb_matr        : number of elementary matrixes
!                       SANS/RCMKs
!
! --------------------------------------------------------------------------------------------------
!
    character(len=14) :: nume_ddl
    character(len=24), pointer :: v_matr(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    nume_ddl = nume
!
! ----- LISTE DES CHARGES
!
    call jeveuo(list_matr, 'L', vk24=v_matr)
!
! ----- CALCUL DE LA NUMEROTATION PROPREMENT DITE :
!
    call numddl(nume_ddl, 'GG', nb_matr, v_matr)
!
! ----- CREATION ET CALCUL DU STOCKAGE MORSE DE LA MATRICE :
!
    call promor(nume_ddl, 'G')
!
!
end subroutine
