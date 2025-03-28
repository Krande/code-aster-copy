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
subroutine nume_ddl_matr(numeDofZ, jvListOfMatrZ, nbMatrElem)
!
    implicit none
!
#include "asterfort/as_deallocate.h"
#include "asterfort/jeveuo.h"
#include "asterfort/numddl.h"
#include "asterfort/promor.h"
!
    character(len=*), intent(in) :: numeDofZ, jvListOfMatrZ
    integer, intent(in) :: nbMatrElem
!
! --------------------------------------------------------------------------------------------------
!
! Factor
!
! Numbering - Create NUME_EQUA objects with matrix
!
! --------------------------------------------------------------------------------------------------
!
! In  numeDof       : name of numeDof object
! In  jvListOfMatr  : name of JEVEUX name for list of elementary matrixes
! In  nbMatrElem    : number of elementary matrixes
!
! --------------------------------------------------------------------------------------------------
!
    character(len=24), parameter :: renumSans = "SANS"
    character(len=14) :: numeDof
    character(len=24), pointer :: listMatrElem(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    numeDof = numeDofZ

! - Get list of matrix
    call jeveuo(jvListOfMatrZ, 'L', vk24=listMatrElem)

! - CALCUL DE LA NUMEROTATION PROPREMENT DITE :
    call numddl(numeDof, renumSans, 'GG', nbMatrElem, listMatrElem)

! - CREATION ET CALCUL DU STOCKAGE MORSE DE LA MATRICE :
    call promor(numeDof, 'G')
!
end subroutine
