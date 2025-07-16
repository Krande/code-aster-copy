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
! person_in_charge: mickael.abbas at edf.fr
!
subroutine comp_meca_exc2(l_cristal, l_pmf, &
                          l_excl, vari_excl)
!
    implicit none
!
#include "asterf_types.h"
!
    aster_logical, intent(in) :: l_cristal, l_pmf
    aster_logical, intent(out) :: l_excl
    character(len=16), intent(out) :: vari_excl
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of comportment (mechanics)
!
! Exception for name of internal variables
!
! --------------------------------------------------------------------------------------------------
!
! In  l_cristal        : .true. if *CRISTAL comportment
! In  l_pmf            : .true. if PMF
! Out l_excl           : .true. if exception case (no names for internal variables)
! Out vari_excl        : name of internal variables if l_excl
!
! --------------------------------------------------------------------------------------------------
!
    l_excl = ASTER_FALSE
    vari_excl = ' '

! - Multiple comportment (PMF)
    if (l_pmf) then
        l_excl = ASTER_TRUE
        vari_excl = '&&MULT_PMF'
    end if

! - Multiple comportment
    if (l_cristal) then
        l_excl = ASTER_TRUE
        vari_excl = '&&MULT_COMP'
    end if
!
end subroutine
