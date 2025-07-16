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
subroutine comp_meca_elas(compElas, l_etat_init)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nocart.h"
#include "asterfort/Behaviour_type.h"
!
    character(len=19), intent(in) :: compElas
    aster_logical, intent(in) :: l_etat_init
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of comportment (mechanics)
!
! Set elastic comportment
!
! --------------------------------------------------------------------------------------------------
!
! In  compElas    : name of ELAS <CARTE> COMPOR
! In  l_etat_init : .true. if initial state is defined
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbCmp = COMPOR_SIZE
    character(len=16), pointer :: compElasValv(:) => null()
!
! --------------------------------------------------------------------------------------------------
!

! - Access to map
    call jeveuo(compElas(1:19)//'.VALV', 'E', vk16=compElasValv)

! - Init <CARTE>
    compElasValv(1:COMPOR_SIZE) = 'VIDE'

! - Set for ELASTIQUE
    compElasValv(RELA_NAME) = 'ELAS'
    compElasValv(NVAR) = '1'
    compElasValv(DEFO) = 'PETIT'
    if (l_etat_init) then
        compElasValv(INCRELAS) = 'COMP_INCR'
    else
        compElasValv(INCRELAS) = 'COMP_ELAS'
    end if
    compElasValv(PLANESTRESS) = 'ANALYTIQUE'
    write (compElasValv(NUME), '(I16)') 1
    write (compElasValv(KIT1_NVAR), '(I16)') 1
    write (compElasValv(KIT2_NVAR), '(I16)') 1
    write (compElasValv(KIT3_NVAR), '(I16)') 1
    write (compElasValv(KIT4_NVAR), '(I16)') 1

! - Create <CARTE>
    call nocart(compElas, 1, nbCmp)
!
end subroutine
