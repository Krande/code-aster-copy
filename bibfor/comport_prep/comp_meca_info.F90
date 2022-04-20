! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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
! aslint: disable=W1003
! person_in_charge: mickael.abbas at edf.fr
!
subroutine comp_meca_info(behaviourPrepPara)
!
use Behaviour_type
!
implicit none
!
#include "asterf_types.h"
#include "asterc/getfac.h"
#include "asterfort/comp_meca_init.h"
!
type(Behaviour_PrepPara), intent(out) :: behaviourPrepPara
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of behaviour (mechanics)
!
! Create datastructure to prepare comportement
!
! --------------------------------------------------------------------------------------------------
!
! Out behaviourPrepPara: datastructure to prepare behaviour
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: factorKeyword = 'COMPORTEMENT'
    integer :: nb_info_comp, nbFactorKeyword
    type(Behaviour_Para) :: behaviourPara
!
! --------------------------------------------------------------------------------------------------
!
    nbFactorKeyword = 0
    call getfac(factorKeyword, nbFactorKeyword)

! - Number of comportement information
    if (nbFactorKeyword .eq. 0) then
        nb_info_comp = 1
    else
        nb_info_comp = nbFactorKeyword
    endif
    behaviourPrepPara%nb_comp = nbFactorKeyword

! - Allocate objects
    allocate(behaviourPrepPara%v_para(nb_info_comp))
    allocate(behaviourPrepPara%v_paraExte(nb_info_comp))

! - If nothing in COMPORTEMENT: all is elastic
    call comp_meca_init(behaviourPara)
    if (nbFactorKeyword .eq. 0) then
        behaviourPrepPara%v_para(1) = behaviourPara
        behaviourPrepPara%v_para(1)%rela_comp = 'ELAS'
        behaviourPrepPara%v_para(1)%defo_comp = 'PETIT'
        behaviourPrepPara%v_para(1)%type_comp = 'COMP_ELAS'
        behaviourPrepPara%v_para(1)%type_cpla = 'ANALYTIQUE'
    endif
!
end subroutine
