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
subroutine comp_meca_info(prepMapCompor)
!
    use BehaviourPrepare_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterc/getfac.h"
#include "asterfort/comp_meca_init.h"
!
    type(BehaviourPrep_MapCompor), intent(out) :: prepMapCompor
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of behaviour (mechanics)
!
! Create datastructure to prepare comportement
!
! --------------------------------------------------------------------------------------------------
!
! Out prepMapCompor    : datastructure to construct COMPOR map
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: factorKeyword = 'COMPORTEMENT'
    integer(kind=8) :: nb_info_comp, nbFactorKeyword
    type(BehaviourPrep_Para) :: prepPara
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
    end if
    prepMapCompor%nb_comp = nbFactorKeyword

! - Allocate objects
    allocate (prepMapCompor%prepPara(nb_info_comp))
    allocate (prepMapCompor%prepExte(nb_info_comp))

! - If nothing in COMPORTEMENT: all is elastic
    call comp_meca_init(prepPara)
    if (nbFactorKeyword .eq. 0) then
        prepMapCompor%prepPara(1) = prepPara
        prepMapCompor%prepPara(1)%rela_comp = 'ELAS'
        prepMapCompor%prepPara(1)%defo_comp = 'PETIT'
        prepMapCompor%prepPara(1)%type_comp = 'COMP_ELAS'
        prepMapCompor%prepPara(1)%type_cpla = 'ANALYTIQUE'
    end if
!
end subroutine
