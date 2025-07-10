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
! aslint: disable=W1403
! person_in_charge: mickael.abbas at edf.fr
!
subroutine comp_meca_init(prepPara)
!
    use BehaviourPrepare_type
!
    implicit none
!
!
    type(BehaviourPrep_Para), intent(out) :: prepPara
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of behaviour (mechanics)
!
! Init datastructure to describe comportement
!
! --------------------------------------------------------------------------------------------------
!
! Out prepPara         : behaviour parameters from user
!
! --------------------------------------------------------------------------------------------------
!
    prepPara%rela_comp = 'VIDE'
    prepPara%defo_comp = 'VIDE'
    prepPara%type_comp = 'VIDE'
    prepPara%type_cpla = 'VIDE'
    prepPara%kit_comp = 'VIDE'
    prepPara%mult_comp = 'VIDE'
    prepPara%post_iter = 'VIDE'
    prepPara%defo_ldc = 'VIDE'
    prepPara%rigi_geom = 'VIDE'
    prepPara%regu_visc = 'VIDE'
    prepPara%post_incr = 'VIDE'
    prepPara%nbVari = 0
    prepPara%numeLaw = 0
    prepPara%nbVariKit = 0
    prepPara%numeLawKit = 0
!
end subroutine
