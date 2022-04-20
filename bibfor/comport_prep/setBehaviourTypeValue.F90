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
! person_in_charge: mickael.abbas at edf.fr
!
subroutine setBehaviourTypeValue(behaviourPara, iFactorKeyword_,&
                                 comporList_, comporMap_)
!
use Behaviour_type
!
implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/comp_meca_l.h"
#include "asterfort/Behaviour_type.h"
!
type(Behaviour_Para), pointer :: behaviourPara(:)
integer, optional, intent(in) :: iFactorKeyword_
character(len=16), intent(out), optional :: comporList_(:)
character(len=16), pointer, optional :: comporMap_(:)
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of constitutive laws (mechanics)
!
! Set values in the map or in list
!
! --------------------------------------------------------------------------------------------------
!
! In  behaviourPara    : parameters for constitutive laws
! In  iFactorKeyword   : index of factor keyword (for map)
! In  comporList       : list for parameters of constitutive laws
! In  comporMap        : map for parameters of constitutive laws
!
! --------------------------------------------------------------------------------------------------
!
    integer :: iFactorKeyword
    aster_logical :: l_pmf, l_kit_thm, l_kit_ddi, l_kit_meta, l_kit_cg, l_exte_comp
!
! --------------------------------------------------------------------------------------------------
!
    iFactorKeyword = 1
    if (present(iFactorKeyword_)) then
        iFactorKeyword = iFactorKeyword_
    endif
!
    call comp_meca_l(behaviourPara(iFactorKeyword)%rela_comp, 'PMF'     , l_pmf)
    call comp_meca_l(behaviourPara(iFactorKeyword)%rela_comp, 'KIT_THM' , l_kit_thm)
    call comp_meca_l(behaviourPara(iFactorKeyword)%rela_comp, 'KIT_DDI' , l_kit_ddi)
    call comp_meca_l(behaviourPara(iFactorKeyword)%rela_comp, 'KIT_META', l_kit_meta)
    call comp_meca_l(behaviourPara(iFactorKeyword)%rela_comp, 'KIT_CG'  , l_kit_cg)
    call comp_meca_l(behaviourPara(iFactorKeyword)%rela_comp, 'EXTE_COMP', l_exte_comp)
!
    if (present(comporMap_)) then
        comporMap_(1:COMPOR_SIZE) = 'VIDE'
        comporMap_(RELA_NAME) = behaviourPara(iFactorKeyword)%rela_comp
        write (comporMap_(NVAR),'(I16)') behaviourPara(iFactorKeyword)%nbVari
        comporMap_(DEFO) = behaviourPara(iFactorKeyword)%defo_comp
        comporMap_(INCRELAS) = behaviourPara(iFactorKeyword)%type_comp
        comporMap_(PLANESTRESS) = behaviourPara(iFactorKeyword)%type_cpla
        if (.not.l_pmf) then
            write (comporMap_(NUME),'(I16)') behaviourPara(iFactorKeyword)%numeLaw
        endif
        comporMap_(MULTCOMP) = behaviourPara(iFactorKeyword)%mult_comp
        comporMap_(POSTITER) = behaviourPara(iFactorKeyword)%post_iter
        comporMap_(DEFO_LDC) = behaviourPara(iFactorKeyword)%defo_ldc
        comporMap_(RIGI_GEOM) = behaviourPara(iFactorKeyword)%rigi_geom
        comporMap_(REGUVISC) = behaviourPara(iFactorKeyword)%regu_visc
        if (l_kit_thm) then
            comporMap_(MECA_NAME) = behaviourPara(iFactorKeyword)%kit_comp(1)
            comporMap_(HYDR_NAME) = behaviourPara(iFactorKeyword)%kit_comp(2)
            comporMap_(THER_NAME) = behaviourPara(iFactorKeyword)%kit_comp(3)
            comporMap_(THMC_NAME) = behaviourPara(iFactorKeyword)%kit_comp(4)
            write (comporMap_(THMC_NUME),'(I16)') behaviourPara(iFactorKeyword)%numeLawKit(1)
            write (comporMap_(THER_NUME),'(I16)') behaviourPara(iFactorKeyword)%numeLawKit(2)
            write (comporMap_(HYDR_NUME),'(I16)') behaviourPara(iFactorKeyword)%numeLawKit(3)
            write (comporMap_(MECA_NUME),'(I16)') behaviourPara(iFactorKeyword)%numeLawKit(4)
            write (comporMap_(THMC_NVAR),'(I16)') behaviourPara(iFactorKeyword)%nbVariKit(1)
            write (comporMap_(THER_NVAR),'(I16)') behaviourPara(iFactorKeyword)%nbVariKit(2)
            write (comporMap_(HYDR_NVAR),'(I16)') behaviourPara(iFactorKeyword)%nbVariKit(3)
            write (comporMap_(MECA_NVAR),'(I16)') behaviourPara(iFactorKeyword)%nbVariKit(4)
        endif
        if (l_kit_ddi) then
            comporMap_(CREEP_NAME) = behaviourPara(iFactorKeyword)%kit_comp(1)
            comporMap_(PLAS_NAME)  = behaviourPara(iFactorKeyword)%kit_comp(2)
            comporMap_(COUPL_NAME) = behaviourPara(iFactorKeyword)%kit_comp(3)
            comporMap_(CPLA_NAME)  = behaviourPara(iFactorKeyword)%kit_comp(4)
            write (comporMap_(CREEP_NUME),'(I16)') behaviourPara(iFactorKeyword)%numeLawKit(1)
            write (comporMap_(PLAS_NUME),'(I16)')  behaviourPara(iFactorKeyword)%numeLawKit(2)
            write (comporMap_(CREEP_NVAR),'(I16)') behaviourPara(iFactorKeyword)%nbVariKit(1)
            write (comporMap_(PLAS_NVAR),'(I16)')  behaviourPara(iFactorKeyword)%nbVariKit(2)
        endif
        comporMap_(KIT1_NAME) = behaviourPara(iFactorKeyword)%kit_comp(1)
        if (l_kit_meta) then
            comporMap_(META_PHAS)  = behaviourPara(iFactorKeyword)%kit_comp(1)
            comporMap_(META_RELA)  = behaviourPara(iFactorKeyword)%kit_comp(2)
            comporMap_(META_GLOB)  = behaviourPara(iFactorKeyword)%kit_comp(3)
        endif
        if (l_kit_cg) then
            comporMap_(CABLE_NAME)   = behaviourPara(iFactorKeyword)%kit_comp(1)
            comporMap_(SHEATH_NAME)  = behaviourPara(iFactorKeyword)%kit_comp(2)
            write (comporMap_(CABLE_NUME),'(I16)') behaviourPara(iFactorKeyword)%numeLawKit(1)
            write (comporMap_(SHEATH_NUME),'(I16)')  behaviourPara(iFactorKeyword)%numeLawKit(2)
            write (comporMap_(CABLE_NVAR),'(I16)') behaviourPara(iFactorKeyword)%nbVariKit(1)
            write (comporMap_(SHEATH_NVAR),'(I16)')  behaviourPara(iFactorKeyword)%nbVariKit(2)
        endif
        if (l_exte_comp) then
            write (comporMap_(MECA_NVAR),'(I16)') behaviourPara(iFactorKeyword)%nbVariKit(4)
        endif
    endif
    if (present(comporList_)) then
        comporList_(1:COMPOR_SIZE) = 'VIDE'
        comporList_(RELA_NAME) = behaviourPara(iFactorKeyword)%rela_comp
        write (comporList_(NVAR),'(I16)') behaviourPara(iFactorKeyword)%nbVari
        comporList_(DEFO) = behaviourPara(iFactorKeyword)%defo_comp
        comporList_(INCRELAS) = behaviourPara(iFactorKeyword)%type_comp
        comporList_(PLANESTRESS) = behaviourPara(iFactorKeyword)%type_cpla
        if (.not.l_pmf) then
            write (comporList_(NUME),'(I16)') behaviourPara(iFactorKeyword)%numeLaw
        endif
        comporList_(MULTCOMP) = behaviourPara(iFactorKeyword)%mult_comp
        comporList_(POSTITER) = behaviourPara(iFactorKeyword)%post_iter
        comporList_(DEFO_LDC) = behaviourPara(iFactorKeyword)%defo_ldc
        comporList_(RIGI_GEOM) = behaviourPara(iFactorKeyword)%rigi_geom
        comporList_(REGUVISC) = behaviourPara(iFactorKeyword)%regu_visc
        if (l_kit_thm) then
            ASSERT(ASTER_FALSE)
        endif
        if (l_kit_ddi) then
            comporList_(CREEP_NAME) = behaviourPara(iFactorKeyword)%kit_comp(1)
            comporList_(PLAS_NAME)  = behaviourPara(iFactorKeyword)%kit_comp(2)
            comporList_(COUPL_NAME) = behaviourPara(iFactorKeyword)%kit_comp(3)
            comporList_(CPLA_NAME)  = behaviourPara(iFactorKeyword)%kit_comp(4)
            write (comporList_(CREEP_NUME),'(I16)') behaviourPara(iFactorKeyword)%numeLawKit(1)
            write (comporList_(PLAS_NUME),'(I16)')  behaviourPara(iFactorKeyword)%numeLawKit(2)
            write (comporList_(CREEP_NVAR),'(I16)') behaviourPara(iFactorKeyword)%nbVariKit(1)
            write (comporList_(PLAS_NVAR),'(I16)')  behaviourPara(iFactorKeyword)%nbVariKit(2)
        endif
        comporList_(KIT1_NAME) = behaviourPara(iFactorKeyword)%kit_comp(1)
        if (l_kit_meta) then
            comporList_(META_PHAS)  = behaviourPara(iFactorKeyword)%kit_comp(1)
            comporList_(META_RELA)  = behaviourPara(iFactorKeyword)%kit_comp(2)
            comporList_(META_GLOB)  = behaviourPara(iFactorKeyword)%kit_comp(3)
        endif
        if (l_kit_cg) then
            comporList_(CABLE_NAME)   = behaviourPara(iFactorKeyword)%kit_comp(1)
            comporList_(SHEATH_NAME)  = behaviourPara(iFactorKeyword)%kit_comp(2)
            write (comporList_(CABLE_NUME),'(I16)') behaviourPara(iFactorKeyword)%numeLawKit(1)
            write (comporList_(SHEATH_NUME),'(I16)')  behaviourPara(iFactorKeyword)%numeLawKit(2)
            write (comporList_(CABLE_NVAR),'(I16)') behaviourPara(iFactorKeyword)%nbVariKit(1)
            write (comporList_(SHEATH_NVAR),'(I16)')  behaviourPara(iFactorKeyword)%nbVariKit(2)
        endif
        if (l_exte_comp) then
            write (comporList_(MECA_NVAR),'(I16)') behaviourPara(iFactorKeyword)%nbVariKit(4)
        endif
    endif
!
end subroutine
