! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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
subroutine comp_meca_cvar(behaviourPrepPara)
!
    use Behaviour_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/comp_nbvari.h"
!
    type(Behaviour_PrepPara), intent(inout) :: behaviourPrepPara
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of constitutive laws (mechanics)
!
! Count all internal state variables
!
! --------------------------------------------------------------------------------------------------
!
! IO  behaviourPrepPara: datastructure to prepare comportement
!
! --------------------------------------------------------------------------------------------------
!
    integer :: iFactorKeyword, nbFactorKeyword
    character(len=16) :: post_iter, extern_addr
    character(len=16) :: rela_comp, defo_comp, mult_comp, kit_comp(4), type_cpla, regu_visc
    integer :: numeLawKit(4), nbVari, nbVariKit(4), numeLaw, nbVariUMAT, model_dim
    integer :: extern_type
!
! --------------------------------------------------------------------------------------------------
!
    nbFactorKeyword = behaviourPrepPara%nb_comp
    do iFactorKeyword = 1, nbFactorKeyword
        nbVari = 0
        numeLaw = 0
        nbVariKit = 0
        numeLawKit = 0

! ----- Get parameters
        rela_comp = behaviourPrepPara%v_para(iFactorKeyword)%rela_comp
        defo_comp = behaviourPrepPara%v_para(iFactorKeyword)%defo_comp
        type_cpla = behaviourPrepPara%v_para(iFactorKeyword)%type_cpla
        kit_comp = behaviourPrepPara%v_para(iFactorKeyword)%kit_comp
        mult_comp = behaviourPrepPara%v_para(iFactorKeyword)%mult_comp
        post_iter = behaviourPrepPara%v_para(iFactorKeyword)%post_iter
        regu_visc = behaviourPrepPara%v_para(iFactorKeyword)%regu_visc
        nbVariUMAT = behaviourPrepPara%v_paraExte(iFactorKeyword)%nbVariUMAT
        extern_addr = behaviourPrepPara%v_paraExte(iFactorKeyword)%extern_addr
        extern_type = behaviourPrepPara%v_paraExte(iFactorKeyword)%extern_type
        model_dim = behaviourPrepPara%v_paraExte(iFactorKeyword)%model_dim

! ----- Count the number of internal state variables and index of behaviours
        call comp_nbvari(rela_comp, defo_comp, type_cpla, kit_comp, &
                         post_iter, mult_comp, regu_visc, &
                         extern_type, extern_addr, model_dim, &
                         nbVariUMAT, &
                         nbVari, numeLaw, nbVariKit, numeLawKit)

! ----- Save informations
        behaviourPrepPara%v_para(iFactorKeyword)%nbVari = nbVari
        behaviourPrepPara%v_para(iFactorKeyword)%nbVariKit = nbVariKit
        behaviourPrepPara%v_para(iFactorKeyword)%numeLaw = numeLaw
        behaviourPrepPara%v_para(iFactorKeyword)%numeLawKit = numeLawKit

        if (behaviourPrepPara%lDebug) then
            WRITE (6, *) "- Occurrence : ", iFactorKeyword
            WRITE (6, *) "--- nbVari : ", behaviourPrepPara%v_para(iFactorKeyword)%nbVari
            WRITE (6, *) "--- nbVariKit : ", behaviourPrepPara%v_para(iFactorKeyword)%nbVariKit
            WRITE (6, *) "--- numeLaw : ", behaviourPrepPara%v_para(iFactorKeyword)%numeLaw
            WRITE (6, *) "--- numeLawKit : ", behaviourPrepPara%v_para(iFactorKeyword)%numeLawKit
        end if

    end do
!
end subroutine
