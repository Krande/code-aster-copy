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
subroutine comp_mfront_modelem(elem_type_name, l_mfront_cp, &
                               model_mfront, &
                               codret, type_cpla)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/teattr.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/BehaviourMGIS_type.h"
!
    character(len=16), intent(in) :: elem_type_name
    aster_logical, intent(in) :: l_mfront_cp
    integer(kind=8), intent(out) :: model_mfront
    integer(kind=8), intent(out) :: codret
    character(len=16), intent(out) :: type_cpla
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of comportment (mechanics)
!
! Select type of modelisation for MFront - On selected element
!
! --------------------------------------------------------------------------------------------------
!
! In  elem_type_name   : type of finite element
! In  l_mfront_cp      : .true. if plane stress is possible for this MFront behaviour
! Out model_mfront     : type of modelisation for MFront
! Out codret           : code for error
!                        0 - OK
!                        1 - Error - Not same finite element
!                        2 - Error - No MFront modelisation allowed on this element
! Out type_cpla        : stress plane hypothesis (for Deborst)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: iret
    character(len=16) :: principal, model_type
!
! --------------------------------------------------------------------------------------------------
!
    codret = 0
    model_mfront = MGIS_MODEL_UNSET
    type_cpla = 'VIDE'

! - Get attributes on finite element
    call teattr('C', 'TYPMOD', model_type, iret, typel=elem_type_name)
    call teattr('C', 'PRINCIPAL', principal, iret, typel=elem_type_name)

! - Select modelisation for MFront
    if (principal .eq. 'OUI') then
        if (model_type .eq. '3D') then
            model_mfront = MGIS_MODEL_TRIDIMENSIONAL
        elseif (model_type .eq. 'C_PLAN') then
            if (l_mfront_cp) then
                model_mfront = MGIS_MODEL_PLANESTRESS
                type_cpla = 'ANALYTIQUE'
            else
                model_mfront = MGIS_MODEL_AXISYMMETRICAL
                type_cpla = 'DEBORST'
            end if
        elseif (model_type .eq. 'D_PLAN') then
            model_mfront = MGIS_MODEL_PLANESTRAIN
        elseif (model_type .eq. 'PLAN') then
            model_mfront = MGIS_MODEL_PLANESTRAIN
        elseif (model_type .eq. 'AXIS') then
            model_mfront = MGIS_MODEL_AXISYMMETRICAL
        elseif (model_type .eq. '1D') then
            model_mfront = MGIS_MODEL_AXISYMMETRICAL
            type_cpla = 'DEBORST'
        else
            codret = 2
        end if
    end if
!
end subroutine
