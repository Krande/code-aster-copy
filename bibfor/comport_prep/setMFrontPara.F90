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
subroutine setMFrontPara(behaviourCrit, iFactorKeyword)
!
    use Behaviour_type
!
    implicit none
!
#include "asterc/mgis_set_double_parameter.h"
#include "asterc/mgis_set_integer_parameter.h"
#include "asterc/mgis_set_outofbounds_policy.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
!
    type(Behaviour_Crit), pointer :: behaviourCrit(:)
    integer, intent(in) :: iFactorKeyword
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of constitutive laws (mechanics)
!
! Set parameters for MFront
!
! --------------------------------------------------------------------------------------------------
!
! In  behaviourCrit    : parameters for integration of constitutive law
! In  iFactorKeyword   : index of factor keyword (for map)
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: iter_inte_maxi, resi_inte_rela
    integer:: extern_type, iveriborne
    character(len=16) :: extern_addr
!
! --------------------------------------------------------------------------------------------------
!
    iveriborne = behaviourCrit(iFactorKeyword)%iveriborne
    resi_inte_rela = behaviourCrit(iFactorKeyword)%resi_inte_rela
    iter_inte_maxi = behaviourCrit(iFactorKeyword)%iter_inte_maxi
    extern_addr = behaviourCrit(iFactorKeyword)%paraExte%extern_addr
    extern_type = behaviourCrit(iFactorKeyword)%extern_type
!
! - Set values
!
    if (extern_type .eq. 1 .or. extern_type .eq. 2) then
        call mgis_set_double_parameter(extern_addr, "epsilon", resi_inte_rela)
        call mgis_set_integer_parameter(extern_addr, "iterMax", int(iter_inte_maxi))
        call mgis_set_outofbounds_policy(extern_addr, iveriborne)
    end if
!
end subroutine
