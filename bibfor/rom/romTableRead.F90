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
subroutine romTableRead(tablReduCoor)
!
    use Rom_Datastructure_type
!
    implicit none
!
#include "asterfort/jeveuo.h"
#include "asterfort/tbexve.h"
!
    type(ROM_DS_TablReduCoor), intent(in) :: tablReduCoor
!
! --------------------------------------------------------------------------------------------------
!
! Model reduction
!
! Read reduced coordinates
!
! --------------------------------------------------------------------------------------------------
!
! In  tablReduCoor     : table for reduced coordinates
!
! --------------------------------------------------------------------------------------------------
!
    character(len=24), parameter :: coorReduObj = '&&COORHR'
!
! --------------------------------------------------------------------------------------------------
!
    if (tablReduCoor%lTablFromUser) then
        call tbexve(tablReduCoor%tablUserName, &
                    tablReduCoor%tablResu%tablSymbName, coorReduObj)
    else
        call tbexve(tablReduCoor%tablResu%tablName, &
                    tablReduCoor%tablResu%tablSymbName, coorReduObj)
    end if
    call jeveuo(coorReduObj, 'L', vr=tablReduCoor%coorRedu)
!
end subroutine
