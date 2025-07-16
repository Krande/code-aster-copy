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
!
subroutine nmetcv(model, fieldRefe, fieldType, &
                  fieldInName, fieldInDisc, fieldOutName, fieldOutDisc)
!
    implicit none
!
#include "asterfort/chpchd.h"
#include "asterfort/copisd.h"
#include "asterfort/utmess.h"
!
    character(len=8), intent(in) :: model
    character(len=24), intent(in) :: fieldRefe
    character(len=24), intent(in) :: fieldType
    character(len=24), intent(in) :: fieldInName, fieldOutName
    character(len=4), intent(in) :: fieldInDisc, fieldOutDisc
!
! --------------------------------------------------------------------------------------------------
!
! *_NON_LINE - Input/output datastructure
!
! Field conversion (discretization)
!
! --------------------------------------------------------------------------------------------------
!
! In  model           : model
! In  fieldtype       : type of field
! In  fieldRefe       : name of a reference field to convert ELGA fields
! In  fieldInName     : name of field to convert
! In  fieldInDisc     : spatial discretization of field to convert
! In  fieldOutName    : name of field converted
! In  fieldOutDisc    : spatial discretization of field converted
!
! --------------------------------------------------------------------------------------------------
!
    character(len=24) :: valk(3)
!
! --------------------------------------------------------------------------------------------------
!

    if (fieldInDisc .eq. fieldOutDisc) then
! ----- Good discretization -> nothing to do
        call copisd('CHAMP_GD', 'V', fieldInName, fieldOutName)
    else
! ----- Not good discretization -> is it possible to convert ?
        valk(1) = fieldType
        valk(2) = fieldInDisc
        valk(3) = fieldOutDisc
        if (fieldOutDisc .eq. 'ELGA') then
            if (fieldRefe .eq. ' ') then
                call utmess('F', 'ETATINIT_52', nk=3, valk=valk)
            else
                call utmess('I', 'ETATINIT_51', nk=3, valk=valk)
            end if
        else
            call utmess('F', 'ETATINIT_52', nk=3, valk=valk)
        end if

! ----- Not good discretization -> try to convert
        call chpchd(fieldInName, fieldOutDisc, fieldRefe, 'NON', 'V', &
                    fieldOutName, model)
    end if
!
end subroutine
