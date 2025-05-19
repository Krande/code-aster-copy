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
! person_in_charge: nicolas.sellenet at edf.fr
!
subroutine numer3_wrap(numeDofZ, base, &
                       modelZ, listLoadZ, &
                       defiContZ, virtContElemZ, verbose)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/numer3.h"
#include "asterfort/infbav.h"
#include "asterfort/infmue.h"
!
    character(len=*), intent(inout) :: numeDofZ
    character(len=2), intent(in) :: base
    character(len=*), intent(in) :: modelZ, listLoadZ, defiContZ, virtContElemZ
    aster_logical, intent(in) :: verbose
!
! --------------------------------------------------------------------------------------------------
!
! Factor
!
! (re)-Numbering
!
! --------------------------------------------------------------------------------------------------
!
! IO  numeDof          : name of numbering object (NUME_DDL)
! In  base             : JEVEUX base to create objects
!                         base(1:1) => NUME_EQUA objects
!                         base(2:2) => NUME_DDL objects
! In  model            : name of model
! In  listLoad         : list of loads
! In  defiCont         : name of datastructure from DEFI_CONT
! In  virtContElem     : name of virtual elements for contact (in solving)
!
! --------------------------------------------------------------------------------------------------
!
    type(NL_DS_Contact) :: ds_contact
!
! --------------------------------------------------------------------------------------------------
!
    if (.not. verbose) then
        call infmue()
    end if
!
! - Prepare datastructure from contact
    if (defiContZ .ne. " ") then
        ds_contact%l_contact = ASTER_TRUE
        ds_contact%iden_rela = " "
        ds_contact%sdcont = defiContZ(1:8)
        ds_contact%ligrel_elem_cont = virtContElemZ
    end if

! - Numbering
    call numer3(modelZ, base, listLoadZ, numeDofZ, ds_contact)
!
    call infbav()
!
end subroutine
