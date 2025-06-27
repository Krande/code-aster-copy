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
subroutine numero_wrap(numeDofZ, base, &
                       modelZ, listLoadZ, &
                       contactZ, verbose)
!
    use listLoad_module
    use contactAlgo_module
!
    implicit none
!
#include "asterfort/addModelLigrel.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/numero.h"
#include "asterfort/infbav.h"
#include "asterfort/infmue.h"
!
    character(len=*), intent(inout) :: numeDofZ
    character(len=2), intent(in) :: base
    character(len=*), intent(in) :: modelZ, listLoadZ, contactZ
    aster_logical, intent(in) :: verbose
!
! --------------------------------------------------------------------------------------------------
!
! Factor
!
! Numbering
!
! --------------------------------------------------------------------------------------------------
!
! IO  numeDof        : name of numbering object (NUME_DDL)
! In  base           : JEVEUX base to create objects
!                      base(1:1) => NUME_EQUA objects
!                      base(2:2) => NUME_DDL objects
! In  numeDofOld     : name of previous numeDof object
! In  modeLoc        : local mode for GRANDEUR numbering
! In  model          : name of model
! In  listLoad       : list of loads
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nbLigr
    character(len=24), pointer :: listLigr(:) => null()
    character(len=24) :: virtualCell
!
! --------------------------------------------------------------------------------------------------
!
    if (.not. verbose) then
        call infmue()
    end if
!
    nbLigr = 0

! - Add LIGREL from model
    call addModelLigrel(modelZ, nbLigr, listLigr)

! - Get list of LIGREL from loads
    call getListLoadLigrel(listLoadZ, nbLigr, listLigr)

! - Get list of LIGREL from contact
    if (contactZ .ne. " ") then
        virtualCell = "None"
        call addContactLigrel(contactZ, virtualCell, nbLigr, listLigr)
    end if

! - Numbering
    call numero(numeDofZ, base, &
                nbLigr, listLigr)
    AS_DEALLOCATE(vk24=listLigr)
!
    call infbav()
!
end subroutine
