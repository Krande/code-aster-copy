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
subroutine nume_ddl_chamElem(numeDofZ, listLigrelJvZ, modeLocZ, verbose)
!
    implicit none
!
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/infbav.h"
#include "asterfort/infmue.h"
#include "asterfort/jelira.h"
#include "asterfort/jeveuo.h"
#include "asterfort/numero.h"
!
    character(len=*), intent(in) :: numeDofZ
    character(len=*), intent(in) :: listLigrelJvZ
    character(len=*), intent(in) :: modeLocZ
    aster_logical, intent(in) :: verbose
!
! ----------------------------------------------------------------------------------------------
!
! Factor
!
! Numbering - Create numbering for a CHAM_ELNO
!
! ----------------------------------------------------------------------------------------------
!
! In  numeDof          : name of numeDof object
! In  listLigrelJv     : list of listLigrelJv
! In  modloc           : name of mode local
! In  model            : name of model
!
! ----------------------------------------------------------------------------------------------
!
    character(len=14) :: numeDof
    character(len=24), pointer :: listLigr(:) => null(), listLigrJv(:) => null()
    integer(kind=8) :: nbLigr, nbLigrJv
!
! ----------------------------------------------------------------------------------------------
!
    if (.not. verbose) then
        call infmue()
    end if
!
    numeDof = numeDofZ

! - Set list of ligrel
    call jeveuo(listLigrelJvZ, 'L', vk24=listLigrJv)
    call jelira(listLigrelJvZ, 'LONUTI', nbLigrJv)
    AS_ALLOCATE(vk24=listLigr, size=nbLigrJv)

! - Copy previous LIGREL
    listLigr(1:nbLigrJv) = listLigrJv(1:nbLigrJv)
    nbLigr = nbLigrJv

! - Numbering
    call numero(numeDof, 'GG', &
                nbLigr, listLigr, &
                modeLocZ_=modeLocZ)

    AS_DEALLOCATE(vk24=listLigr)
!
    call infbav()
!
end subroutine
