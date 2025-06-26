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
subroutine op0028()
!
    implicit none
!
#include "asterc/getres.h"
#include "asterfort/dfllad.h"
#include "asterfort/dflldb.h"
#include "asterfort/dfllec.h"
#include "asterfort/dfllty.h"
#include "asterfort/getvid.h"
#include "asterfort/infmaj.h"
#include "asterfort/infniv.h"
#include "asterfort/wkvect.h"
!
! --------------------------------------------------------------------------------------------------
!
! DEFI_LIST_INST
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8) :: sdlist
    character(len=16) :: k16bid
    character(len=16) :: list_method
    real(kind=8) :: dtmin
    integer(kind=8) :: ifm, niv
    integer(kind=8):: nocc
    character(len=8)::model
    character(len=8), pointer:: sdlistModel(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call infmaj()
    call infniv(ifm, niv)

! - Get result datastructure
    call getres(sdlist, k16bid, k16bid)

! - Read model if present
    call getvid(' ', 'MODELE', scal=model, nbret=nocc)
    if (nocc .eq. 0) then
        model = ' '
    end if
    call wkvect(sdlist//'.MODELE', 'G V K8', 1, vk8=sdlistModel)
    sdlistModel(1) = model

! - Read parameters for keyword DEFI_LIST
    call dfllty(sdlist, list_method, dtmin)

! - Read parameters for keyword ECHEC
    call dfllec(sdlist, dtmin)

! - Read parameters for keyword ADAPTATION
    if (list_method .eq. 'AUTO') then
        call dfllad(sdlist)
    end if

! - Print debug
    if (niv .ge. 1) then
        call dflldb(sdlist)
    end if

end subroutine
