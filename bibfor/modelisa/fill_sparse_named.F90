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
!
subroutine fill_sparse_named(coll, names, sizes, groupsOfCells)
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jeecra.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
!
    character(len=*), intent(in) :: coll, names, sizes, groupsOfCells
!
    integer :: nbGroups, ichk, i, j, idx
    character(len=24), pointer :: vnam(:) => null()
    integer, pointer :: vsiz(:) => null(), vgrp(:) => null(), vnum(:) => null()

    call jemarq()
    call jelira(names, "LONUTI", ival=nbGroups)
    call jelira(sizes, "LONUTI", ival=ichk)
    ASSERT(nbGroups .eq. ichk)
    call jeveuo(names, "L", vk24=vnam)
    call jeveuo(sizes, "L", vi=vsiz)
    call jeveuo(groupsOfCells, "L", vi=vgrp)

    idx = 1
    do i = 1, nbGroups
        ! print *, vnam(i), ":", vsiz(i)
        call jecroc(jexnom(coll, vnam(i)))
        call jeecra(jexnom(coll, vnam(i)), "LONMAX", vsiz(i))
        call jeecra(jexnom(coll, vnam(i)), "LONUTI", vsiz(i))
        call jeveuo(jexnom(coll, vnam(i)), "E", vi=vnum)
        do j = 1, vsiz(i)
            vnum(j) = vgrp(idx+j-1)
        end do
        idx = idx+vsiz(i)
    end do

    call jedema()
end subroutine fill_sparse_named
