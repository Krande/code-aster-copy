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
subroutine caonde(load, mesh, valeType, nbOcc)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/getvr8.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nocart.h"
#include "asterfort/char_crea_cart.h"
#include "asterfort/getelem.h"
!
    character(len=8), intent(in) :: load, mesh
    character(len=4), intent(in) :: valeType
    integer(kind=8), intent(in) :: nbOcc
!
! --------------------------------------------------------------------------------------------------
!
! Loads affectation
!
! Treatment of load ONDE_FLUI
!
! --------------------------------------------------------------------------------------------------
!
! In  load      : load
! In  mesh      : mesh
! In  valeType  : affected value type (real, complex or function)
! In  nbOcc     : number of factor keywords
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: keywordfact = 'ONDE_FLUI'
    character(len=24), parameter :: listCell = '&&CAONDE.LIST_ELEM'
    integer(kind=8) :: iocc, n, jvalv, nbCell, jvCell
    character(len=19) :: map(LOAD_MAP_NBMAX)
    integer(kind=8) :: nbMap, nbCmp(LOAD_MAP_NBMAX)
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
! - Initializations
!
    ASSERT(valeType .eq. 'REEL')
!
! - Creation and initialization to zero of <CARTE>
!
    call char_crea_cart('MECANIQUE', keywordfact, load, mesh, valeType, &
                        nbMap, map, nbCmp)
    ASSERT(nbMap .eq. 1)
    call jeveuo(map(1)//'.VALV', 'E', jvalv)
!
! - Loop on factor keyword
!
    do iocc = 1, nbOcc
! ----- Read mesh affectation
        call getelem(mesh, keywordfact, iocc, 'A', listCell, nbCell)
        if (nbCell .ne. 0) then
            call getvr8(keywordfact, 'PRES', iocc=iocc, scal=zr(jvalv), nbret=n)
            call jeveuo(listCell, 'L', jvCell)
            call nocart(map(1), 3, nbCmp(1), mode='NUM', nma=nbCell, &
                        limanu=zi(jvCell))
            call jedetr(listCell)
        end if
    end do
!
    call jedema()
end subroutine
