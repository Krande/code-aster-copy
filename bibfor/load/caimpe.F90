! --------------------------------------------------------------------
! Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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
subroutine caimpe(load, mesh, nbOcc)
!
implicit none
!
#include "asterf_types.h"
#include "LoadTypes_type.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/char_crea_cart.h"
#include "asterfort/getelem.h"
#include "asterfort/getvc8.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nocart.h"
!
character(len=8), intent(in) :: load, mesh
integer, intent(in) :: nbOcc
!
! --------------------------------------------------------------------------------------------------
!
! Loads affectation - Acoustic
!
! Treatment of load IMPE_FACE
!
! --------------------------------------------------------------------------------------------------
!
! In  load             : load
! In  mesh             : mesh
! In  nbOcc            : number of factor keywords
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: keywFact = 'IMPE_FACE'
    character(len=4), parameter :: valeType = 'COMP'
    character(len=24), parameter :: listCell = '&&CAIMPE.LIST_CELL'
    integer :: nbCell, jvCell
    integer :: iocc, nbRet
    complex(kind=8) :: impe
    complex(kind=8), pointer :: valv(:) => null()
    character(len=19) :: map(LOAD_MAP_NBMAX)
    integer :: nbMap, nbCmp(LOAD_MAP_NBMAX)
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Creation and initialization to zero of <CARTE>
    call char_crea_cart('ACOUSTIQUE', keywFact, load , mesh, valeType,&
                        nbMap       , map     , nbCmp)
    ASSERT(nbMap .eq. 1)
    call jeveuo(map(1)//'.VALV', 'E', vc = valv)

! - Loop on factor keyword
    do iocc = 1, nbOcc
! ----- Read mesh affectation
        call getelem(mesh, keywFact, iocc, 'A', listCell, nbCell)

! ----- Get parameters
        call getvc8(keywFact, 'IMPE', iocc=iocc, scal = impe, nbret=nbRet)
        valv(1) = impe

! ----- Set parameter in field
        if (nbCell .ne. 0) then
            call jeveuo(listCell, 'L', jvCell)
            call nocart(map(1), 3, nbCmp(1), mode='NUM', nma=nbCell, limanu=zi(jvCell))
            call jedetr(listCell)
        endif

    end do
!
    call jedema()
end subroutine
