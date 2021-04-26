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
subroutine caimpe(phenom, load, mesh, valeType, nbOcc)
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
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nocart.h"
!
character(len=16), intent(in) :: phenom
character(len=8), intent(in) :: load, mesh
character(len=4), intent(in) :: valeType
integer, intent(in) :: nbOcc
!
! --------------------------------------------------------------------------------------------------
!
! Loads affectation - Acoustic and mechanic
!
! Treatment of load IMPE_FACE
!
! --------------------------------------------------------------------------------------------------
!
! In  phenom           : phenomenon (MECANIQUE/THERMIQUE/ACOUSTIQUE)
! In  load             : load
! In  mesh             : mesh
! In  valeType         : affected value type (real, complex or function)
! In  nbOcc            : number of factor keywords
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: keywFact = 'IMPE_FACE'
    character(len=24), parameter :: listCell = '&&CAIMPE.LIST_CELL'
    integer :: nbCell, jvCell
    integer :: iocc, nbRet
    complex(kind=8) :: impeCplx
    real(kind=8) :: impeReal
    character(len=8) :: impeFunc
    integer :: jvValv
    character(len=19) :: map(LOAD_MAP_NBMAX)
    integer :: nbMap, nbCmp(LOAD_MAP_NBMAX)
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Creation and initialization to zero of <CARTE>
    call char_crea_cart(phenom, keywFact, load , mesh, valeType,&
                        nbMap , map     , nbCmp)
    ASSERT(nbMap .eq. 1)
    call jeveuo(map(1)//'.VALV', 'E', jvValv)

! - Loop on factor keyword
    do iocc = 1, nbOcc
! ----- Read mesh affectation
        call getelem(mesh, keywFact, iocc, 'A', listCell, nbCell)

! ----- Get parameter
        if (valeType .eq. 'REEL') then
            call getvr8(keywFact, 'IMPE', iocc=iocc, scal = impeReal, nbret=nbRet)
            zr(jvValv-1+1) = impeReal
        elseif (valeType .eq. 'COMP') then
            call getvc8(keywFact, 'IMPE', iocc=iocc, scal = impeCplx, nbret=nbRet)
            zc(jvValv-1+1) = impeCplx
        elseif (valeType .eq. 'FONC') then
            call getvid(keywFact, 'IMPE', iocc=iocc, scal = impeFunc, nbret=nbRet)
            zk8(jvValv-1+1) = impeFunc
        else
            ASSERT(ASTER_FALSE)
        endif

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
