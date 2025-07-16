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
subroutine caethr(load, mesh, model, valeType)
!
    implicit none
!
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterfort/char_crea_cart.h"
#include "asterfort/assert.h"
#include "asterfort/exixfe.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nocart.h"
#include "asterfort/getelem.h"
#include "asterfort/utmess.h"
!
    character(len=8), intent(in) :: load, mesh, model
    character(len=4), intent(in) :: valeType
!
! --------------------------------------------------------------------------------------------------
!
! Loads affectation
!
! Treatment of load ECHANGE_THM
!
! --------------------------------------------------------------------------------------------------
!
! In  load      : load
! In  mesh      : mesh
! In  model     : model
! In  valeType  : affected value type (real, complex or function)
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: keywordfact = 'ECHANGE_THM_HR'
    character(len=24), parameter :: listCell = '&&CAETHR.LIST_ELEM'
    integer(kind=8) :: jnfis, jvalv, jvCell
    integer(kind=8) :: nbCell, nbOcc(3), nfiss, nech
    integer(kind=8) :: iret, iocc
    character(len=19) :: map(LOAD_MAP_NBMAX)
    integer(kind=8) :: nbMap, nbCmp(LOAD_MAP_NBMAX)
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
! - Initializations
!
    nbocc(:) = 0
    call exixfe(model, iret)
    nfiss = 0
    if (iret .ne. 0) then
        call jeveuo(model//'.NFIS', 'L', jnfis)
        nfiss = zi(jnfis)
    end if
    call getfac(keywordFact, nech)
    if (nech .eq. 0) goto 99
    if (nfiss .ne. 0) then
        call utmess('F', 'XFEM_48')
    end if
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
    do iocc = 1, nech
! ----- Read mesh affectation
        call getelem(mesh, keywordfact, iocc, 'A', listCell, nbCell)

        if (nbCell .ne. 0) then
            if (valeType .eq. 'REEL') then
                call getvr8(keywordFact, 'ALPHA', iocc=iocc, scal=zr(jvalv), nbret=nbOcc(1))
                call getvr8(keywordFact, 'PVAP_SAT', iocc=iocc, scal=zr(jvalv+1), nbret=nbOcc(2))
                call getvr8(keywordFact, 'HR_EXT', iocc=iocc, scal=zr(jvalv+2), nbret=nbOcc(3))
            elseif (valeType .eq. 'FONC') then
                call getvid(keywordFact, 'ALPHA', iocc=iocc, scal=zk8(jvalv), nbret=nbOcc(1))
                call getvid(keywordFact, 'PVAP_SAT', iocc=iocc, scal=zk8(jvalv+1), nbret=nbOcc(2))
                call getvid(keywordFact, 'HR_EXT', iocc=iocc, scal=zk8(jvalv+2), nbret=nbOcc(3))
            else
                ASSERT(ASTER_FALSE)
            end if
            call jeveuo(listCell, 'L', jvCell)
            call nocart(map(1), 3, nbCmp(1), mode='NUM', nma=nbCell, &
                        limanu=zi(jvCell))
            call jedetr(listCell)
        end if
    end do
99  continue
!
    call jedema()
end subroutine
