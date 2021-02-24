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
! person_in_charge: sylvie.granet at edf.fr
!
subroutine caethm(load, mesh, ligrmo, valeType)
!
implicit none
!
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterfort/alcart.h"
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
character(len=8), intent(in) :: load, mesh
character(len=19), intent(in) :: ligrmo
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
! In  ligrmo    : model <LIGREL>
! In  valeType  : affected value type (real, complex or function)
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: keywordfact = 'ECHANGE_THM'
    character(len=24), parameter :: listCell = '&&CAETHM.LIST_ELEM'
    integer :: jnfis, jvalv, jvCell
    integer :: nbCell, nbOcc, nfiss, nech
    integer :: iret, iocc
    character(len=8) :: model
    character(len=19) :: map(LOAD_MAP_NBMAX)
    integer :: nbMap, nbCmp(LOAD_MAP_NBMAX)
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
! - Initializations
!
    model = ligrmo(1:8)
    call exixfe(ligrmo(1:8), iret)
    nfiss = 0
    if (iret .ne. 0) then
        call jeveuo(model//'.NFIS', 'L', jnfis)
        nfiss = zi(jnfis)
    endif
    call getfac(keywordFact, nech)
    if (nech .eq. 0) goto 99
    if (nfiss .ne. 0) then
       call utmess('F', 'XFEM_48')
    endif
!
! - Creation and initialization to zero of <CARTE>
!
    call char_crea_cart('MECANIQUE', keywordfact, load, mesh, valeType,&
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
                call getvr8(keywordFact, 'COEF_11', iocc=iocc, scal=zr(jvalv), nbret=nbOcc)
                call getvr8(keywordFact, 'COEF_12', iocc=iocc, scal=zr(jvalv+1), nbret=nbOcc)
                call getvr8(keywordFact, 'COEF_21', iocc=iocc, scal=zr(jvalv+2), nbret=nbOcc)
                call getvr8(keywordFact, 'COEF_22', iocc=iocc, scal=zr(jvalv+3), nbret=nbOcc)
                call getvr8(keywordFact, 'PRE1_EXT', iocc=iocc, scal=zr(jvalv+4), nbret=nbOcc)
                call getvr8(keywordFact, 'PRE2_EXT', iocc=iocc, scal=zr(jvalv+5), nbret=nbOcc)
            elseif (valeType .eq. 'FONC') then
                call getvid(keywordFact, 'COEF_11', iocc=iocc, scal=zk8(jvalv), nbret=nbOcc)
                call getvid(keywordFact, 'COEF_12', iocc=iocc, scal=zk8(jvalv+1), nbret=nbOcc)
                call getvid(keywordFact, 'COEF_21', iocc=iocc, scal=zk8(jvalv+2), nbret=nbOcc)
                call getvid(keywordFact, 'COEF_22', iocc=iocc, scal=zk8(jvalv+3), nbret=nbOcc)
                call getvid(keywordFact, 'PRE1_EXT', iocc=iocc, scal=zk8(jvalv+4), nbret=nbOcc)
                call getvid(keywordFact, 'PRE2_EXT', iocc=iocc, scal=zk8(jvalv+5), nbret=nbOcc)
            else
                ASSERT(ASTER_FALSE)
            endif
            call jeveuo(listCell, 'L', jvCell)
            call nocart(map(1), 3, nbCmp(1), mode='NUM', nma=nbCell,&
                        limanu=zi(jvCell))
            call jedetr(listCell)
        endif
    end do
 99 continue
!
    call jedema()
end subroutine
