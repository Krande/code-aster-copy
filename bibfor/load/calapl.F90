! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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
subroutine calapl(load, mesh, model, nbOcc)
!
implicit none
!
#include "jeveux.h"
#include "LoadTypes_type.h"
#include "asterfort/assert.h"
#include "asterfort/alcart.h"
#include "asterfort/codent.h"
#include "asterfort/char_crea_cart.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvr8.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/nocart.h"
#include "asterfort/getelem.h"
#include "asterfort/wkvect.h"
!
character(len=8), intent(in) :: load, mesh, model
integer, intent(in) :: nbOcc
!
! --------------------------------------------------------------------------------------------------
!
! Loads affectation
!
! Keyword = 'INTE_ELEC'
!
! --------------------------------------------------------------------------------------------------
!
! In  load             : load
! In  model            : model
! In  mesh             : mesh
! In  nbOcc            : number of factor keywords
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: keywordfact = 'INTE_ELEC'
    character(len=24), parameter :: listCell1 = '&&CALAPL.LISTEL1'
    character(len=24), parameter :: listCell2 = '&&CALAPL.LISTEL2'
    character(len=4), parameter :: valeType = 'REEL'
    integer :: iCell2, iocc
    integer :: nbCell1, nbCell2, ntra, nsym
    integer :: jma1, jma2, jtran, jno, jnuma
    integer :: cellNume
    character(len=16) :: listma, ltrans
    character(len=24) :: connex
    character(len=16), pointer :: valv(:) => null()
    character(len=8), pointer :: ncmp(:) => null()
    character(len=8) :: physQuantity(LOAD_MAP_NBMAX), cmpName(LOAD_MAP_NBMAX, LOAD_MAP_NBCMPMAX)
    character(len=19) :: map(LOAD_MAP_NBMAX)
    integer :: nbMap, nbCmp(LOAD_MAP_NBMAX)
    aster_logical, parameter :: createMap = ASTER_FALSE
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Initializations
    ASSERT(nbOcc .le. 99)
!
! - Creation and initialization to zero of <CARTE>
!
    call char_crea_cart('MECANIQUE', keywordfact , load , mesh, valeType,&
                        nbMap      , map         , nbCmp,&
                        createMap  , physQuantity, cmpName)
    ASSERT(nbMap .eq. 1)
    ASSERT(nbCmp(1) .eq. 2)
!
! - Access to mesh
!
    connex = mesh//'.CONNEX'
!
! - Loop on factor keyword
!
    do iocc = 1, nbOcc

! ----- Generate names
        listma(1:14) = load//'.LISMA'
        ltrans(1:14) = load//'.TRANS'
        call codent(iocc, 'D0', map(1)(18:19))
        call codent(iocc, 'D0', listma(15:16))
        call codent(iocc, 'D0', ltrans(15:16))

! ----- Create map
        call alcart('G', map(1), mesh, physQuantity(1))
        call jeveuo(map(1)//'.NCMP', 'E', vk8=ncmp)
        call jeveuo(map(1)//'.VALV', 'E', vk16=valv)

! ----- Zero everywhere
        ncmp(1) = cmpName(1, 1)
        ncmp(2) = cmpName(1, 2)
        valv(1) = ' '
        valv(2) = ' '
        call nocart(map(1), 1, nbCmp(1))

! ----- Create object
        call wkvect(ltrans, 'G V R', 6, jtran)

! ----- Set values in map
        ncmp(1) = cmpName(1, 1)
        ncmp(2) = cmpName(1, 2)
        valv(1) = listma
        valv(2) = ltrans
!
        call getvr8(keywordfact, 'TRANS', iocc=iocc, nbval=0, nbret=ntra)
        call getvr8(keywordfact, 'SYME', iocc=iocc, nbval=0, nbret=nsym)
        ntra = -ntra
        nsym = -nsym
        if (ntra .ne. 0) then
            call getvr8(keywordfact, 'TRANS', iocc=iocc, nbval=ntra, vect=zr(jtran))
        endif
        if (nsym .ne. 0) then
            call getvr8(keywordfact, 'SYME', iocc=iocc, nbval=nsym, vect=zr(jtran))
        endif

! ----- GEOMETRIE DU CONDUCTEUR SECONDAIRE
        call getelem(mesh, keywordfact, iocc, 'A', listCell2, nbCell2, suffix = '_2', model=model)
        if (nbCell2 .eq. 0) cycle
        call jeveuo(listCell2, 'L', jma2)

! ----- GEOMETRIE DU CONDUCTEUR PRINCIPAL
        call getelem(mesh, keywordfact, iocc, 'A', listCell1, nbCell1, model=model)
        if (nbCell1 .eq. 0) cycle
        call jeveuo(listCell1, 'L', jma1)

        call wkvect(listma, 'G V I', 2*nbCell2, jnuma)
        do iCell2 = 1, nbCell2
            cellNume = zi(jma2+iCell2-1)
            call jeveuo(jexnum(connex, cellNume), 'L', jno)
            zi(jnuma+2*iCell2-2) = zi(jno-1+1)
            zi(jnuma+2*iCell2-1) = zi(jno-1+2)
        end do

! ----- Affect map
        call nocart(map(1), 3, nbCmp(1), mode='NUM', nma=nbCell1,&
                    limanu=zi(jma1))

        call jedetr(listCell1)
        call jedetr(listCell2)
!
    end do
!
    call jedema()
end subroutine
