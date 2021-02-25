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

subroutine cafotu(load, ligrmo, mapAlreadyCreated, mesh, ndim, valeType, nbOcc)
!
implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/char_affe_neum.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/char_crea_cart.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
!
aster_logical, intent(in) :: mapAlreadyCreated
character(len=8), intent(in) :: load, mesh
character(len=19), intent(in) :: ligrmo
character(len=4), intent(in) :: valeType
integer, intent(in) :: nbOcc, ndim
!
! --------------------------------------------------------------------------------------------------
!
! Loads affectation
!
! Treatment of load FORCE_TUYAU
!
! --------------------------------------------------------------------------------------------------
!
! In  load             : load
! In  ligrmo           : list of elements in model
! In  mapAlreadyCreated : flag when maps already created
! In  mesh             : mesh
! In  ndim             : space dimension
! In  valeType         : affected value type (real, complex or function)
! In  nbOcc            : number of factor keywords
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: keywordfact = 'FORCE_TUYAU'
    integer :: jvalv, iocc, nbret
    character(len=8) :: model
    character(len=19) :: map(LOAD_MAP_NBMAX)
    integer :: nbMap, nbCmp(LOAD_MAP_NBMAX)
    aster_logical :: createMap
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
! - Does create map ?
!
    createMap = ASTER_TRUE
    if (mapAlreadyCreated) then
        createMap = ASTER_FALSE
    endif
!
! - Initializations
!
    call dismoi('NOM_MODELE', ligrmo, 'LIGREL', repk=model)
    ASSERT(valeType .eq. 'REEL' .or. valeType .eq. 'FONC')
!
! - Creation and initialization to zero of <CARTE>
!
    call char_crea_cart('MECANIQUE', keywordfact, load , mesh, valeType,&
                        nbMap      , map        , nbCmp, createMap)
    ASSERT(nbMap .eq. 1)
    call jeveuo(map(1)//'.VALV', 'E', jvalv)
!
! - Loop on factor keyword
!
    do iocc = 1, nbOcc

! ----- Get values of load
        if (valeType .eq. 'REEL') then
            call getvr8(keywordfact, 'PRES', iocc=iocc, scal=zr(jvalv), nbret=nbret)
        else
            call getvid(keywordfact, 'PRES', iocc=iocc, scal=zk8(jvalv), nbret=nbret)
        endif

! ----- Affect values of load
        call char_affe_neum(model      , mesh, ndim,&
                            keywordfact, iocc,&
                            nbMap      , map , nbCmp)

    end do
!
    call jedema()
end subroutine
