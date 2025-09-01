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
subroutine comp_meta_pvar(model, comporMeta, comporMetaInfo)
!
    use Metallurgy_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterc/lccree.h"
#include "asterc/lcdiscard.h"
#include "asterc/lcinfo.h"
#include "asterc/lcvari.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/etenca.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jeecra.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/wkvect.h"
!
    character(len=8), intent(in) :: model
    character(len=24), intent(in) :: comporMeta
    character(len=19), intent(in) :: comporMetaInfo
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of comportment (metallurgy)
!
! Prepare informations about internal variables
!
! --------------------------------------------------------------------------------------------------
!
! In  model             : name of model
! In  comporMetaInfo    : name of object for information about internal variables and comportement
! In  comporMeta        : name of map for behaviour in metallurgy
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: l_zone_read
    character(len=8) :: mesh
    character(len=19) :: modelLigrel
    integer(kind=8), pointer :: comporInfoInfo(:) => null()
    integer(kind=8), pointer :: comporInfoZone(:) => null()
    integer(kind=8), pointer :: zoneRead(:) => null()
    integer(kind=8), pointer :: modelCell(:) => null()
    character(len=16), pointer :: comporInfoVari(:) => null()
    character(len=16), pointer :: comporInfoRela(:) => null()
    character(len=16), pointer :: comporMetaVale(:) => null()
    integer(kind=8), pointer :: comporMetaDesc(:) => null()
    integer(kind=8), pointer :: comporMetaPtma(:) => null()
    integer(kind=8) :: iret, iDummy, idummy2
    integer(kind=8) :: nbVale, mapNbCmpMax, mapNbZone, nbZoneActi
    integer(kind=8) :: mapZoneNume, iCellMesh, nbCellMesh, numeComp
    integer(kind=8) :: nbVari, ntVari, nbVariMaxi, nbPhase
    character(len=16) :: metaType, metaLaw
    integer(kind=8) :: nbCompElem
    character(len=16) :: compElem(2), compCodePY, metaCodePY
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Initializations
    nbZoneActi = 0

! - Access to map
    call jeveuo(comporMeta(1:19)//'.DESC', 'L', vi=comporMetaDesc)
    call jeveuo(comporMeta(1:19)//'.VALE', 'L', vk16=comporMetaVale)
    call jelira(comporMeta(1:19)//'.VALE', 'LONMAX', nbVale)
    mapNbZone = comporMetaDesc(3)
    mapNbCmpMax = nbVale/comporMetaDesc(2)
    call dismoi('NOM_MAILLA', comporMeta, 'CARTE', repk=mesh)
    call dismoi('NB_MA_MAILLA', mesh, 'MAILLAGE', repi=nbCellMesh)
    call dismoi('NOM_LIGREL', model, 'MODELE', repk=modelLigrel)
    call jeveuo(modelLigrel//'.TYFE', 'L', vi=modelCell)
    call etenca(comporMeta, modelLigrel, iret)
    call jeveuo(comporMeta(1:19)//'.PTMA', 'L', vi=comporMetaPtma)

! - Create list of zones: for each zone (in CARTE), how many elements ?
    call wkvect(comporMetaInfo(1:19)//'.ZONE', 'V V I', mapNbZone, vi=comporInfoZone)

! - Count number of elements by zone (in CARTE)
    do iCellMesh = 1, nbCellMesh
        mapZoneNume = comporMetaPtma(iCellMesh)
        if (mapZoneNume .ne. 0 .and. modelCell(iCellMesh) .ne. 0) then
            comporInfoZone(mapZoneNume) = comporInfoZone(mapZoneNume)+1
        end if
    end do

! - Count total of internal variables
    ntVari = 0
    nbVariMaxi = 0
    do mapZoneNume = 1, mapNbZone
        read (comporMetaVale(mapNbCmpMax*(mapZoneNume-1)+2), '(I16)') nbVari
        ntVari = ntVari+nbVari
        nbVariMaxi = max(nbVariMaxi, nbVari)
    end do
    AS_ALLOCATE(vi=zoneRead, size=mapNbZone)
    if (ntVari .eq. 0) then
        goto 99
    end if

! - Create list of comportment information
    call wkvect(comporMetaInfo(1:19)//'.RELA', 'V V K16', 3*mapNbZone, vk16=comporInfoRela)

! - Create list of internal variables names
    call jecrec(comporMetaInfo(1:19)//'.VARI', 'V V K16', 'NU', 'DISPERSE', 'VARIABLE', mapNbZone)
    do mapZoneNume = 1, mapNbZone
        call jecroc(jexnum(comporMetaInfo(1:19)//'.VARI', mapZoneNume))
    end do
!
    do iCellMesh = 1, nbCellMesh
! ----- Get current zone
        mapZoneNume = comporMetaPtma(iCellMesh)
        if (mapZoneNume .eq. 0) then
            l_zone_read = ASTER_TRUE
        else
            ASSERT(mapZoneNume .ne. 0)
            l_zone_read = zoneRead(mapZoneNume) .eq. 1
        end if
        if (.not. l_zone_read) then
! --------- Get parameters
            metaType = comporMetaVale(mapNbCmpMax*(mapZoneNume-1)+1)
            metaLaw = comporMetaVale(mapNbCmpMax*(mapZoneNume-1)+3)
            read (comporMetaVale(mapNbCmpMax*(mapZoneNume-1)+2), '(I16)') nbVari
            read (comporMetaVale(mapNbCmpMax*(mapZoneNume-1)+4), '(I16)') numeComp

! --------- Create composite comportment
            nbCompElem = 2
            compElem(1) = metaType
            compElem(2) = metaLaw
            call lccree(nbCompElem, compElem, compCodePY)

! --------- Get number of phases
            nbCompElem = 1
            compElem(1) = metaType
            call lccree(nbCompElem, compElem, metaCodePY)

! --------- Get number of variables and index of behaviour
            call lcinfo(metaCodePY, idummy, nbPhase, idummy2)

! --------- Save names of relation
            comporInfoRela(3*(mapZoneNume-1)+1) = metaType
            comporInfoRela(3*(mapZoneNume-1)+2) = metaLaw
            write (comporInfoRela(3*(mapZoneNume-1)+3), '(I16)') nbPhase

! --------- Save name of internal variables
            call jeecra(jexnum(comporMetaInfo(1:19)//'.VARI', mapZoneNume), 'LONMAX', nbVari)
            call jeveuo(jexnum(comporMetaInfo(1:19)//'.VARI', mapZoneNume), &
                        'E', vk16=comporInfoVari)
            call lcvari(compCodePY, nbVari, comporInfoVari)
            call lcdiscard(compCodePY)
            call lcdiscard(metaCodePY)

! --------- Save current zone
            zoneRead(mapZoneNume) = 1
            nbZoneActi = nbZoneActi+1
        end if
    end do
!
99  continue

! - Save general information
    call wkvect(comporMetaInfo(1:19)//'.INFO', 'V V I', 5, vi=comporInfoInfo)
    comporInfoInfo(1) = nbCellMesh
    comporInfoInfo(2) = mapNbZone
    comporInfoInfo(3) = nbVariMaxi
    comporInfoInfo(4) = ntVari
    comporInfoInfo(5) = nbZoneActi
!
    AS_DEALLOCATE(vi=zoneRead)
!
    call jedema()
!
end subroutine
