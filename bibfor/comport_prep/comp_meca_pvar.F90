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
! person_in_charge: mickael.abbas at edf.fr
!
subroutine comp_meca_pvar(ligrel_, comporMap_, comporList_, comporInfo)
!
    use BehaviourPrepare_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/comp_meca_exc2.h"
#include "asterfort/comp_meca_l.h"
#include "asterfort/comp_meca_name.h"
#include "asterfort/comp_ntvari.h"
#include "asterfort/dismoi.h"
#include "asterfort/etenca.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jeecra.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/lteatt.h"
#include "asterfort/wkvect.h"
!
    character(len=19), optional, intent(in) :: ligrel_
    character(len=19), optional, intent(in) :: comporMap_
    character(len=16), optional, intent(in) :: comporList_(COMPOR_SIZE)
    character(len=19), intent(in) :: comporInfo
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of constitutive laws (mechanics)
!
! Prepare informations about internal state variables
!
! --------------------------------------------------------------------------------------------------
!
! In  ligrel           : ligrel
! In  comporInfo       : object for information about internal state variables and behaviour
! In  comporList       : list for parameters of constitutive laws
! In  comporMap        : map for parameters of constitutive laws
!
!    ComporInfo:
!       INFO.INFO = global parameters
!         comporInfoInfo(1) = nbCellMesh
!          => total number of elements in mesh
!         comporInfoInfo(2) = mapNbZone
!          => total number of zone in CARTE
!         comporInfoInfo(3) = nb_vari_maxi
!          => maximum number of internal variables
!         comporInfoInfo(4) = nt_vari
!          => total number of internal variables
!       INFO.VARI = Collection of mapNbZone (from CARTE) x Vecteur_Info
!       For each zone   : Vector_Info is list of nbVari name of internal variables (K16)
!       INFO.ZONE = list on mapNbZone (from CARTE)
!       For each zone   : number of elements with this comportement
!       INFO.RELA = list on mapNbZone (from CARTE) * 8
!       For each zone   : some information from comprotement (name of RELATION, DEFORMATION, ...)
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: l_excl, l_kit_meta, l_cristal, l_pmf, l_kit_thm
    aster_logical :: l_mfront_proto, l_mfront_offi
    aster_logical :: l_zone_read
    character(len=8) :: mesh
    integer(kind=8), pointer :: comporInfoInfo(:) => null()
    integer(kind=8), pointer :: comporInfoZone(:) => null()
    integer(kind=8), pointer :: zoneRead(:) => null()
    integer(kind=8), pointer :: modelCell(:) => null()
    character(len=16), pointer :: comporInfoVari(:) => null()
    character(len=16), pointer :: comporInfoRela(:) => null()
    character(len=16), pointer :: comporVale(:) => null()
    integer(kind=8), pointer :: comporDesc(:) => null()
    integer(kind=8), pointer :: comporPtma(:) => null()
    integer(kind=8) :: nbVale, mapNbCmpMax, mapNbZone, nbVari, nt_vari, nb_vari_maxi, nb_zone_acti
    integer(kind=8) :: mapZoneNume, iCellMesh, nbCellMesh, iret, nutyel, nbVariMeca, nb_zone2
    character(len=16) :: post_iter, vari_excl, regu_visc, post_incr
    character(len=16) :: rela_comp, defo_comp, kit_comp(4), type_cpla, type_comp
    character(len=255) :: libr_name, subr_name
    character(len=16) :: extern_addr, notype
    integer(kind=8) :: extern_type
    type(BehaviourPrep_Exte), pointer :: prepExte(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Initializations
    nb_zone_acti = 0

! - Access to map
    if (present(comporMap_)) then
        ASSERT(present(ligrel_))
        call jeveuo(comporMap_//'.DESC', 'L', vi=comporDesc)
        call jeveuo(comporMap_//'.VALE', 'L', vk16=comporVale)
        call jelira(comporMap_//'.VALE', 'LONMAX', nbVale)
        mapNbZone = comporDesc(3)
        mapNbCmpMax = nbVale/comporDesc(2)
        call dismoi('NOM_MAILLA', comporMap_, 'CARTE', repk=mesh)
        call dismoi('NB_MA_MAILLA', mesh, 'MAILLAGE', repi=nbCellMesh)
        call jeveuo(ligrel_//'.TYFE', 'L', vi=modelCell)
        call etenca(comporMap_, ligrel_, iret)
        call jeveuo(comporMap_//'.PTMA', 'L', vi=comporPtma)
    end if

! - Parameters if list
    if (present(comporList_)) then
        mapNbZone = 1
        mapNbCmpMax = 0
        nbCellMesh = 1
    end if

! - Create list of zones: for each zone (in CARTE), how many elements ?
    call wkvect(comporInfo(1:19)//'.ZONE', 'V V I', mapNbZone, vi=comporInfoZone)

! - Count number of cells by zone
    if (present(comporMap_)) then
        do iCellMesh = 1, nbCellMesh
            mapZoneNume = comporPtma(iCellMesh)
            if (mapZoneNume .ne. 0 .and. modelCell(iCellMesh) .ne. 0) then
                comporInfoZone(mapZoneNume) = comporInfoZone(mapZoneNume)+1
            end if
        end do
    else
        comporInfoZone(1) = 1
    end if

! - Count total of internal state variables
    if (present(comporMap_)) then
        call comp_ntvari(ligrel_=ligrel_, comporMap_=comporMap_, comporInfo=comporInfo, &
                         nt_vari=nt_vari, nb_vari_maxi=nb_vari_maxi, &
                         mapNbZone=nb_zone2, prepExte=prepExte)
    elseif (present(comporList_)) then
        call comp_ntvari(comporList_=comporList_, comporInfo=comporInfo, &
                         nt_vari=nt_vari, nb_vari_maxi=nb_vari_maxi, &
                         mapNbZone=nb_zone2, prepExte=prepExte)
    else
        ASSERT(ASTER_FALSE)
    end if
    ASSERT(nb_zone2 .eq. mapNbZone)
    AS_ALLOCATE(vi=zoneRead, size=mapNbZone)
!
! - No internal variables names
!
    if (nt_vari .eq. 0) then
        goto 99
    end if

! - Create list of comportment information (RELATION, DEFORMATION, etc.)
    call wkvect(comporInfo(1:19)//'.RELA', 'V V K16', 5*mapNbZone, vk16=comporInfoRela)
!
! - Create list of internal variables names
!
    call jecrec(comporInfo(1:19)//'.VARI', 'V V K16', 'NU', 'DISPERSE', 'VARIABLE', mapNbZone)
    do mapZoneNume = 1, mapNbZone
        call jecroc(jexnum(comporInfo(1:19)//'.VARI', mapZoneNume))
    end do
!
    do iCellMesh = 1, nbCellMesh
! ----- Get current zone
        if (present(comporMap_)) then
            mapZoneNume = comporPtma(iCellMesh)
            if (mapZoneNume .eq. 0) then
                l_zone_read = ASTER_TRUE
            else
                ASSERT(mapZoneNume .ne. 0)
                l_zone_read = zoneRead(mapZoneNume) .eq. 1
            end if
        else
            mapZoneNume = 1
            l_zone_read = ASTER_FALSE
        end if
        if (.not. l_zone_read) then
! --------- Get parameters
            if (present(comporMap_)) then
                rela_comp = comporVale(mapNbCmpMax*(mapZoneNume-1)+RELA_NAME)
                defo_comp = comporVale(mapNbCmpMax*(mapZoneNume-1)+DEFO)
                type_comp = comporVale(mapNbCmpMax*(mapZoneNume-1)+INCRELAS)
                type_cpla = comporVale(mapNbCmpMax*(mapZoneNume-1)+PLANESTRESS)
                kit_comp(1) = comporVale(mapNbCmpMax*(mapZoneNume-1)+KIT1_NAME)
                kit_comp(2) = comporVale(mapNbCmpMax*(mapZoneNume-1)+KIT2_NAME)
                kit_comp(3) = comporVale(mapNbCmpMax*(mapZoneNume-1)+KIT3_NAME)
                kit_comp(4) = comporVale(mapNbCmpMax*(mapZoneNume-1)+KIT4_NAME)
                post_iter = comporVale(mapNbCmpMax*(mapZoneNume-1)+POSTITER)
                read (comporVale(mapNbCmpMax*(mapZoneNume-1)+NVAR), '(I16)') nbVari
                nbVariMeca = 0
                if (comporVale(mapNbCmpMax*(mapZoneNume-1)+MECA_NVAR) .ne. 'VIDE') then
                    read (comporVale(mapNbCmpMax*(mapZoneNume-1)+MECA_NVAR), '(I16)') nbVariMeca
                end if
                regu_visc = comporVale(mapNbCmpMax*(mapZoneNume-1)+REGUVISC)
                post_incr = comporVale(mapNbCmpMax*(mapZoneNume-1)+POSTINCR)
            else
                rela_comp = comporList_(RELA_NAME)
                defo_comp = comporList_(DEFO)
                type_cpla = comporList_(PLANESTRESS)
                type_comp = comporList_(INCRELAS)
                kit_comp(1) = comporList_(KIT1_NAME)
                kit_comp(2) = comporList_(KIT2_NAME)
                kit_comp(3) = comporList_(KIT3_NAME)
                kit_comp(4) = comporList_(KIT4_NAME)
                post_iter = comporList_(POSTITER)
                read (comporList_(NVAR), '(I16)') nbVari
                nbVariMeca = 0
                if (comporList_(MECA_NVAR) .ne. 'VIDE') then
                    read (comporList_(MECA_NVAR), '(I16)') nbVariMeca
                end if
                regu_visc = comporList_(REGUVISC)
                post_incr = comporList_(POSTINCR)
            end if
! --------- Detection of specific cases
            call comp_meca_l(rela_comp, 'KIT_THM', l_kit_thm)
            call comp_meca_l(rela_comp, 'KIT_META', l_kit_meta)
            call comp_meca_l(rela_comp, 'CRISTAL', l_cristal)
            if (present(comporList_)) then
                l_pmf = ASTER_FALSE
            else
                nutyel = modelCell(iCellMesh)
                if (nutyel .eq. 0) then
                    l_pmf = ASTER_FALSE
                else
                    call jenuno(jexnum('&CATA.TE.NOMTE', nutyel), notype)
                    l_pmf = lteatt('TYPMOD2', 'PMF', typel=notype)
                end if
            end if
! --------- Parameters for external constitutive laws
            l_mfront_proto = prepExte(mapZoneNume)%l_mfront_proto
            l_mfront_offi = prepExte(mapZoneNume)%l_mfront_offi
            subr_name = prepExte(mapZoneNume)%subr_name
            libr_name = prepExte(mapZoneNume)%libr_name
            extern_addr = prepExte(mapZoneNume)%extern_addr
            extern_type = prepExte(mapZoneNume)%extern_type

! --------- Exception for name of internal state variables
            call comp_meca_exc2(l_cristal, l_pmf, &
                                l_excl, vari_excl)

! --------- Save names of relation
            comporInfoRela(5*(mapZoneNume-1)+1) = rela_comp
            comporInfoRela(5*(mapZoneNume-1)+2) = defo_comp
            comporInfoRela(5*(mapZoneNume-1)+3) = type_cpla
            comporInfoRela(5*(mapZoneNume-1)+4) = regu_visc
            comporInfoRela(5*(mapZoneNume-1)+5) = post_incr

! --------- Get names of internal state variables
            call jeecra(jexnum(comporInfo(1:19)//'.VARI', mapZoneNume), 'LONMAX', nbVari)
            call jeveuo(jexnum(comporInfo(1:19)//'.VARI', mapZoneNume), 'E', vk16=comporInfoVari)
            call comp_meca_name(nbVari, nbVariMeca, l_excl, vari_excl, l_kit_meta, &
                                rela_comp, defo_comp, kit_comp, type_cpla, post_iter, &
                                regu_visc, post_incr, &
                                extern_addr, extern_type, comporInfoVari)

! --------- Save current zone
            zoneRead(mapZoneNume) = 1
            nb_zone_acti = nb_zone_acti+1
        end if
    end do
!
99  continue
!
! - Save general information
!
    call wkvect(comporInfo(1:19)//'.INFO', 'V V I', 5, vi=comporInfoInfo)
    comporInfoInfo(1) = nbCellMesh
    comporInfoInfo(2) = mapNbZone
    comporInfoInfo(3) = nb_vari_maxi
    comporInfoInfo(4) = nt_vari
    comporInfoInfo(5) = nb_zone_acti
!
    deallocate (prepExte)
    AS_DEALLOCATE(vi=zoneRead)
!
    call jedema()
!
end subroutine
