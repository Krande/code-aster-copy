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
! person_in_charge: mickael.abbas at edf.fr
!
subroutine comp_ntvari(model_, comporMap_, comporList_, comporInfo, &
                       nt_vari, nb_vari_maxi, mapNbZone, behaviourParaExte)
!
    use Behaviour_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jelira.h"
#include "asterfort/jenuno.h"
#include "asterfort/jexnum.h"
#include "asterfort/jexatr.h"
#include "asterfort/getExternalBehaviourPara.h"
#include "asterfort/teattr.h"
#include "asterfort/Behaviour_type.h"
!
    character(len=8), optional, intent(in) :: model_
    character(len=19), optional, intent(in) :: comporMap_
    character(len=16), optional, intent(in) :: comporList_(COMPOR_SIZE)
    character(len=19), intent(in) :: comporInfo
    integer, intent(out) :: nt_vari, nb_vari_maxi, mapNbZone
    type(Behaviour_ParaExte), pointer :: behaviourParaExte(:)
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of constitutive laws (mechanics)
!
! Count total of internal state variables
!
! --------------------------------------------------------------------------------------------------
!
! In  model            : model
! In  comporList       : list for parameters of constitutive laws
! In  comporMap        : map for parameters of constitutive laws
! In  comporInfo       : object for information about internal state variables and behaviour
! Out nt_vari          : total number of internal variables (on all <CARTE> COMPOR)
! Out nb_vari_maxi     : maximum number of internal variables on all comportments"
! Out mapNbZone        : number of affected zones
! Out behaviourParaExte: pointer to external behaviours parameters
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: l_comp_external
    integer, pointer :: modelCell(:) => null()
    character(len=16), pointer :: comporVale(:) => null()
    integer, pointer :: comporInfoZone(:) => null()
    integer, pointer :: comporDesc(:) => null()
    integer, pointer :: comporLima(:) => null()
    integer, pointer :: comporLimaCumu(:) => null()
    integer :: nbVale, mapNbCmpMax, nb_vari, nbCell, nbCellMesh
    integer :: iMapZone, iret, iCell, posit
    integer :: affeZoneType, affeZoneNume, cellTypeNume, cellNume, model_dim
    character(len=16) :: cellTypeName
    character(len=16) :: post_iter
    character(len=16) :: rela_comp, defo_comp, mult_comp, kit_comp(4), type_cpla
    character(len=16) :: model_mfront
    character(len=255) :: libr_name, subr_name
    character(len=16) :: principal
    character(len=8) :: mesh
!
! --------------------------------------------------------------------------------------------------
!
    nt_vari = 0
    nb_vari_maxi = 0
    mapNbZone = 0
    behaviourParaExte => null()
    if (present(model_)) then
        call jeveuo(model_//'.MAILLE', 'L', vi=modelCell)
    end if

! - Access to map
    if (present(comporMap_)) then
        call jeveuo(comporMap_//'.DESC', 'L', vi=comporDesc)
        call jeveuo(comporMap_//'.VALE', 'L', vk16=comporVale)
        call jelira(comporMap_//'.VALE', 'LONMAX', nbVale)
        call jeveuo(jexnum(comporMap_//'.LIMA', 1), 'L', vi=comporLima)
        call jeveuo(jexatr(comporMap_//'.LIMA', 'LONCUM'), 'L', vi=comporLimaCumu)
        mapNbZone = comporDesc(3)
        mapNbCmpMax = nbVale/comporDesc(2)
        call dismoi('NOM_MAILLA', comporMap_, 'CARTE', repk=mesh)
        call dismoi('NB_MA_MAILLA', mesh, 'MAILLAGE', repi=nbCellMesh)
    end if

! - Parameters if list
    if (present(comporList_)) then
        mapNbZone = 1
        mapNbCmpMax = 0
        nbCellMesh = 1
    end if

! - Create list of zones: for each zone (in CARTE), how many elements
    call jeveuo(comporInfo(1:19)//'.ZONE', 'L', vi=comporInfoZone)

! - Prepare objects for external constitutive laws
    allocate (behaviourParaExte(mapNbZone))

! - Count internal variables by comportment
    do iMapZone = 1, mapNbZone
        subr_name = ' '
        libr_name = ' '
        model_mfront = ' '
        model_dim = 0

! ----- Get parameters
        if (present(comporMap_)) then
            rela_comp = comporVale(mapNbCmpMax*(iMapZone-1)+RELA_NAME)
            defo_comp = comporVale(mapNbCmpMax*(iMapZone-1)+DEFO)
            type_cpla = comporVale(mapNbCmpMax*(iMapZone-1)+PLANESTRESS)
            mult_comp = comporVale(mapNbCmpMax*(iMapZone-1)+MULTCOMP)
            kit_comp(1) = comporVale(mapNbCmpMax*(iMapZone-1)+KIT1_NAME)
            kit_comp(2) = comporVale(mapNbCmpMax*(iMapZone-1)+KIT2_NAME)
            kit_comp(3) = comporVale(mapNbCmpMax*(iMapZone-1)+KIT3_NAME)
            kit_comp(4) = comporVale(mapNbCmpMax*(iMapZone-1)+KIT4_NAME)
            post_iter = comporVale(mapNbCmpMax*(iMapZone-1)+POSTITER)
        else
            rela_comp = comporList_(RELA_NAME)
            defo_comp = comporList_(DEFO)
            type_cpla = comporList_(PLANESTRESS)
            mult_comp = comporList_(MULTCOMP)
            kit_comp(1) = comporList_(KIT1_NAME)
            kit_comp(2) = comporList_(KIT2_NAME)
            kit_comp(3) = comporList_(KIT3_NAME)
            kit_comp(4) = comporList_(KIT4_NAME)
            post_iter = comporList_(POSTITER)
        end if

! ----- Find right TYPELEM
        if (present(comporMap_)) then
            affeZoneType = comporDesc(1+3+(iMapZone-1)*2)
            affeZoneNume = comporDesc(1+4+(iMapZone-1)*2)
            if (affeZoneType .eq. 3) then
                nbCell = comporLimaCumu(1+affeZoneNume)-comporLimaCumu(affeZoneNume)
                posit = comporLimaCumu(affeZoneNume)
            elseif (affeZoneType .eq. 1) then
                nbCell = nbCellMesh
                posit = 0
            else
                ASSERT(ASTER_FALSE)
            end if
        else
            affeZoneType = 0
            nbCell = 1
            ASSERT(iMapZone .eq. 1)
        end if
        do iCell = 1, nbCell
            if (affeZoneType .eq. 3) then
                cellNume = comporLima(posit+iCell-1)
            elseif (affeZoneType .eq. 1) then
                cellNume = iCell
            elseif (affeZoneType .eq. 0) then
                cellNume = 1
            else
                ASSERT(.false.)
            end if
            if (cellNume .ne. 0 .and. affeZoneType .gt. 0) then
                cellTypeNume = modelCell(cellNume)
                if (cellTypeNume .ne. 0) then
                    call jenuno(jexnum('&CATA.TE.NOMTE', cellTypeNume), cellTypeName)
                    call teattr('C', 'PRINCIPAL', principal, iret, typel=cellTypeName)
                    if (principal .eq. 'OUI') then
                        goto 20
                    end if
                end if
            end if
        end do
20      continue

! ----- Get parameters for external programs (MFRONT/UMAT)
        call getExternalBehaviourPara(mesh, modelCell, &
                                      rela_comp, kit_comp, &
                                      l_comp_external, behaviourParaExte(iMapZone), &
                                      elem_type_=cellTypeNume, &
                                      type_cpla_in_=type_cpla)

! ----- Get number of internal variables
        if (present(comporMap_)) then
            read (comporVale(mapNbCmpMax*(iMapZone-1)+NVAR), '(I16)') nb_vari
        else
            read (comporList_(NVAR), '(I16)') nb_vari
        end if
        nt_vari = nt_vari+nb_vari
        nb_vari_maxi = max(nb_vari_maxi, nb_vari)
    end do
!
end subroutine
