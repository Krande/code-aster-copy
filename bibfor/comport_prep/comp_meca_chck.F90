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
subroutine comp_meca_chck(model, mesh, chmate, &
                          fullElemField, lInitialState, behaviourPrepPara)
!
    use Behaviour_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterc/lccree.h"
#include "asterc/lcdiscard.h"
#include "asterc/lctest.h"
#include "asterfort/assert.h"
#include "asterfort/compMecaChckModel.h"
#include "asterfort/compMecaChckStrain.h"
#include "asterfort/compMecaSelectPlaneStressAlgo.h"
#include "asterfort/comp_read_mesh.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvtx.h"
#include "asterfort/utmess.h"
#include "asterfort/nmvcd2.h"
!
#include "asterc/asmpi_comm.h"
#include "asterfort/asmpi_info.h"
!
    character(len=8), intent(in) :: model, mesh, chmate
    character(len=19), intent(in) :: fullElemField
    aster_logical, intent(in) :: lInitialState
    type(Behaviour_PrepPara), intent(inout) :: behaviourPrepPara
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of behaviour (mechanics)
!
! Check with Comportement.py
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : name of mesh
! In  model            : name of model
! In  chmate           : material field
! In  fullElemField    : <CHELEM_S> of FULL_MECA option
! In  lInitialState    : .true. if initial state is defined
! IO  behaviourPrepPara: datastructure to prepare behaviour
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: factorKeyword = 'COMPORTEMENT'
    character(len=24), parameter :: cellAffe = '&&COMPMECASAVE.LIST'
    aster_logical :: lAllCellAffe
    integer :: nbCellAffe
    integer :: iFactorKeyword, nbFactorKeyword, exteDefo, lctestIret
    character(len=16) :: defoComp, relaComp, typeCpla, typeComp, reguVisc, postIncr
    character(len=16) :: relaCompPY, defoCompPY
    character(len=19) :: partit, answer
    character(len=24) :: modelLigrel
    mpi_int :: nbCPU, mpiCurr
    aster_logical :: lElasByDefault, lNeedDeborst, lMfront, lDistParallel
    aster_logical :: lIncoUpo, lExistVarc,exis_temp,exis_sech
!
! --------------------------------------------------------------------------------------------------
!
    nbFactorKeyword = behaviourPrepPara%nb_comp
    lNeedDeborst = ASTER_FALSE
    lElasByDefault = ASTER_FALSE
    lDistParallel = ASTER_FALSE

! - MPI initialisation
    call asmpi_comm('GET', mpiCurr)
    call asmpi_info(mpiCurr, size=nbCPU)

! - Generic properties
    call dismoi('EXI_VARC', chmate, 'CHAM_MATER', repk=answer)
    lExistVarc = answer .eq. 'OUI'
    call dismoi('NOM_LIGREL', model, 'MODELE', repk=modelLigrel)
    call dismoi('PARTITION', modelLigrel, 'LIGREL', repk=partit)

! - Distributed parallelism
    lDistParallel = partit .ne. ' ' .and. nbCPU .gt. 1

! - Loop on occurrences of COMPORTEMENT
    do iFactorKeyword = 1, nbFactorKeyword

! ----- Get list of cells where behaviour is defined
        call comp_read_mesh(mesh, factorKeyword, iFactorKeyword, &
                            cellAffe, lAllCellAffe, nbCellAffe)

! ----- Get main parameters for this behaviour
        relaComp = behaviourPrepPara%v_para(iFactorKeyword)%rela_comp
        defoComp = behaviourPrepPara%v_para(iFactorKeyword)%defo_comp
        typeComp = behaviourPrepPara%v_para(iFactorKeyword)%type_comp
        reguVisc = behaviourPrepPara%v_para(iFactorKeyword)%regu_visc
        lMfront = behaviourPrepPara%v_paraExte(iFactorKeyword)%l_mfront_offi .or. &
                  behaviourPrepPara%v_paraExte(iFactorKeyword)%l_mfront_proto
        exteDefo = behaviourPrepPara%v_paraExte(iFactorKeyword)%strain_model
        postIncr = behaviourPrepPara%v_para(iFactorKeyword)%post_incr

! ----- Coding comportment (Python)
        call lccree(1, relaComp, relaCompPY)
        call lccree(1, defoComp, defoCompPY)

! ----- Check the consistency of the modelization with the behaviour
        call compMecaChckModel(iFactorKeyword, &
                               model, fullElemField, &
                               lAllCellAffe, cellAffe, nbCellAffe, &
                               relaComp, relaCompPY, chmate, typeComp, &
                               lElasByDefault, lNeedDeborst, lIncoUpo)

! ----- Select plane stress algorithm
        typeCpla = behaviourPrepPara%v_para(iFactorKeyword)%type_cpla
        call compMecaSelectPlaneStressAlgo(lNeedDeborst, typeCpla)
        behaviourPrepPara%v_para(iFactorKeyword)%type_cpla = typeCpla

! ----- Check the consistency of the strain model with the behaviour
        call compMecaChckStrain(iFactorKeyword, &
                                model, fullElemField, &
                                lAllCellAffe, cellAffe, nbCellAffe, &
                                lMfront, exteDefo, &
                                defoComp, defoCompPY, &
                                relaComp, relaCompPY)

! ----- Check REGU_VISC
        if (reguVisc .ne. 'VIDE') then
            call lctest(relaCompPY, 'REGU_VISC', reguVisc, lctestIret)
            if (lctestIret .eq. 0) then
                call utmess('F', 'COMPOR1_33', nk=2, valk=[reguVisc, relaComp])
            end if
        end if

! ----- Check POST_INCR
        if (postIncr .eq. 'REST_ECRO') then
            call lctest(relaCompPY, 'post_incr', 'REST_ECRO', lctestIret)
            if (lctestIret .eq. 0) then
                call utmess('F', 'COMPOR1_90', nk=1, valk=relaComp)
            end if
        end if

! ----- No Deborst allowed with large strains models
        if (lNeedDeborst .and. defoComp .eq. 'GDEF_LOG') then
            call utmess('F', 'COMPOR1_13')
        end if
        if (lNeedDeborst .and. defoComp .eq. 'SIMO_MIEHE') then
            call utmess('F', 'COMPOR1_13')
        end if

! ----- No INCO_UPO modelization with GDEF_LOG
        if (lIncoUpo .and. defoComp .eq. 'GDEF_LOG') then
            call utmess('F', 'COMPOR1_16')
        end if

! ----- No ENDO_HETEROGENE whith distributed parallelism
        if (relaComp .eq. 'ENDO_HETEROGENE') then
            if (lDistParallel) then
                call utmess('F', 'COMPOR5_25')
            end if
        end if

! ----- Temperature and drying with_BETON_BURGER
        if (relaComp .eq. 'BETON_BURGER') then
            call nmvcd2('TEMP', chmate, exis_temp)
            call nmvcd2('SECH', chmate, exis_sech)
            if  (.not.(exis_temp)) then
                call utmess('F', 'COMPOR6_8')
            end if
            if  (.not.(exis_sech)) then
                call utmess('F', 'COMPOR6_9')
            end if
        end if

! ----- Warning if ELASTIC comportment and initial state
        if (lInitialState .and. typeComp .eq. 'COMP_ELAS') then
            call utmess('A', 'COMPOR1_61')
        end if

! ----- Coding comportment (Python)
        call lcdiscard(relaCompPY)
        call lcdiscard(defoCompPY)

    end do

! - General
    if (lNeedDeborst) then
        call utmess('I', 'COMPOR5_20')
    end if
    if (lElasByDefault) then
        call utmess('I', 'COMPOR5_21')
    end if
    if (lExistVarc .and. behaviourPrepPara%lTotalStrain) then
        call utmess('A', 'COMPOR4_17')
    end if
    if (behaviourPrepPara%nb_comp .eq. 0) then
        call utmess('I', 'COMPOR4_64')
    end if
    if (behaviourPrepPara%nb_comp .ge. 99999) then
        call utmess('A', 'COMPOR4_65')
    end if
!
end subroutine
