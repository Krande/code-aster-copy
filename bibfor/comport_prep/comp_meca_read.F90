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
subroutine comp_meca_read(l_etat_init, behaviourPrepPara, model)
!
    use Behaviour_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/getvtx.h"
#include "asterfort/assert.h"
#include "asterfort/jeveuo.h"
#include "asterfort/compGetMecaPart.h"
#include "asterfort/compGetRelation.h"
#include "asterfort/comp_meca_incr.h"
#include "asterfort/comp_meca_deflc.h"
#include "asterfort/getExternalBehaviourPara.h"
#include "asterfort/comp_meca_rkit.h"
#include "asterfort/comp_meca_l.h"
!
    aster_logical, intent(in) :: l_etat_init
    type(Behaviour_PrepPara), intent(inout) :: behaviourPrepPara
    character(len=8), intent(in), optional :: model
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of constitutive laws (mechanics)
!
! Read from command file
!
! --------------------------------------------------------------------------------------------------
!
! In  l_etat_init      : .true. if initial state is defined
! IO  behaviourPrepPara: datastructure to prepare behaviour
! In  model            : model
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter:: factorKeyword = 'COMPORTEMENT'
    character(len=8) :: mesh
    integer :: iFactorKeyword, nbFactorKeyword, iret
    character(len=16) :: defo_comp, rela_comp, type_cpla, mult_comp, type_comp, meca_comp
    character(len=16) :: post_iter, defo_ldc, rigi_geom, regu_visc
    character(len=16) :: kit_comp(4), answer
    aster_logical :: l_cristal, l_kit, lTotalStrain
    integer, pointer :: modelCell(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    nbFactorKeyword = behaviourPrepPara%nb_comp
    mesh = ' '
    lTotalStrain = ASTER_FALSE

! - Pointer to list of elements in model
    if (present(model)) then
        call jeveuo(model//'.MAILLE', 'L', vi=modelCell)
        call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)
    end if

! - Read informations
    do iFactorKeyword = 1, nbFactorKeyword
! ----- Get RELATION from command file
        rela_comp = 'VIDE'
        call compGetRelation(factorKeyword, iFactorKeyword, rela_comp)

! ----- Detection of specific cases
        call comp_meca_l(rela_comp, 'KIT', l_kit)
        call comp_meca_l(rela_comp, 'CRISTAL', l_cristal)

! ----- Get DEFORMATION from command file
        defo_comp = 'VIDE'
        call getvtx(factorKeyword, 'DEFORMATION', iocc=iFactorKeyword, scal=defo_comp)

! ----- Get RIGI_GEOM from command file
        rigi_geom = ' '
        call getvtx(factorKeyword, 'RIGI_GEOM', iocc=iFactorKeyword, &
                    scal=rigi_geom, nbret=iret)
        if (iret .eq. 0) then
            rigi_geom = 'VIDE'
        end if

! ----- Damage post-treatment
        post_iter = 'VIDE'
        call getvtx(factorKeyword, 'POST_ITER', iocc=iFactorKeyword, &
                    scal=post_iter, nbret=iret)
        if (iret .eq. 0) then
            post_iter = 'VIDE'
        end if

! ----- Viscuous regularization
        regu_visc = 'VIDE'
        call getvtx(factorKeyword, 'REGU_VISC', iocc=iFactorKeyword, scal=answer, nbret=iret)
        if (iret .eq. 1) then
            if (answer .eq. 'OUI') then
                regu_visc = 'REGU_VISC_ELAS'
            elseif (answer .eq. 'NON') then
                regu_visc = 'VIDE'
            else
                ASSERT(ASTER_FALSE)
            end if
        end if

! ----- For KIT
        kit_comp = 'VIDE'
        if (l_kit) then
            call comp_meca_rkit(factorKeyword, iFactorKeyword, rela_comp, kit_comp, l_etat_init)
        end if

! ----- Get mechanical part of behaviour
        meca_comp = 'VIDE'
        call compGetMecaPart(rela_comp, kit_comp, meca_comp)

! ----- Get multi-material *CRISTAL
        mult_comp = 'VIDE'
        if (l_cristal) then
            call getvid(factorKeyword, 'COMPOR', iocc=iFactorKeyword, scal=mult_comp)
        end if

! ----- Get parameters for external programs (MFRONT/UMAT)
        type_cpla = 'VIDE'
        call getExternalBehaviourPara(mesh, modelCell, rela_comp, defo_comp, kit_comp, &
                                      behaviourPrepPara%v_paraExte(iFactorKeyword), &
                                      factorKeyword, iFactorKeyword, &
                                      type_cpla_out_=type_cpla)

! ----- Select type of behaviour (incremental or total)
        type_comp = 'VIDE'
        call comp_meca_incr(rela_comp, defo_comp, type_comp, l_etat_init)

! ----- Select type of strain (mechanical or total) from catalog
        defo_ldc = 'VIDE'
        call comp_meca_deflc(rela_comp, defo_comp, defo_ldc)
        if (defo_ldc .eq. 'TOTALE') then
            lTotalStrain = ASTER_TRUE
        end if

! ----- Save parameters
        behaviourPrepPara%v_para(iFactorKeyword)%rela_comp = rela_comp
        behaviourPrepPara%v_para(iFactorKeyword)%meca_comp = meca_comp
        behaviourPrepPara%v_para(iFactorKeyword)%defo_comp = defo_comp
        behaviourPrepPara%v_para(iFactorKeyword)%type_comp = type_comp
        behaviourPrepPara%v_para(iFactorKeyword)%type_cpla = type_cpla
        behaviourPrepPara%v_para(iFactorKeyword)%kit_comp = kit_comp
        behaviourPrepPara%v_para(iFactorKeyword)%mult_comp = mult_comp
        behaviourPrepPara%v_para(iFactorKeyword)%post_iter = post_iter
        behaviourPrepPara%v_para(iFactorKeyword)%defo_ldc = defo_ldc
        behaviourPrepPara%v_para(iFactorKeyword)%rigi_geom = rigi_geom
        behaviourPrepPara%v_para(iFactorKeyword)%regu_visc = regu_visc

    end do

    if (behaviourPrepPara%lDebug) then
        WRITE (6, *) "Données lues: ", nbFactorKeyword, " occurrences."
        do iFactorKeyword = 1, nbFactorKeyword
            WRITE (6, *) "- Occurrence : ", iFactorKeyword
            WRITE (6, *) "--- rela_comp : ", behaviourPrepPara%v_para(iFactorKeyword)%rela_comp
            WRITE (6, *) "--- meca_comp : ", behaviourPrepPara%v_para(iFactorKeyword)%meca_comp
            WRITE (6, *) "--- defo_comp : ", behaviourPrepPara%v_para(iFactorKeyword)%defo_comp
            WRITE (6, *) "--- type_comp : ", behaviourPrepPara%v_para(iFactorKeyword)%type_comp
            WRITE (6, *) "--- type_cpla : ", behaviourPrepPara%v_para(iFactorKeyword)%type_cpla
            WRITE (6, *) "--- kit_comp  : ", behaviourPrepPara%v_para(iFactorKeyword)%kit_comp
            WRITE (6, *) "--- mult_comp : ", behaviourPrepPara%v_para(iFactorKeyword)%mult_comp
            WRITE (6, *) "--- post_iter : ", behaviourPrepPara%v_para(iFactorKeyword)%post_iter
            WRITE (6, *) "--- defo_ldc  : ", behaviourPrepPara%v_para(iFactorKeyword)%defo_ldc
            WRITE (6, *) "--- rigi_geom : ", behaviourPrepPara%v_para(iFactorKeyword)%rigi_geom
            WRITE (6, *) "--- regu_visc : ", behaviourPrepPara%v_para(iFactorKeyword)%regu_visc
        end do
    end if

! - Is at least ONE behaviour is not incremental ?
    behaviourPrepPara%lTotalStrain = lTotalStrain
!
end subroutine
