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
subroutine comp_meca_read(l_etat_init, prepMapCompor, model)
!
    use BehaviourPrepare_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/comp_meca_deflc.h"
#include "asterfort/comp_meca_incr.h"
#include "asterfort/comp_meca_l.h"
#include "asterfort/comp_meca_rkit.h"
#include "asterfort/compGetMecaPart.h"
#include "asterfort/compGetRelation.h"
#include "asterfort/dismoi.h"
#include "asterfort/getExternalBehaviourPara.h"
#include "asterfort/getvid.h"
#include "asterfort/getvtx.h"
#include "asterfort/jeveuo.h"
!
    aster_logical, intent(in) :: l_etat_init
    type(BehaviourPrep_MapCompor), intent(inout) :: prepMapCompor
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
! IO  prepMapCompor    : datastructure to construct COMPOR map
! In  model            : model
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter:: factorKeyword = 'COMPORTEMENT'
    character(len=8) :: mesh
    integer(kind=8) :: iFactorKeyword, nbFactorKeyword, iret
    character(len=16) :: defo_comp, rela_comp, type_cpla, mult_comp, type_comp, meca_comp
    character(len=16) :: post_iter, defo_ldc, rigi_geom, regu_visc, post_incr
    character(len=16) :: kit_comp(4), answer
    character(len=19) :: modelLigrel
    aster_logical :: l_cristal, l_kit, lTotalStrain
    integer(kind=8), pointer :: modelCell(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    nbFactorKeyword = prepMapCompor%nb_comp
    mesh = ' '
    lTotalStrain = ASTER_FALSE

! - Pointer to list of elements in model
    if (present(model)) then
        call dismoi('NOM_LIGREL', model, 'MODELE', repk=modelLigrel)
        call jeveuo(modelLigrel//'.TYFE', 'L', vi=modelCell)
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

! ----- Post-treatment at each Newton iteration
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

! ----- Post-treatment at each time step
        post_incr = "VIDE"
        call getvtx(factorKeyword, 'POST_INCR', iocc=iFactorKeyword, &
                    scal=post_incr, nbret=iret)
        if (iret .eq. 0) then
            post_incr = 'VIDE'
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
                                      prepMapCompor%prepExte(iFactorKeyword), &
                                      factorKeyword, iFactorKeyword, &
                                      type_cpla_out_=type_cpla)

! ----- Select type of behaviour (incremental or total)
        type_comp = 'VIDE'
        call comp_meca_incr(rela_comp, defo_comp, type_comp, l_etat_init)

! ----- Select type of strain (mechanical or total) from catalog
        defo_ldc = 'VIDE'
        call comp_meca_deflc(rela_comp, defo_comp, defo_ldc)
        lTotalStrain = defo_ldc .eq. 'TOTALE'

! ----- Save parameters
        prepMapCompor%prepPara(iFactorKeyword)%rela_comp = rela_comp
        prepMapCompor%prepPara(iFactorKeyword)%meca_comp = meca_comp
        prepMapCompor%prepPara(iFactorKeyword)%defo_comp = defo_comp
        prepMapCompor%prepPara(iFactorKeyword)%type_comp = type_comp
        prepMapCompor%prepPara(iFactorKeyword)%type_cpla = type_cpla
        prepMapCompor%prepPara(iFactorKeyword)%kit_comp = kit_comp
        prepMapCompor%prepPara(iFactorKeyword)%mult_comp = mult_comp
        prepMapCompor%prepPara(iFactorKeyword)%post_iter = post_iter
        prepMapCompor%prepPara(iFactorKeyword)%defo_ldc = defo_ldc
        prepMapCompor%prepPara(iFactorKeyword)%rigi_geom = rigi_geom
        prepMapCompor%prepPara(iFactorKeyword)%regu_visc = regu_visc
        prepMapCompor%prepPara(iFactorKeyword)%post_incr = post_incr
        prepMapCompor%prepPara(iFactorKeyword)%lTotalStrain = lTotalStrain
    end do

    if (prepMapCompor%lDebug) then
        WRITE (6, *) "Données lues: ", nbFactorKeyword, " occurrences."
        do iFactorKeyword = 1, nbFactorKeyword
            WRITE (6, *) "- Occurrence : ", iFactorKeyword
            WRITE (6, *) "--- rela_comp : ", prepMapCompor%prepPara(iFactorKeyword)%rela_comp
            WRITE (6, *) "--- meca_comp : ", prepMapCompor%prepPara(iFactorKeyword)%meca_comp
            WRITE (6, *) "--- defo_comp : ", prepMapCompor%prepPara(iFactorKeyword)%defo_comp
            WRITE (6, *) "--- type_comp : ", prepMapCompor%prepPara(iFactorKeyword)%type_comp
            WRITE (6, *) "--- type_cpla : ", prepMapCompor%prepPara(iFactorKeyword)%type_cpla
            WRITE (6, *) "--- kit_comp  : ", prepMapCompor%prepPara(iFactorKeyword)%kit_comp
            WRITE (6, *) "--- mult_comp : ", prepMapCompor%prepPara(iFactorKeyword)%mult_comp
            WRITE (6, *) "--- post_iter : ", prepMapCompor%prepPara(iFactorKeyword)%post_iter
            WRITE (6, *) "--- defo_ldc  : ", prepMapCompor%prepPara(iFactorKeyword)%defo_ldc
            WRITE (6, *) "--- rigi_geom : ", prepMapCompor%prepPara(iFactorKeyword)%rigi_geom
            WRITE (6, *) "--- regu_visc : ", prepMapCompor%prepPara(iFactorKeyword)%regu_visc
            WRITE (6, *) "--- post_incr : ", prepMapCompor%prepPara(iFactorKeyword)%post_incr
            WRITE (6, *) "--- total strain : ", prepMapCompor%prepPara(iFactorKeyword)%lTotalStrain
        end do
    end if
!
end subroutine
