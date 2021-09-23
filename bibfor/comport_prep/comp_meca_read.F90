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
! person_in_charge: mickael.abbas at edf.fr
!
subroutine comp_meca_read(l_etat_init, ds_compor_prep, model)
!
use Behaviour_type
!
implicit none
!
#include "asterf_types.h"
#include "asterc/getexm.h"
#include "asterc/mfront_get_nbvari.h"
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
type(Behaviour_PrepPara), intent(inout) :: ds_compor_prep
character(len=8), intent(in), optional :: model
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of behaviour (mechanics)
!
! Read informations from command file and catalog
!
! --------------------------------------------------------------------------------------------------
!
! In  l_etat_init      : .true. if initial state is defined
! IO  ds_compor_prep   : datastructure to prepare comportement
! In  model            : name of model
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter:: keywordfact = 'COMPORTEMENT'
    character(len=8) :: mesh
    integer :: i_comp, nb_comp, iret
    character(len=16) :: defo_comp, rela_comp, type_cpla, mult_comp, type_comp, meca_comp
    character(len=16) :: post_iter, defo_ldc, rigi_geom, regu_visc
    character(len=16) :: kit_comp(4), answer
    aster_logical :: l_cristal, l_kit, lNonIncr
    aster_logical :: l_comp_external
    integer, pointer :: v_model_elem(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    nb_comp  = ds_compor_prep%nb_comp
    mesh     = ' '
    lNonIncr = ASTER_FALSE
!
! - Pointer to list of elements in model
!
    if ( present(model) ) then
        call jeveuo(model//'.MAILLE', 'L', vi = v_model_elem)
        call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)
    endif
!
! - Read informations
!
    do i_comp = 1, nb_comp
! ----- Get RELATION from command file
        rela_comp = 'VIDE'
        call compGetRelation(i_comp, rela_comp)

! ----- Detection of specific cases
        call comp_meca_l(rela_comp, 'KIT'    , l_kit)
        call comp_meca_l(rela_comp, 'CRISTAL', l_cristal)

! ----- Get DEFORMATION from command file
        defo_comp = 'VIDE'
        call getvtx(keywordfact, 'DEFORMATION', iocc = i_comp, scal = defo_comp)

! ----- Get RIGI_GEOM from command file
        rigi_geom = ' '
        if (getexm(keywordfact,'RIGI_GEOM') .eq. 1) then
            call getvtx(keywordfact, 'RIGI_GEOM', iocc = i_comp, scal=rigi_geom, nbret=iret)
            if (iret .eq. 0) then
                rigi_geom = 'VIDE'
            end if
        end if

! ----- Damage post-treatment
        post_iter = 'VIDE'
        if (getexm(keywordfact,'POST_ITER') .eq. 1) then
            call getvtx(keywordfact, 'POST_ITER', iocc = i_comp, scal=post_iter, nbret=iret)
            if (iret .eq. 0) then
                post_iter = 'VIDE'
            endif
        endif
! ----- Viscuous regularization
        regu_visc = 'VIDE'
        if (getexm(keywordfact,'REGU_VISC') .eq. 1) then
            call getvtx(keywordfact, 'REGU_VISC', iocc = i_comp, scal=answer)
            if (answer .eq. 'OUI') then
                regu_visc = 'REGU_VISC_ELAS'
            elseif (answer .eq. 'NON') then
                regu_visc = 'VIDE'
            else
                ASSERT(ASTER_FALSE)
            endif
        endif

! ----- For KIT
        kit_comp = 'VIDE'
        if (l_kit) then
            call comp_meca_rkit(keywordfact, i_comp, rela_comp, kit_comp, l_etat_init)
        endif

! ----- Get mechanical part of behaviour
        meca_comp = 'VIDE'
        call compGetMecaPart(rela_comp, kit_comp, meca_comp)

! ----- Get multi-comportment *CRISTAL
        mult_comp = 'VIDE'
        if (l_cristal) then
            call getvid(keywordfact, 'COMPOR', iocc = i_comp, scal = mult_comp)
        endif

! ----- Get parameters for external programs (MFRONT/UMAT)
        type_cpla = 'VIDE'
        call getExternalBehaviourPara(mesh           , v_model_elem, rela_comp, kit_comp,&
                                      l_comp_external, ds_compor_prep%v_paraExte(i_comp),&
                                      keywordfact    , i_comp,&
                                      type_cpla_out_ = type_cpla)

! ----- Select type of behaviour (incremental or total)
        type_comp = 'VIDE'
        call comp_meca_incr(rela_comp, defo_comp, type_comp, l_etat_init)
        if (type_comp .eq. 'COMP_ELAS') then
            lNonIncr = ASTER_TRUE
        endif

! ----- Select type of strain (mechanical or total) from catalog
        defo_ldc = 'VIDE'
        call comp_meca_deflc(rela_comp, defo_comp, defo_ldc)

! ----- Save parameters
        ds_compor_prep%v_para(i_comp)%rela_comp = rela_comp
        ds_compor_prep%v_para(i_comp)%meca_comp = meca_comp
        ds_compor_prep%v_para(i_comp)%defo_comp = defo_comp
        ds_compor_prep%v_para(i_comp)%type_comp = type_comp
        ds_compor_prep%v_para(i_comp)%type_cpla = type_cpla
        ds_compor_prep%v_para(i_comp)%kit_comp  = kit_comp
        ds_compor_prep%v_para(i_comp)%mult_comp = mult_comp
        ds_compor_prep%v_para(i_comp)%post_iter = post_iter
        ds_compor_prep%v_para(i_comp)%defo_ldc  = defo_ldc
        ds_compor_prep%v_para(i_comp)%rigi_geom = rigi_geom
        ds_compor_prep%v_para(i_comp)%regu_visc = regu_visc
    end do

! - Is at least ONE behaviour is not incremental ?
    ds_compor_prep%lNonIncr = lNonIncr
!
end subroutine
