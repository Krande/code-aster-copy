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
subroutine comp_meca_save(model, mesh, chmate, compor, prepMapCompor)
!
    use BehaviourPrepare_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/comp_meca_l.h"
#include "asterfort/comp_read_mesh.h"
#include "asterfort/dismoi.h"
#include "asterfort/isParallelMesh.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nmdpmf.h"
#include "asterfort/nocart.h"
#include "asterfort/setBehaviourTypeValue.h"
#include "asterfort/utmess.h"
!
    character(len=8), intent(in) :: model, mesh, chmate
    character(len=19), intent(in) :: compor
    type(BehaviourPrep_MapCompor), intent(in) :: prepMapCompor
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of behaviour (mechanics)
!
! Save informations in COMPOR <CARTE>
!
! --------------------------------------------------------------------------------------------------
!
! In  model            : model
! In  mesh             : mesh
! In  chmate           : material field
! In  compor           : map for parameters of constitutive laws
! In  prepMapCompor    : datastructure to construct COMPOR map
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: factorKeyword = 'COMPORTEMENT'
    character(len=24), parameter :: list_elem_affe = '&&COMPMECASAVE.LIST'
    aster_logical :: l_affe_all
    integer(kind=8), parameter :: nbCmp = COMPOR_SIZE
    integer(kind=8) :: nb_elem_affe, nb_model_affe
    integer(kind=8), pointer :: v_elem_affe(:) => null()
    integer(kind=8), pointer :: modelCell(:) => null()
    integer(kind=8) :: i_elem_affe
    integer(kind=8) :: iFactorKeyword, nbFactorKeyword
    character(len=16) :: rela_comp
    character(len=16), pointer :: comporValv(:) => null()
    character(len=19) :: ligrel
    aster_logical :: l_cristal, l_pmf, l_is_pmf, l_parallel_mesh
    integer(kind=8) :: elem_nume
!
! --------------------------------------------------------------------------------------------------
!
    nbFactorKeyword = prepMapCompor%nb_comp
    l_is_pmf = ASTER_FALSE
    l_parallel_mesh = isParallelMesh(mesh)

! - Access to MODEL
    call dismoi('NOM_LIGREL', model, 'MODELE', repk=ligrel)
    call jeveuo(ligrel//'.TYFE', 'L', vi=modelCell)

! - Access map
    call jeveuo(compor//'.VALV', 'E', vk16=comporValv)

! - Loop on occurrences of COMPORTEMENT
    do iFactorKeyword = 1, nbFactorKeyword
! ----- Detection of specific cases
        rela_comp = prepMapCompor%prepPara(iFactorKeyword)%rela_comp
        call comp_meca_l(rela_comp, 'CRISTAL', l_cristal)
        call comp_meca_l(rela_comp, 'PMF', l_pmf)

! ----- Multifiber beams
        if (l_pmf) then
            l_is_pmf = ASTER_TRUE
        end if

! ----- Get elements
        call comp_read_mesh(mesh, factorKeyword, iFactorKeyword, &
                            list_elem_affe, l_affe_all, nb_elem_affe)

! ----- Check if elements belong to model
        nb_model_affe = 0
        if (nb_elem_affe .ne. 0) then
            call jeveuo(list_elem_affe, 'L', vi=v_elem_affe)
            do i_elem_affe = 1, nb_elem_affe
                elem_nume = v_elem_affe(i_elem_affe)
                if (modelCell(elem_nume) .ne. 0) then
                    nb_model_affe = nb_model_affe+1
                end if
            end do
        end if
        if (.not. l_affe_all) then
            if (nb_model_affe .eq. 0) then
                call utmess('A', 'COMPOR4_72', si=iFactorKeyword)
            end if
        end if

! ----- Save informations in the field <COMPOR>
        call setBehaviourTypeValue(prepMapCompor, iFactorKeyword, &
                                   comporMap_=comporValv)

! ----- Affect in <CARTE>
        if (l_affe_all) then
            call nocart(compor, 1, nbCmp)
        else
            if (nb_elem_affe > 0 .or. .not. l_parallel_mesh) then
                call jeveuo(list_elem_affe, 'L', vi=v_elem_affe)
                call nocart(compor, 3, nbCmp, mode='NUM', nma=nb_elem_affe, &
                            limanu=v_elem_affe)
            end if
            call jedetr(list_elem_affe)
        end if
    end do

! - Compor <CARTE> fusing for multifiber beams
    if (l_is_pmf) then
        call nmdpmf(compor, chmate)
    end if
!
    call jedetr(compor//'.NCMP')
    call jedetr(compor//'.VALV')
!
end subroutine
