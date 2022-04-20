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
! person_in_charge: mickael.abbas at edf.fr
!
subroutine comp_meca_save(model, mesh, chmate, compor, behaviourPrepPara)
!
use Behaviour_type
!
implicit none
!
#include "asterf_types.h"
#include "asterc/getexm.h"
#include "asterfort/assert.h"
#include "asterfort/comp_meca_l.h"
#include "asterfort/comp_read_mesh.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/nmdpmf.h"
#include "asterfort/nocart.h"
#include "asterfort/utmess.h"
#include "asterfort/setBehaviourTypeValue.h"
#include "asterfort/Behaviour_type.h"
!
character(len=8), intent(in) :: model, mesh, chmate
character(len=19), intent(in) :: compor
type(Behaviour_PrepPara), intent(in) :: behaviourPrepPara
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
! In  behaviourPrepPara: datastructure to prepare comportement
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: factorKeyword = 'COMPORTEMENT'
    character(len=24), parameter :: list_elem_affe = '&&COMPMECASAVE.LIST'
    aster_logical :: l_affe_all
    integer, parameter :: nbCmp = COMPOR_SIZE
    integer :: nb_elem_affe, nb_model_affe
    integer, pointer :: v_elem_affe(:) => null()
    integer, pointer :: modelCell(:) => null()
    integer :: i_elem_affe
    integer :: iFactorKeyword, nbFactorKeyword
    character(len=16) :: rela_comp
    character(len=16), pointer :: comporValv(:) => null()
    aster_logical :: l_cristal, l_pmf, l_is_pmf
    integer :: elem_nume
!
! --------------------------------------------------------------------------------------------------
!
    nbFactorKeyword = behaviourPrepPara%nb_comp
    l_is_pmf = ASTER_FALSE

! - Access to MODEL
    call jeveuo(model//'.MAILLE', 'L', vi = modelCell)

! - Access map
    call jeveuo(compor//'.VALV', 'E', vk16 = comporValv)

! - Loop on occurrences of COMPORTEMENT
    do iFactorKeyword = 1, nbFactorKeyword
! ----- Detection of specific cases
        rela_comp = behaviourPrepPara%v_para(iFactorKeyword)%rela_comp
        call comp_meca_l(rela_comp, 'CRISTAL', l_cristal)
        call comp_meca_l(rela_comp, 'PMF', l_pmf)

! ----- Multifiber beams
        if (l_pmf) then
            l_is_pmf = ASTER_TRUE
        endif

! ----- Get elements
        call comp_read_mesh(mesh, factorKeyword, iFactorKeyword,&
                            list_elem_affe, l_affe_all, nb_elem_affe)

! ----- Check if elements belong to model
        nb_model_affe = 0
        if (nb_elem_affe .ne. 0) then
            call jeveuo(list_elem_affe, 'L', vi = v_elem_affe)
            do i_elem_affe = 1, nb_elem_affe
                elem_nume = v_elem_affe(i_elem_affe)
                if (modelCell(elem_nume) .ne. 0) then
                    nb_model_affe = nb_model_affe + 1
                endif
            end do
        endif
        if (.not.l_affe_all) then
            if (nb_model_affe.eq.0) then
                call utmess('A', 'COMPOR4_72', si = iFactorKeyword)
            endif
        endif

! ----- Save informations in the field <COMPOR>
        call setBehaviourTypeValue(behaviourPrepPara%v_para, iFactorKeyword,&
                                   comporMap_ = comporValv)

! ----- Affect in <CARTE>
        if (l_affe_all) then
            call nocart(compor, 1, nbCmp)
        else
            call jeveuo(list_elem_affe, 'L', vi = v_elem_affe)
            call nocart(compor, 3, nbCmp, mode = 'NUM', nma = nb_elem_affe,&
                        limanu = v_elem_affe)
            call jedetr(list_elem_affe)
        endif
    enddo

! - Compor <CARTE> fusing for multifiber beams
    if (l_is_pmf) then
        call nmdpmf(compor, chmate)
    endif
!
    call jedetr(compor//'.NCMP')
    call jedetr(compor//'.VALV')
!
end subroutine
