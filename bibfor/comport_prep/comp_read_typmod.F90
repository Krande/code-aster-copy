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
subroutine comp_read_typmod(mesh, v_model_elem, elem_type, &
                            keywf, i_comp, rela_comp, type_cpla_in, &
                            model_mfront, type_cpla_out)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/BehaviourMGIS_type.h"
#include "asterfort/comp_mfront_modelem.h"
#include "asterfort/comp_read_mesh.h"
#include "asterfort/dismoi.h"
#include "asterfort/getMFrontPlaneStress.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
!
    character(len=8), intent(in) :: mesh
    integer(kind=8), pointer :: v_model_elem(:)
    integer(kind=8), intent(in) :: elem_type
    character(len=16), intent(in) :: keywf
    integer(kind=8), intent(in) :: i_comp
    character(len=16), intent(in) :: rela_comp
    character(len=16), intent(in) :: type_cpla_in
    integer(kind=8), intent(out) :: model_mfront
    character(len=16), intent(out) :: type_cpla_out
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of comportment (mechanics)
!
! Find dimension and type of modelisation for MFront
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : name of mesh
! In  v_model_elem     : pointer to list of elements in model
! In  elem_type        : type of element
!                         0 -  Get from affectation
! In  keywf            : factor keyword to read (COMPORTEMENT)
! In  i_comp           : factor keyword index
! In  rela_comp        : RELATION comportment
! In  type_cpla_in     : stress plane hypothesis if known
! Out model_mfront     : type of modelisation MFront
! Out type_cpla_out    : stress plane hypothesis (for Deborst)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nb_elem_affe, nb_elem, i_elem, elem_nume, model_save
    integer(kind=8) :: elem_type_nume, codret
    aster_logical :: l_affe_all, l_mfront_cp
    character(len=24) :: list_elem_affe
    character(len=16) :: elem_type_name
    integer(kind=8), pointer :: v_elem_affe(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    model_mfront = MGIS_MODEL_UNSET
    model_save = MGIS_MODEL_UNSET
    list_elem_affe = '&&COMPMECASAVE.LIST'
    type_cpla_out = 'VIDE'
!
! - Get current element number
!
    if (elem_type .eq. 0) then
        call comp_read_mesh(mesh, keywf, i_comp, &
                            list_elem_affe, l_affe_all, nb_elem_affe)
        if (l_affe_all) then
            call dismoi('NB_MA_MAILLA', mesh, 'MAILLAGE', repi=nb_elem)
        else
            call jeveuo(list_elem_affe, 'L', vi=v_elem_affe)
            nb_elem = nb_elem_affe
        end if
    else
        nb_elem = 1
    end if
!
! - For plane stress hypothesis
!
    if (i_comp .eq. 0) then
        l_mfront_cp = type_cpla_in .eq. 'ANALYTIQUE'
    else
        call getMFrontPlaneStress(keywf, i_comp, rela_comp, l_mfront_cp)
    end if
!
! - Loop on elements
!
    do i_elem = 1, nb_elem
! ----- Current element
        if (elem_type .eq. 0) then
            if (l_affe_all) then
                elem_nume = i_elem
            else
                elem_nume = v_elem_affe(i_elem)
            end if
            elem_type_nume = v_model_elem(elem_nume)
        else
            elem_type_nume = elem_type
        end if
! ----- Select type of modelisation for MFront
        if (elem_type_nume .ne. 0) then
            call jenuno(jexnum('&CATA.TE.NOMTE', elem_type_nume), elem_type_name)
            call comp_mfront_modelem(elem_type_name, l_mfront_cp, &
                                     model_mfront, &
                                     codret, type_cpla_out)
            if (model_mfront .ne. MGIS_MODEL_UNSET) then
                if (model_save .eq. MGIS_MODEL_UNSET) then
                    model_save = model_mfront
                else
                    if ((model_save .ne. model_mfront)) then
                        codret = 1
                    end if
                end if
            end if
            if (codret .eq. 1) then
                call utmess('F', 'COMPOR4_13', ni=2, &
                            vali=[model_save, model_mfront], &
                            sk="MGISBehaviourFort.h")
            end if
            if (codret .eq. 2) then
                call utmess('F', 'COMPOR4_14', si=model_mfront, &
                            sk="MGISBehaviourFort.h")
            end if
        end if
    end do
!
! - Final model for MFront
!
    model_mfront = model_save
!
end subroutine
