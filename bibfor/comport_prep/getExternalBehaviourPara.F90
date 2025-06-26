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
subroutine getExternalBehaviourPara(mesh, v_model_elem, rela_comp, defo_comp, &
                                    kit_comp, prepExte, keywf_, i_comp_, &
                                    elem_type_, type_cpla_in_, type_cpla_out_)
!
    use BehaviourPrepare_type
!
    implicit none
!
#include "asterc/mgis_debug.h"
#include "asterc/mgis_load_library.h"
#include "asterc/umat_get_function.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/BehaviourMGIS_type.h"
#include "asterfort/comp_meca_l.h"
#include "asterfort/comp_read_exte.h"
#include "asterfort/comp_read_mfront.h"
#include "asterfort/comp_read_typmod.h"
#include "asterfort/getExternalStrainModel.h"
!
    character(len=8), intent(in) :: mesh
    integer(kind=8), pointer :: v_model_elem(:)
    character(len=16), intent(in) :: rela_comp
    character(len=16), intent(in) :: defo_comp
    character(len=16), intent(in) :: kit_comp(4)
    type(BehaviourPrep_Exte), intent(inout) :: prepExte
    character(len=16), optional, intent(in) :: keywf_
    integer(kind=8), optional, intent(in) :: i_comp_
    integer(kind=8), optional, intent(in) :: elem_type_
    character(len=16), optional, intent(in) :: type_cpla_in_
    character(len=16), optional, intent(out) :: type_cpla_out_
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of behaviour (mechanics)
!
! Get parameters for external programs (MFRONT/UMAT)
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : name of mesh
! In  v_model_elem     : pointer to list of elements in model
! In  elem_type        : type of element
!                         0 -  Get from affectation
! In  rela_comp        : RELATION comportment
! In  kit_comp         : KIT comportment
! IO  prepExte         : external behaviours parameters
! In  keywf            : factor keyword to read (COMPORTEMENT)
! In  i_comp           : factor keyword index
! In  type_cpla_in     : stress plane hypothesis if known
! Out type_cpla_out    : stress plane hypothesis (for Deborst)
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: l_kit_thm
    character(len=16) :: relaMeca, keywf
    integer(kind=8) :: i_comp
    aster_logical :: l_mfront_proto, l_mfront_offi, l_umat
    character(len=255) :: libr_name, subr_name
    integer(kind=8) :: model_mfront, nbVariUMAT
    character(len=16) :: type_cpla_out, type_cpla_in, extern_addr
    integer(kind=8) :: extern_ptr, extern_type, elem_type, strain_model
!
! --------------------------------------------------------------------------------------------------
!
    l_umat = ASTER_FALSE
    l_mfront_proto = ASTER_FALSE
    l_mfront_offi = ASTER_FALSE
    keywf = 'None'
    i_comp = 0
    libr_name = ' '
    subr_name = ' '
    nbVariUMAT = 0
    model_mfront = MGIS_MODEL_UNSET
    type_cpla_out = 'VIDE'
    if (present(type_cpla_in_)) then
        type_cpla_in = type_cpla_in_
    end if
    elem_type = 0
    if (present(elem_type_)) then
        elem_type = elem_type_
    end if
    strain_model = 0
!
! - Read from command file or not ?
!
    if (present(keywf_)) then
        keywf = keywf_
        i_comp = i_comp_
    end if

    ! - Get mechanical part of behaviour (required for KIT_THM)
    call comp_meca_l(rela_comp, 'KIT_THM', l_kit_thm)
    if (l_kit_thm) then
        relaMeca = kit_comp(1)
    else
        relaMeca = rela_comp
    end if
    !
    ! - Detect type
    !
    call comp_meca_l(relaMeca, 'UMAT', l_umat)
    call comp_meca_l(relaMeca, 'MFRONT_OFFI', l_mfront_offi)
    call comp_meca_l(relaMeca, 'MFRONT_PROTO', l_mfront_proto)
    extern_type = 0
    if (l_mfront_offi) then
        extern_type = 1
    elseif (l_mfront_proto) then
        extern_type = 2
    elseif (l_umat) then
        extern_type = 4
    end if
!
! - Get parameters for external programs (MFRONT/UMAT)
!
    extern_addr = " "
    extern_ptr = 0
    if (l_umat) then
        if (i_comp .ne. 0) then
            call comp_read_exte(keywf, i_comp, libr_name, subr_name, nbVariUMAT)
            call umat_get_function(libr_name, subr_name, extern_ptr)
        end if
!
! - Get pointer and model for MFRONT
!
    elseif (l_mfront_proto .or. l_mfront_offi) then
        ASSERT(i_comp .ne. 0 .or. prepExte%extern_addr .ne. ' ')
        if (i_comp .ne. 0) then
            call comp_read_mfront(keywf, i_comp, extern_addr)
        else
            extern_addr = prepExte%extern_addr
        end if

        if (extern_addr .ne. " ") then

            if (associated(v_model_elem)) then
!               For *_NON_LINE cases
                call comp_read_typmod(mesh, v_model_elem, elem_type, &
                                      keywf, i_comp, rela_comp, type_cpla_in, &
                                      model_mfront, type_cpla_out)
            else
!               For CALC_POINT_MAT case
                model_mfront = MGIS_MODEL_TRIDIMENSIONAL
            end if
!           Get model of strains and load library
            call getExternalStrainModel(defo_comp, strain_model)

            call mgis_load_library(extern_addr, model_mfront, strain_model)
            ! call mgis_debug(extern_addr, "Loaded behaviour:")

        end if
    end if
!
! - Save
!
    if (present(type_cpla_out_)) then
        type_cpla_out_ = type_cpla_out
    end if
    prepExte%extern_type = extern_type
    prepExte%extern_addr = extern_addr
    prepExte%libr_name = libr_name
    prepExte%subr_name = subr_name
    prepExte%extern_ptr = extern_ptr
    prepExte%model_mfront = model_mfront
    prepExte%nbVariUMAT = nbVariUMAT
    prepExte%l_umat = l_umat
    prepExte%l_mfront_proto = l_mfront_proto
    prepExte%l_mfront_offi = l_mfront_offi
    prepExte%strain_model = strain_model
!
end subroutine
