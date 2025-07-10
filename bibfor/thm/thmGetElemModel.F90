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
subroutine thmGetElemModel(ds_thm, l_axi_, l_vf_, ndim_, type_elem_)
!
    use THM_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/lteatt.h"
!
    type(THM_DS), intent(inout) :: ds_thm
    aster_logical, optional, intent(out) :: l_axi_, l_vf_
    integer(kind=8), optional, intent(out) :: ndim_
    character(len=8), optional, intent(out) :: type_elem_(2)
!
! --------------------------------------------------------------------------------------------------
!
! THM - Parameters
!
! Get model of finite element
!
! --------------------------------------------------------------------------------------------------
!
! IO  ds_thm           : datastructure for THM
! Out l_axi            : flag is axisymmetric model
! Out l_vf             : flag for finite volume
! Out l_steady         : .true. for steady state
! Out ndim             : dimension of element (2 ou 3)
! Out type_elem        : type of element
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: l_axi, l_vf
    integer(kind=8) :: ndim
    character(len=8) :: type_elem(2)
!
! --------------------------------------------------------------------------------------------------
!
    l_axi = ASTER_FALSE
    ndim = 0
    type_elem = ' '
    l_vf = ASTER_FALSE

! - Get dof in finite element
    ds_thm%ds_elem%l_dof_ther = lteatt('THER', 'OUI')
    ds_thm%ds_elem%l_dof_meca = lteatt('MECA', 'OUI')
    ds_thm%ds_elem%l_dof_pre1 = .not. lteatt('HYDR1', '0')
    ds_thm%ds_elem%l_dof_pre2 = .not. lteatt('HYDR2', '0')
    ds_thm%ds_elem%l_dof_2nd = lteatt('DIL', 'OUI')
    ds_thm%ds_elem%nb_phase(1) = 0
    ds_thm%ds_elem%nb_phase(2) = 0
    if (lteatt('HYDR1', '1')) then
        ds_thm%ds_elem%nb_phase(1) = 1
    end if
    if (lteatt('HYDR1', '2')) then
        ds_thm%ds_elem%nb_phase(1) = 2
    end if
    if (lteatt('HYDR2', '1')) then
        ds_thm%ds_elem%nb_phase(2) = 1
    end if
    if (lteatt('HYDR2', '2')) then
        ds_thm%ds_elem%nb_phase(2) = 2
    end if

! - Get general model
    l_axi = lteatt('AXIS', 'OUI')
    if (l_axi) then
        type_elem(1) = 'AXIS'
        ndim = 2
    else if (lteatt('D_PLAN', 'OUI')) then
        type_elem(1) = 'D_PLAN'
        ndim = 2
    else
        type_elem(1) = '3D'
        ndim = 3
    end if
    type_elem(2) = 'THM'

! - Finite volume ?
    l_vf = lteatt('TYPMOD3', 'SUSHI')

! - Copy
    if (present(l_axi_)) then
        l_axi_ = l_axi
    end if
    if (present(l_vf_)) then
        l_vf_ = l_vf
    end if
    if (present(ndim_)) then
        ndim_ = ndim
    end if
    if (present(type_elem_)) then
        type_elem_ = type_elem
    end if
!
end subroutine
