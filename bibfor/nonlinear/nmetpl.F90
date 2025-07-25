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

subroutine nmetpl(ds_inout, sd_suiv, sd_obsv)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/jeveuo.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    type(NL_DS_InOut), intent(inout) :: ds_inout
    character(len=24), intent(in) :: sd_suiv
    character(len=19), intent(in) :: sd_obsv
!
! --------------------------------------------------------------------------------------------------
!
! *_NON_LINE - Input/output datastructure
!
! Update name of field in algorithm
!
! This utiliy is required for "hat" variables
!
! --------------------------------------------------------------------------------------------------
!
! IO  ds_inout         : datastructure for input/output management
! In  sd_suiv          : datastructure for dof monitor parameters
! In  sd_obsv          : datastructure for observation parameters
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nb_field, i_field
    integer(kind=8) :: i_keyw_fact, nb_keyw_fact
    character(len=24) :: algo_name_old, algo_name_new
    character(len=6) :: hat_type, hat_vari
    character(len=24) :: field_type
    character(len=14) :: sdextr_suiv, sdextr_obsv
    character(len=24) :: extr_info
    integer(kind=8), pointer :: v_extr_info(:) => null()
    character(len=24) :: extr_field
    character(len=24), pointer :: v_extr_field(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    nb_field = ds_inout%nb_field
!
! - For storing
!
    do i_field = 1, nb_field
        algo_name_old = ds_inout%field(i_field)%algo_name
        if (algo_name_old(1:3) .eq. '#H#') then
            hat_type = algo_name_old(4:9)
            hat_vari = algo_name_old(11:16)
            if (hat_type .eq. 'VALINC') then
                if (hat_vari .eq. 'TEMP') then
                    ASSERT(.false.)
                end if
                algo_name_new = algo_name_old
                algo_name_new(14:16) = 'PLU'
                ds_inout%field(i_field)%algo_name = algo_name_new
            end if
        end if
    end do
!
! - For DOF monitoring
!
    sdextr_suiv = sd_suiv(1:14)
    extr_info = sdextr_suiv(1:14)//'     .INFO'
    call jeveuo(extr_info, 'L', vi=v_extr_info)
    nb_keyw_fact = v_extr_info(1)
    nb_field = v_extr_info(6)
    extr_field = sdextr_suiv(1:14)//'     .CHAM'
    if (nb_keyw_fact .ne. 0) then
        call jeveuo(extr_field, 'E', vk24=v_extr_field)
    end if
    do i_keyw_fact = 1, nb_keyw_fact
        i_field = v_extr_info(7+7*(i_keyw_fact-1)+7)
        field_type = v_extr_field(4*(i_field-1)+1)
        if (field_type .ne. 'NONE') then
            algo_name_old = v_extr_field(4*(i_field-1)+4) (1:19)
            if (algo_name_old(13:15) .eq. 'MOI') then
                algo_name_new = algo_name_old
                algo_name_new(13:15) = 'PLU'
                v_extr_field(4*(i_field-1)+4) = algo_name_new
            end if
        end if
    end do
!
! - For observation
!
    sdextr_obsv = sd_obsv(1:14)
    extr_info = sdextr_obsv(1:14)//'     .INFO'
    call jeveuo(extr_info, 'L', vi=v_extr_info)
    nb_keyw_fact = v_extr_info(1)
    nb_field = v_extr_info(6)
    extr_field = sdextr_obsv(1:14)//'     .CHAM'
    if (nb_keyw_fact .ne. 0) then
        call jeveuo(extr_field, 'E', vk24=v_extr_field)
    end if
    do i_keyw_fact = 1, nb_keyw_fact
        i_field = v_extr_info(7+7*(i_keyw_fact-1)+7)
        field_type = v_extr_field(4*(i_field-1)+1)
        if (field_type .ne. 'NONE') then
            algo_name_old = v_extr_field(4*(i_field-1)+4) (1:19)
            if (algo_name_old(13:15) .eq. 'MOI') then
                algo_name_new = algo_name_old
                algo_name_new(13:15) = 'PLU'
                v_extr_field(4*(i_field-1)+4) = algo_name_new
            end if
        end if
    end do
!
end subroutine
