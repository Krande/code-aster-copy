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

subroutine load_neut_spec(type_ther, type_calc, model, time_curr, time, &
                          load_name, load_nume, i_type_neum, nb_type_neumz, nb_in_maxi, &
                          nb_in_prep, lchin, lpain, nb_in_add, lpaout, &
                          load_ligrel, load_option, &
                          time_move_)
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/load_neut_data.h"
#include "asterfort/load_neut_evol.h"
#include "asterfort/load_neut_iden.h"
#include "asterfort/exixfe.h"
#include "asterfort/xajcin.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    character(len=4), intent(in) :: type_ther
    character(len=4), intent(in) :: type_calc
    character(len=24), intent(in) :: model
    real(kind=8), intent(in) :: time_curr
    character(len=24), intent(in) :: time
    character(len=8), intent(in) :: load_name
    integer, intent(in) :: load_nume
    integer, intent(in) :: i_type_neum
    integer, intent(in) :: nb_type_neumz
    integer, intent(in) :: nb_in_maxi
    integer, intent(in) :: nb_in_prep
    character(len=*), intent(inout) :: lpain(nb_in_maxi)
    character(len=*), intent(inout) :: lchin(nb_in_maxi)
    integer, intent(out) :: nb_in_add
    character(len=8), intent(out) :: lpaout
    character(len=19), intent(out) :: load_ligrel
    character(len=16), intent(out) :: load_option
    character(len=24), optional, intent(in) :: time_move_
!
! --------------------------------------------------------------------------------------------------
!
! Neumann loads computation - Thermic
!
! Get information about load (Neumann)
!
! --------------------------------------------------------------------------------------------------
!
! In  type_ther        : type of thermics
!                        'MOVE' for moving sources
!                        'STAT' if not
! In  type_calc        : type of option to compute
!                        '2MBR' for second member (vector)
!                        'RESI' for residual (vector)
!                        'MRIG' for rigidity (matrix)
!                        'MTAN' for tangent matrix
! In  model            : name of the model
! In  time_curr        : current time
! In  time             : time (<CARTE>)
! In  time_move        : modified time (<CARTE>) for THER_NON_LINE_MO
! In  load_name        : name of current load
! In  load_nume        : identification of load type
! In  i_type_neum      : index for Neumann load type
! In  nb_type_neumz    : maximum number of Neumann load type
! In  nb_in_maxi       : maximum number of input fields
! In  nb_in_prep       : number of input fields before specific ones
! IO  lpain            : list of input parameters
! IO  lchin            : list of input fields
! Out nb_in_add        : number of input fields which been added
! Out lpaout           : name of output parameter
! Out load_ligrel      : name of LIGREL for current load
! Out load_option      : name of option for current load
!
! --------------------------------------------------------------------------------------------------
!
    integer :: nb_type_neum
    parameter(nb_type_neum=11)
    aster_logical :: list_load_keyw(nb_type_neum)
!
    integer :: i_field_in
    character(len=19) :: ligrel_load, ligrel_model, cart_name
    logical :: l_constant, l_fonct_0, l_fonct_t
    integer :: ier, nb_obje, i_obje
    logical :: l_xfem
    character(len=6) :: load_type_ligr
    character(len=16) :: load_opti_r, load_opti_f
    character(len=8)  :: load_para_r(2), load_para_f(2), load_para(2)
    character(len=24) :: load_keyw
    character(len=10) :: load_obje(2)
    character(len=19) :: load_name_evol(2)
!
! --------------------------------------------------------------------------------------------------
!
    ligrel_model = model(1:8)//'.MODELE'
    ligrel_load = load_name(1:8)//'.CHTH.LIGRE'
    load_ligrel = ' '
    call exixfe(model, ier)
    l_xfem = ier .ne. 0
    load_option = 'No_Load'
    l_constant = .false.
    l_fonct_0 = .false.
    l_fonct_t = .false.
    i_field_in = nb_in_prep
    ASSERT(i_type_neum .le. nb_type_neum)
    ASSERT(nb_type_neumz .eq. nb_type_neum)
!
! - Identify current load
!
    call load_neut_iden(nb_type_neum, load_name, list_load_keyw)
!
! - Get information about load
!
    if (list_load_keyw(i_type_neum)) then
!
! ----- Value type
!
        if (load_nume .eq. 1) then
            l_constant = .true.
        else if (load_nume .eq. 2) then
            l_fonct_0 = .true.
        else if (load_nume .eq. 3) then
            l_fonct_t = .true.
        else
            ASSERT(.false.)
        end if
!
! ----- Get information about load (Neumann)
!
        call load_neut_data(i_type_neum, nb_type_neumz, type_calc, &
                            load_type_ligr, load_opti_r, load_opti_f, load_para_r, &
                            load_para_f, load_keyw, load_obje, nb_obje)
!
! ----- Input fields - All keyword except EVOL_CHAR
!
        if (load_keyw .ne. 'EVOL_CHAR') then
!
! --------- Set name of option and input parameters
!
            if (l_constant) then
                load_option = load_opti_r
                load_para(1) = load_para_r(1)
                load_para(2) = load_para_r(2)
            else if (l_fonct_0 .or. l_fonct_t) then
                load_option = load_opti_f
                load_para(1) = load_para_f(1)
                load_para(2) = load_para_f(2)
            else
                ASSERT(.false.)
            end if
!
! --------- Set name of input fields
!
            do i_obje = 1, nb_obje
                cart_name = load_name(1:8)//'.CHTH'//load_obje(i_obje)
                i_field_in = i_field_in+1
                lchin(i_field_in) = cart_name(1:19)
                lpain(i_field_in) = load_para(i_obje)
            end do
        end if
!
! ----- Input fields for EVOL_CHAR
!
        if (load_keyw .eq. 'EVOL_CHAR') then
            call load_neut_evol(nb_type_neumz, type_calc, time_curr, load_name, load_type_ligr, &
                                load_opti_r, load_para_r, load_name_evol, nb_obje)
            ASSERT(l_constant)
!
! --------- Set name of option and input parameters
!
            load_option = load_opti_r
            load_para(1) = load_para_r(1)
            if (nb_obje .eq. 2) load_para(2) = load_para_r(2)
!
! --------- Set name of input fields
!
            i_field_in = i_field_in+1
            lchin(i_field_in) = load_name_evol(1)
            lpain(i_field_in) = load_para_r(1)
            if (nb_obje .eq. 2) then
                i_field_in = i_field_in+1
                lchin(i_field_in) = load_name_evol(2)
                lpain(i_field_in) = load_para_r(2)
            end if
        end if
!
! ----- Select time for ECHANGE_PAROI load
!
        i_field_in = i_field_in+1
        lpain(i_field_in) = 'PTEMPSR'
        lchin(i_field_in) = time
        if (load_keyw .eq. 'ECHANGE_PAROI') then
            if (type_ther .eq. 'MOVE') then
                lchin(i_field_in) = time_move_
            end if
        end if
!
! ----- XFEM fields
!
        if (type_calc .eq. '2MBR' .or. type_calc .eq. 'MRIG') then
            if (l_xfem) then
                call xajcin(model, load_option, nb_in_maxi, lchin, lpain, &
                            i_field_in)
            end if
        end if
!
! ----- Select LIGREL
!
        if (l_xfem) then
            load_ligrel = ligrel_model
        else
            if (load_type_ligr .eq. 'Load') then
                load_ligrel = ligrel_load
            elseif (load_type_ligr .eq. 'Model') then
                load_ligrel = ligrel_model
            else
                ASSERT(.false.)
            end if
        end if
!
! ----- Ouput parameter
!
        if (type_calc .eq. '2MBR') then
            lpaout = 'PVECTTR'
        else if (type_calc .eq. 'RESI') then
            lpaout = 'PRESIDU'
        else if (type_calc .eq. 'MRIG') then
            lpaout = 'PMATTTR'
        else if (type_calc .eq. 'MTAN') then
            lpaout = 'PMATTTR'
        else
            ASSERT(.false.)
        end if
!
! ----- Number of input fields which been added
!
        nb_in_add = i_field_in-nb_in_prep
        ASSERT(i_field_in .le. nb_in_maxi)

    end if
!
end subroutine
