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

subroutine load_neut_excl(command, &
                          lload_name_, lload_info_, &
                          list_load_, &
                          list_nbload_, list_name_)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/utmess.h"
#include "asterfort/load_list_info.h"
#include "asterfort/load_neut_iden.h"
#include "asterfort/load_neut_data.h"
!
!
    character(len=*), intent(in) :: command
    character(len=19), optional, intent(in) :: list_load_
    character(len=*), optional, intent(in) :: lload_name_
    character(len=*), optional, intent(in) :: lload_info_
    character(len=*), optional, target, intent(in) :: list_name_(*)
    integer, optional, intent(in) :: list_nbload_
!
! --------------------------------------------------------------------------------------------------
!
! Neumann loads - Thermic
!
! Exclusion of loads in commands
!
! --------------------------------------------------------------------------------------------------
!
! In  command          : current command
! In  list_load        : list of loads
! In  lload_name       : name of JEVEUX object in list_load datastructure for loads name
! In  lload_info       : name of JEVEUX object in list_load datastructure for loads info
! In  list_nbload      : length of list_name
! In  list_name        : list of loads name
!
! --------------------------------------------------------------------------------------------------
!
    integer :: nb_type_neum
    parameter(nb_type_neum=11)
    aster_logical :: list_load_keyw(nb_type_neum)
    integer :: nb_excl_tnlm, nb_excl_cve, nb_excl_cme
    parameter(nb_excl_tnlm=1, nb_excl_cve=1, nb_excl_cme=1)
    integer :: nb_excl_ce, nb_excl_cma, nb_excl_cfa
    parameter(nb_excl_ce=1, nb_excl_cma=1, nb_excl_cfa=1)
    character(len=24) :: list_excl_tnlm(nb_excl_tnlm)
    character(len=24) :: list_excl_cve(nb_excl_cve)
    character(len=24) :: list_excl_cme(nb_excl_cme)
    character(len=24) :: list_excl_ce(nb_excl_ce)
    character(len=24) :: list_excl_cfa(nb_excl_cfa)
    character(len=24) :: list_excl_cma(nb_excl_cma)
!
    integer :: i_load, i_type_neum, nb_load, i_excl
    character(len=8) :: load_name
    character(len=24) :: load_keyw
    aster_logical :: load_empty
    character(len=24), pointer :: v_load_name(:) => null()
    integer, pointer :: v_load_info(:) => null()
! - For THER_NON_LINE_MO
    data list_excl_tnlm/'EVOL_CHAR'/
! - For CALC_ERREUR
    data list_excl_ce/'EVOL_CHAR'/
! - For CALC_VECT_ELEM
    data list_excl_cve/'EVOL_CHAR'/
! - For CALC_MATR_ELEM
    data list_excl_cme/'EVOL_CHAR'/
! - For CALC_FORC_AJOU
    data list_excl_cfa/'EVOL_CHAR'/
! - For CALC_MATR_AJOU
    data list_excl_cma/'EVOL_CHAR'/
!
! --------------------------------------------------------------------------------------------------
!

!
! - Informations about loads
!
    if (present(list_load_)) then
        call load_list_info(load_empty, nb_load, v_load_name, v_load_info, &
                            list_load_=list_load_)
    elseif (present(lload_name_)) then
        call load_list_info(load_empty, nb_load, v_load_name, v_load_info, &
                            lload_name_, lload_info_)
    elseif (present(list_name_)) then
        nb_load = list_nbload_
    else
        ASSERT(.false.)
    end if
!
! - Look for excluded loads
!
    do i_load = 1, nb_load
        if (present(list_name_)) then
            load_name = list_name_(i_load) (1:8)
        else
            load_name = v_load_name(i_load) (1:8)
        end if
        call load_neut_iden(nb_type_neum, load_name, list_load_keyw)
        do i_type_neum = 1, nb_type_neum
            if (list_load_keyw(i_type_neum)) then
                call load_neut_data(i_type_neum, nb_type_neum, load_keyw_=load_keyw)
                if (command .eq. 'THER_NON_LINE_MO') then
                    do i_excl = 1, nb_excl_tnlm
                        if (load_keyw .eq. list_excl_tnlm(i_excl)) then
                            call utmess('F', 'CHARGES6_3')
                        end if
                    end do
                end if
                if (command .eq. 'CALC_VECT_ELEM') then
                    do i_excl = 1, nb_excl_cve
                        if (load_keyw .eq. list_excl_cve(i_excl)) then
                            call utmess('F', 'CHARGES6_3')
                        end if
                    end do
                end if
                if (command .eq. 'CALC_MATR_ELEM') then
                    do i_excl = 1, nb_excl_cme
                        if (load_keyw .eq. list_excl_cme(i_excl)) then
                            call utmess('F', 'CHARGES6_3')
                        end if
                    end do
                end if
                if (command .eq. 'CALC_FORC_AJOU') then
                    do i_excl = 1, nb_excl_cfa
                        if (load_keyw .eq. list_excl_cfa(i_excl)) then
                            call utmess('F', 'CHARGES6_3')
                        end if
                    end do
                end if
                if (command .eq. 'CALC_MATR_AJOU') then
                    do i_excl = 1, nb_excl_cma
                        if (load_keyw .eq. list_excl_cma(i_excl)) then
                            call utmess('F', 'CHARGES6_3')
                        end if
                    end do
                end if
                if (command .eq. 'CALC_ERREUR') then
                    do i_excl = 1, nb_excl_ce
                        if (load_keyw .eq. list_excl_ce(i_excl)) then
                            call utmess('F', 'CHARGES6_3')
                        end if
                    end do
                end if
            end if
        end do
    end do
end subroutine
