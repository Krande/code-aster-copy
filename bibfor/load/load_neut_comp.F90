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

subroutine load_neut_comp(type_calc, stop_calc, model, time_curr, time, &
                          load_name, load_nume, nb_in_maxi, nb_in_prep, lpain, &
                          lchin, base, resu_elem, matr_vect_elem, l_stat, &
                          time_move_, i_load_)
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/calcul.h"
#include "asterfort/reajre.h"
#include "asterfort/gcnco2.h"
#include "asterfort/corich.h"
#include "asterfort/ntdomt.h"
#include "asterfort/multResuElem.h"
#include "asterfort/load_neut_spec.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    character(len=4), intent(in) :: type_calc
    character(len=1), intent(in) :: stop_calc
    character(len=24), intent(in) :: model
    real(kind=8), intent(in) :: time_curr
    character(len=24), intent(in) :: time
    character(len=8), intent(in) :: load_name
    integer, intent(in) :: load_nume
    integer, intent(in) :: nb_in_maxi
    integer, intent(in) :: nb_in_prep
    character(len=*), intent(inout) :: lpain(nb_in_maxi)
    character(len=*), intent(inout) :: lchin(nb_in_maxi)
    character(len=19), intent(inout) :: resu_elem
    character(len=19), intent(in) :: matr_vect_elem
    character(len=1), intent(in) :: base
    aster_logical, intent(in) :: l_stat
    character(len=24), optional, intent(in) :: time_move_
    integer, optional, intent(in) :: i_load_
!
! --------------------------------------------------------------------------------------------------
!
! Neumann loads computation - Thermic
!
! Elementary (on one load)
!
! --------------------------------------------------------------------------------------------------
!
! In  type_calc        : type of option to compute
!                        '2MBR' for second member (vector)
!                        'RESI' for residual (vector)
!                        'MRIG' for rigidity (matrix)
!                        'MTAN' for tangent matrix
! In  stop_calc        : CALCUL subroutine comportement
! In  model            : name of the model
! In  time_curr        : current time
! In  time             : time (<CARTE>)
! In  load_name        : name of current load
! In  load_nume        : identification of load type
! In  nb_in_maxi       : maximum number of input fields
! In  nb_in_prep       : number of input fields before specific ones
! IO  lpain            : list of input parameters
! IO  lchin            : list of input fields
! IO  resu_elem        : name of resu_elem
! In  matr_vect_elem   : name of matr_elem or vect_elem
! In  base             : JEVEUX base to create vect_elem
! In  time_move        : modified time (<CARTE>) for THER_NON_LINE_MO
! In  i_load           : index of current load
!
! --------------------------------------------------------------------------------------------------
!
    integer :: nb_type_neum
    parameter(nb_type_neum=11)
!
    integer :: i_type_neum, nb_in_add
    character(len=16) :: load_option
    character(len=24) :: load_ligrel
    integer :: nbout, nbin
    character(len=8) :: lpaout, newnom
    real(kind=8) :: theta
!
! --------------------------------------------------------------------------------------------------
!
    call ntdomt(theta)
!
    do i_type_neum = 1, nb_type_neum
!
! ----- Get information about load
!
        if (present(time_move_)) then
            call load_neut_spec('MOVE', type_calc, model, time_curr, time, &
                                load_name, load_nume, i_type_neum, nb_type_neum, nb_in_maxi, &
                                nb_in_prep, lchin, lpain, nb_in_add, lpaout, &
                                load_ligrel, load_option, &
                                time_move_=time_move_)
        else
            call load_neut_spec('STAT', type_calc, model, time_curr, time, &
                                load_name, load_nume, i_type_neum, nb_type_neum, nb_in_maxi, &
                                nb_in_prep, lchin, lpain, nb_in_add, lpaout, &
                                load_ligrel, load_option)
        end if
!
        if (load_option .ne. 'No_Load') then
!
! --------- Generate new RESU_ELEM name
!
            newnom = resu_elem(9:16)
            call gcnco2(newnom)
            resu_elem(10:16) = newnom(2:8)
!
! --------- Attach load to RESU_ELEM
!
            if (present(i_load_)) then
                call corich('E', resu_elem, ichin_=i_load_)
            else
                call corich('E', resu_elem, ichin_=-1)
            end if
!
! --------- Number of fields
!
            nbin = nb_in_prep+nb_in_add
            nbout = 1
!
! --------- Computation
!
            call calcul(stop_calc, load_option, load_ligrel, nbin, lchin, &
                        lpain, nbout, resu_elem, lpaout, base, &
                        'OUI')

            if (type_calc .ne. "2MBR" .and. .not. l_stat) then
                call multResuElem(resu_elem, theta)
            end if
!
! --------- Copying output field
!
            call reajre(matr_vect_elem, resu_elem, base)
        end if
    end do

end subroutine
