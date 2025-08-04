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
#include "asterf_types.h"
!
interface
    subroutine calcCalcMeca(nb_option, list_option, &
                            l_elem_nonl, &
                            listLoadZ, modelZ, caraElemZ, &
                            ds_constitutive, ds_material, ds_system, &
                            hval_incr, hval_algo, &
                            vediri, vefnod, &
                            nb_obje_maxi, obje_name, obje_sdname, nb_obje, &
                            l_pred)
        use NonLin_Datastructure_type
        integer(kind=8), intent(in) :: nb_option
        character(len=16), intent(in) :: list_option(:)
        aster_logical, intent(in) :: l_elem_nonl
        character(len=*), intent(in) :: listLoadZ, modelZ, caraElemZ
        type(NL_DS_Constitutive), intent(in) :: ds_constitutive
        type(NL_DS_Material), intent(in) :: ds_material
        type(NL_DS_System), intent(in) :: ds_system
        character(len=19), intent(in) :: hval_incr(:), hval_algo(:)
        character(len=19), intent(in) :: vediri, vefnod
        integer(kind=8), intent(in) :: nb_obje_maxi
        character(len=16), intent(inout) :: obje_name(nb_obje_maxi)
        character(len=24), intent(inout) :: obje_sdname(nb_obje_maxi)
        integer(kind=8), intent(out) ::  nb_obje
        aster_logical, intent(in) :: l_pred
    end subroutine calcCalcMeca
end interface
