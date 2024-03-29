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
!
module contact_type
!
!
    implicit none
!
    private
!
#include "asterf_types.h"
#include "contact_module.h"
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Type
!
! Generic type for contact
!
! --------------------------------------------------------------------------------------------------
! --------------------------------------------------------------------------------------------------
!
    type ContactParameters
    !! Contact parameters
        integer                             :: algo_cont = 0
        integer                             :: type_cont = 0
        integer                             :: vari_cont = 0
        real(kind=8)                        :: vari_cont_coef = 0.d0
        real(kind=8), dimension(9)          :: coef_cont = 0.d0

    !! Friction paramaters
        aster_logical                       :: l_fric = ASTER_FALSE
        integer                             :: algo_fric = 0
        integer                             :: type_fric = 0
        real(kind=8), dimension(9)          :: coef_fric = 0.d0
        real(kind=8), dimension(9)          :: threshold = 0.d0
        real(kind=8)                        :: threshold_given = 0.d0

    !! Other
        real(kind=8)                        :: proj_tole = 0.d0
        integer                             :: cont_init = PAIR_CONT_INTE
    end type
!
    type ContactGeom
    !! Slave side parameters
        integer                             :: nb_node_slav = 0
        character(len=8)                    :: elem_slav_code = " "
        real(kind=8), dimension(3, 9)        :: coor_slav_init = 0.d0
        real(kind=8), dimension(3, 9)        :: coor_slav_prev = 0.d0
        real(kind=8), dimension(3, 9)        :: coor_slav_curr = 0.d0
        real(kind=8), dimension(3, 9)        :: coor_slav_pair = 0.d0
        real(kind=8), dimension(3, 9)        :: depl_slav_curr = 0.d0
        real(kind=8), dimension(4)          :: lagc_slav_curr = 0.d0
        real(kind=8), dimension(2, 4)        :: lagf_slav_curr = 0.d0
        integer                             :: nb_lagr_c = 0
        integer, dimension(9)               :: indi_lagc = 0

    !! Slave cell volu
        integer                             :: nb_node_volu = 0
        character(len=8)                    :: elem_volu_code = " "
        real(kind=8), dimension(3, 27)        :: coor_volu_init = 0.d0
        real(kind=8), dimension(3, 27)        :: coor_volu_prev = 0.d0
        real(kind=8), dimension(3, 27)        :: coor_volu_curr = 0.d0
        real(kind=8), dimension(3, 27)        :: coor_volu_pair = 0.d0
        real(kind=8), dimension(3, 27)        :: depl_volu_curr = 0.d0
        integer, dimension(9)                :: mapVolu2Surf = 0

    !! Master side paramaters
        integer                             :: nb_node_mast = 0
        character(len=8)                    :: elem_mast_code = " "
        real(kind=8), dimension(3, 9)        :: coor_mast_init = 0.d0
        real(kind=8), dimension(3, 9)        :: coor_mast_prev = 0.d0
        real(kind=8), dimension(3, 9)        :: coor_mast_curr = 0.d0
        real(kind=8), dimension(3, 9)        :: coor_mast_pair = 0.d0
        real(kind=8), dimension(3, 9)        :: depl_mast_curr = 0.d0

    !! Time
        real(kind=8) :: time_prev = 0.d0, time_curr = 0.d0

    !! Other
        integer                             :: elem_dime = 0
        aster_logical                       :: l_axis = ASTER_FALSE
        integer                             :: nb_dofs = 0
    end type
!
!===================================================================================================
!
!===================================================================================================
!
    public :: ContactParameters, ContactGeom
!
contains
!
!===================================================================================================
!
!===================================================================================================
!
end module
