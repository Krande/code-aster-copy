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
!
subroutine laelem(nomte         , elem_dime     ,&
                  l_axis        , &
                  nb_dof        , nb_lagr_c       , indi_lagc   ,&
                  elem_slav_code, nb_node_slav,&
                  elem_mast_code, nb_node_mast)
!
implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/lteatt.h"
!
character(len=16), intent(in) :: nomte
integer, intent(out) :: elem_dime
aster_logical, intent(out) :: l_axis
integer, intent(out) :: nb_dof
integer, intent(out) :: nb_lagr_c
integer, intent(out) :: indi_lagc(9)
character(len=8), intent(out) :: elem_slav_code
integer, intent(out) :: nb_node_slav
character(len=8), intent(out) :: elem_mast_code
integer, intent(out) :: nb_node_mast
!
! --------------------------------------------------------------------------------------------------
!
! Contact (augmented Lagrangian method) - Elementary computations
!
! Get informations about contact element
!
! --------------------------------------------------------------------------------------------------
!
! In  nomte            : type of finite element
! Out elem_dime        : dimension of elements
! Out l_axis           : .true. for axisymmetric element
! Out nb_dof           : total number of dof on contact element
! Out nb_lagr_c          : total number of Lagrangian dof on contact element
! Out indi_lagc        : node where Lagrangian dof is present (1) or not (0)
! Out elem_slav_code   : code element for slave side from contact element
! Out nb_node_slav     : number of nodes of for slave side from contact element
! Out elem_mast_code   : code element for master side from contact element
! Out nb_node_mast     : number of nodes of for master side from contact element
!
! --------------------------------------------------------------------------------------------------
!
    l_axis          = lteatt('AXIS','OUI')
    elem_dime       = 0
    indi_lagc       = 0
    nb_dof          = 0
    nb_lagr_c         = 0
    elem_slav_code  = ' '
    nb_node_slav    = 0
    elem_mast_code  = ' '
    nb_node_mast    = 0
!
    if(nomte(1:2) .ne. "CM") then
        ASSERT(ASTER_FALSE)
    end if
!
! - Slave side
!
    select case (nomte(3:4))
        case ("S2")
            elem_dime      = 2
            elem_slav_code = 'SE2'
            nb_node_slav   = 2
            nb_lagr_c      = 2
        case ("S3")
            elem_dime      = 2
            elem_slav_code = 'SE3'
            nb_node_slav   = 3
            nb_lagr_c      = 3
        case ("Q4")
            elem_dime      = 3
            elem_slav_code = 'QU4'
            nb_node_slav   = 4
            nb_lagr_c      = 4
        case ("Q8")
            elem_dime      = 3
            elem_slav_code = 'QU8'
            nb_node_slav   = 8
            nb_lagr_c      = 4
        case ("Q9")
            elem_dime      = 3
            elem_slav_code = 'QU9'
            nb_node_slav   = 9
            nb_lagr_c      = 4
        case ("T3")
            elem_dime      = 3
            elem_slav_code = 'TR3'
            nb_node_slav   = 3
            nb_lagr_c      = 3
        case ("T6")
            elem_dime      = 3
            elem_slav_code = 'TR6'
            nb_node_slav   = 6
            nb_lagr_c      = 3
        case ("P1")
            elem_slav_code = 'PO1'
            nb_node_slav   = 1
        case default
            ASSERT(ASTER_FALSE)
    end select
!
! - Master side
!
    select case (nomte(5:6))
        case ("S2")
            ASSERT(elem_dime == 2)
            elem_mast_code = 'SE2'
            nb_node_mast   = 2
        case ("S3")
            ASSERT(elem_dime == 2)
            elem_mast_code = 'SE3'
            nb_node_mast   = 3
        case ("Q4")
            ASSERT(elem_dime == 3)
            elem_mast_code = 'QU4'
            nb_node_mast   = 4
        case ("Q8")
            ASSERT(elem_dime == 3)
            elem_mast_code = 'QU8'
            nb_node_mast   = 8
        case ("Q9")
            ASSERT(elem_dime == 3)
            elem_mast_code = 'QU9'
            nb_node_mast   = 9
        case ("T3")
            ASSERT(elem_dime == 3)
            elem_mast_code = 'TR3'
            nb_node_mast   = 3
        case ("T6")
            ASSERT(elem_dime == 3)
            elem_mast_code = 'TR6'
            nb_node_mast   = 6
        case ("L2")
            elem_dime      = 2
            elem_mast_code = 'LAGR'
            nb_node_mast   = 0
            nb_lagr_c      = 1
        case ("N2")
            elem_dime      = 2
            elem_mast_code = 'NOLAGR'
            nb_node_mast   = 0
            nb_lagr_c      = 0
        case ("L3")
            elem_dime      = 3
            elem_mast_code = 'LAGR'
            nb_node_mast   = 0
            nb_lagr_c      = 1
        case ("N3")
            elem_dime      = 3
            elem_mast_code = 'NOLAGR'
            nb_node_mast   = 0
            nb_lagr_c      = 0
        case default
            ASSERT(ASTER_FALSE)
    end select
!
    if(l_axis) then
        if(nomte(7:7) .ne. "A") then
            ASSERT(ASTER_FALSE)
        end if
    end if
!
    indi_lagc(1:nb_lagr_c) = 1
    nb_dof = nb_node_mast*elem_dime + nb_node_slav*elem_dime+nb_lagr_c
!
    ASSERT(nb_node_slav .le. 9)
    ASSERT(nb_node_mast .le. 9)
    ASSERT(nb_lagr_c .le. 4)
    ASSERT(nb_dof .le. 58)
    ASSERT((elem_dime .eq. 2).or.(elem_dime .eq. 3))
!
end subroutine
