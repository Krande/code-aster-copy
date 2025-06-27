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
subroutine niElem(nomte, geom, param)
!
    use contact_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/lteatt.h"
#include "asterfort/niQuantities.h"
#include "contact_module.h"
!
    character(len=16), intent(in) :: nomte
    type(ContactGeom), intent(inout) :: geom
    type(ContactParameters), intent(inout) :: param
!
! --------------------------------------------------------------------------------------------------
!
! Contact (Nitsche method) - Elementary computations
!
! Get informations about contact element
!
! --------------------------------------------------------------------------------------------------
!
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
!
    geom%l_axis = lteatt('AXIS', 'OUI')
!
    if (nomte(1:2) .ne. "CN" .and. nomte(1:2) .ne. "FN") then
        ASSERT(ASTER_FALSE)
    end if
!
! - Slave side
!
    select case (nomte(3:4))
    case ("Q4")
        geom%elem_dime = 2
        geom%elem_volu_code = 'QU4'
        geom%nb_node_volu = 4
        geom%elem_slav_code = 'SE2'
        geom%nb_node_slav = 2
    case ("Q8")
        geom%elem_dime = 2
        geom%elem_volu_code = 'QU8'
        geom%nb_node_volu = 8
        geom%elem_slav_code = 'SE3'
        geom%nb_node_slav = 3
    case ("Q9")
        geom%elem_dime = 2
        geom%elem_volu_code = 'QU9'
        geom%nb_node_volu = 9
        geom%elem_slav_code = 'SE3'
        geom%nb_node_slav = 3
    case ("T3")
        geom%elem_dime = 2
        geom%elem_volu_code = 'TR3'
        geom%nb_node_volu = 3
        geom%elem_slav_code = 'SE2'
        geom%nb_node_slav = 2
    case ("T6")
        geom%elem_dime = 2
        geom%elem_volu_code = 'TR6'
        geom%nb_node_volu = 6
        geom%elem_slav_code = 'SE3'
        geom%nb_node_slav = 3
    case default
        ASSERT(ASTER_FALSE)
    end select
!
! - Master side
!
    select case (nomte(5:6))
    case ("S2")
        ASSERT(geom%elem_dime == 2)
        geom%elem_mast_code = 'SE2'
        geom%nb_node_mast = 2
    case ("S3")
        ASSERT(geom%elem_dime == 2)
        geom%elem_mast_code = 'SE3'
        geom%nb_node_mast = 3
    case ("Q4")
        ASSERT(geom%elem_dime == 3)
        geom%elem_mast_code = 'QU4'
        geom%nb_node_mast = 4
    case ("Q8")
        ASSERT(geom%elem_dime == 3)
        geom%elem_mast_code = 'QU8'
        geom%nb_node_mast = 8
    case ("Q9")
        ASSERT(geom%elem_dime == 3)
        geom%elem_mast_code = 'QU9'
        geom%nb_node_mast = 9
    case ("T3")
        ASSERT(geom%elem_dime == 3)
        geom%elem_mast_code = 'TR3'
        geom%nb_node_mast = 3
    case ("T6")
        ASSERT(geom%elem_dime == 3)
        geom%elem_mast_code = 'TR6'
        geom%nb_node_mast = 6
    case default
        ASSERT(ASTER_FALSE)
    end select
!
    if (geom%l_axis) then
        if (nomte(7:7) .ne. "A") then
            ASSERT(ASTER_FALSE)
        end if
    end if
!
    geom%nb_dofs = (geom%nb_node_mast+geom%nb_node_volu)*geom%elem_dime
!
    ASSERT(geom%nb_node_slav .le. 27)
    ASSERT(geom%nb_node_mast .le. 9)
    ASSERT(geom%nb_lagr_c == 0)
    ASSERT(geom%nb_dofs .le. MAX_NITS_DOFS)
    ASSERT((geom%elem_dime .eq. 2) .or. (geom%elem_dime .eq. 3))
!
    call niQuantities(geom, param)
!
end subroutine
