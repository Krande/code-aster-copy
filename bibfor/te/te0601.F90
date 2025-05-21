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
subroutine te0601(option, nomte)

    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/jevech.h"
#include "jeveux.h"

!
    character(len=16) :: nomte, option
!
!
! --------------------------------------------------------------------------------------------------
!
! SHELL-3D link 
!
! Link elementary matrices
! --------------------------------------------------------------------------------------------------
!
    integer :: jv_geom
    integer :: nb_co_nodes, nb_3d_nodes, dim, index
    integer :: i_node_co, i_node_3d, i
    real(kind=8) :: coor(3)

    call jevech('PGEOMER', 'L', jv_geom)
    dim = 3
    index = 0

    do i_node_co = 1, 2
        do i= 1, dim
            coor(i) = zr(jv_geom-1+index+i)
        end do
        !write(*,*) "les coordonnées", coor(1), " ", coor(2), " ", coor(3)
        index = index + dim
    end do

    do i_node_3d = 1, 3
        do i= 1, dim
            coor(i) = zr(jv_geom-1+index+i)
        end do
        !write(*,*) "les coordonnées", coor(1), " ", coor(2), " ", coor(3)
        index = index + dim
    end do


end subroutine