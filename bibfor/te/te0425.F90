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
subroutine te0425(option, nomte)
!
    use contact_module
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/apnorm.h"
#include "asterfort/assert.h"
#include "asterfort/dctest.h"
#include "asterfort/gapGetParamCoor.h"
#include "asterfort/jevech.h"
#include "asterfort/lcgeominit.h"
#include "asterfort/mmelem.h"
#include "asterfort/mmnewd.h"
#include "asterfort/projInsideCell.h"
#include "asterfort/reerel.h"
#include "jeveux.h"
#include "contact_module.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!  Compute Geometric Gap for Mortar methods
! --------------------------------------------------------------------------------------------------
!
    integer :: jgeom, jgap, jstat, index
    integer :: nnl, nbcps, nbdm, elem_dime, nb_dof
    integer :: nb_node_slav, nb_node_mast
    integer :: iret_, lin_sub(1, 4), lin_nbsub
    integer :: i_node, i_dime, lin_mast_nbnode(1)
    aster_logical :: laxis, leltf
    character(len=8) :: elem_slav_code, elem_mast_code, lin_mast_code
    real(kind=8) :: elem_mast_coor(3, 9), elem_slav_coor(3, 9), para_coor(2, 9)
    real(kind=8) :: tau1(3), tau2(3), iNodeCoorReal(3)
    real(kind=8) ::  ksi_ma(2), ksi_line(2), slav_norm(3), coor_ma_re(3)
    real(kind=8) :: dist, ksi1, ksi2, pair_tole, gap
!
    pair_tole = PROJ_TOLE
!
    call jevech('PGEOMCR', 'L', jgeom)
    call jevech('PVECGAP', 'E', jgap)
    call jevech('PVEIGAP', 'E', jstat)
!
! - Get parameters
!
    call mmelem(nomte, elem_dime, nb_dof, &
                elem_slav_code, nb_node_slav, &
                elem_mast_code, nb_node_mast, &
                nnl, nbcps, nbdm, laxis, leltf)
!
! - No values on master side
!
    zr(jgap-1+1:jgap-1+nb_node_slav) = 0.0
    zr(jstat-1+1:jstat-1+nb_node_slav) = 0.0
!
    if (elem_slav_code == "POI1") goto 999
!
! - Get coordinates
!
    elem_mast_coor = 0.d0
    elem_slav_coor = 0.d0
    index = 0
!
    do i_node = 1, nb_node_slav
        do i_dime = 1, elem_dime
            elem_slav_coor(i_dime, i_node) = zr(jgeom-1+index+i_dime)
        end do
        index = index+elem_dime
    end do
!
    do i_node = 1, nb_node_mast
        do i_dime = 1, elem_dime
            elem_mast_coor(i_dime, i_node) = zr(jgeom-1+index+i_dime)
        end do
        index = index+elem_dime
    end do
!
! - Get parametric slave coordinates
!
    call gapGetParamCoor(elem_slav_code, para_coor)
!
! - Linearize master cell
!
    call dctest(elem_mast_code, lin_sub, lin_mast_nbnode, lin_nbsub, lin_mast_code)
!
    iNodeCoorReal = 0.d0
!
    do i_node = 1, nb_node_slav
        ksi1 = 0.d0
        ksi2 = 0.d0
!
! ----- Get node coordinates in real space (ndim)
!
        iNodeCoorReal(1:elem_dime) = elem_slav_coor(1:elem_dime, i_node)
!
! ----- Compute slave normal
!
        ksi1 = para_coor(1, i_node)
        if (elem_dime .eq. 3) then
            ksi2 = para_coor(2, i_node)
        end if

        call apnorm(nb_node_slav, elem_slav_code, elem_dime, elem_slav_coor, &
                    ksi1, ksi2, slav_norm)
!
! ----- Projection of node on linear master cell to know
!                    if it projected inside master cell
!
        call mmnewd(lin_mast_code, lin_mast_nbnode(1), elem_dime, elem_mast_coor, &
                    iNodeCoorReal, 75, pair_tole, slav_norm, &
                    ksi_line(1), ksi_line(2), tau1, tau2, &
                    iret_)
        ASSERT(iret_ == 0)
        call projInsideCell(pair_tole, elem_dime, lin_mast_code, ksi_line, iret_)
!
        if (iret_ == 1) then
            go to 99
        end if
!
! ----- Projection of node on master cell
!
        call mmnewd(elem_mast_code, nb_node_mast, elem_dime, elem_mast_coor, &
                    iNodeCoorReal, 75, pair_tole, slav_norm, &
                    ksi_ma(1), ksi_ma(2), tau1, tau2, &
                    iret_, dist, ksi_line(1), ksi_line(2))
        ASSERT(iret_ == 0)
!
        call projInsideCell(pair_tole, elem_dime, lin_mast_code, ksi_ma, iret_)
        ASSERT(iret_ == 0)
!
! ----- Return in real master space
!
        coor_ma_re = 0.d0
        call reerel(elem_mast_code, nb_node_mast, 3, elem_mast_coor, ksi_ma, coor_ma_re)
!
        gap = gapEval(iNodeCoorReal, coor_ma_re, slav_norm)
!
! ----- Save distance and status
!
        zr(jgap-1+i_node) = gap
        zr(jstat-1+i_node) = 1.d0
!
99      continue
    end do
!
999 continue
!
end subroutine
