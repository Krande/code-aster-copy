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
subroutine thmGetElemDime(ndim, nnos, nnom, &
                          mecani, press1, press2, tempe, second, &
                          nddls, nddlm, &
                          nddl_meca, nddl_p1, nddl_p2, nddl_2nd, &
                          dimdep, dimdef, dimcon, dimuel)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/thmGetGeneDime.h"
!
    integer(kind=8), intent(in) :: ndim, nnos, nnom
    integer(kind=8), intent(in) :: mecani(5), press1(7), press2(7), tempe(5), second(5)
    integer(kind=8), intent(out) :: nddls, nddlm
    integer(kind=8), intent(out) :: nddl_meca, nddl_p1, nddl_p2, nddl_2nd
    integer(kind=8), intent(out) :: dimdep, dimdef, dimcon, dimuel
!
! --------------------------------------------------------------------------------------------------
!
! THM - Parameters
!
! Get dimensions about element
!
! --------------------------------------------------------------------------------------------------
!
! In  ndim             : dimension of element (2 ou 3)
! In  nnos             : number of nodes (not middle ones)
! In  nnom             : number of nodes (middle ones)
! In  mecani           : parameters for mechanic
! In  press1           : parameters for hydraulic (capillary pressure)
! In  press2           : parameters for hydraulic (gaz pressure)
! In  tempe            : parameters for thermic
! In  second           : parameters for second gradient
! Out nddls            : number of dof at nodes (not middle ones)
! Out nddlm            : number of dof at nodes (middle ones)
! Out nddl_meca        : number of dof for mechanical quantity
! Out nddl_p1          : number of dof for first hydraulic quantity
! Out nddl_p2          : number of dof for second hydraulic quantity
! Out nddl_2nd         : number of dof for second gradient
! Out dimdep           : dimension of generalized displacement vector
! Out dimdef           : dimension of generalized strains vector
! Out dimcon           : dimension of generalized stresses vector
! Out dimuel           : total number of dof for element
!
! --------------------------------------------------------------------------------------------------
!
    dimuel = 0
    dimdep = 0
    dimdef = 0
    dimcon = 0
    nddl_meca = 0
    nddl_p1 = 0
    nddl_p2 = 0
    nddl_2nd = 0
    nddls = 0
    nddlm = 0
    if (mecani(1) .eq. 1) then
        nddl_meca = ndim
    end if
    if (press1(1) .eq. 1) then
        nddl_p1 = 1
    end if
    if (press2(1) .eq. 1) then
        nddl_p2 = 1
    end if
    if (second(1) .eq. 1) then
        nddl_2nd = 2
    end if
!
! - Get dimensions of generalized vectors
!
    call thmGetGeneDime(ndim, &
                        mecani, press1, press2, tempe, second, &
                        dimdep, dimdef, dimcon)
!
! - Count dof
!
    nddls = mecani(1)*nddl_meca+press1(1)+press2(1)+tempe(1)+second(1)*nddl_2nd
    nddlm = mecani(1)*nddl_meca
    dimuel = nnos*nddls+nnom*nddlm
!
end subroutine
