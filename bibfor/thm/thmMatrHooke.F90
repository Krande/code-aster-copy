! --------------------------------------------------------------------
! Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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
subroutine thmMatrHooke(ds_thm, angl_naut)
!
use THM_type
!
implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/matrHooke3d.h"
#include "asterfort/separ_RI_elas_3D.h"
!
type(THM_DS), intent(inout) :: ds_thm
real(kind=8), intent(in) :: angl_naut(3)
!
! --------------------------------------------------------------------------------------------------
!
! THM
!
! Compute Hooke elastic matrix
!
! --------------------------------------------------------------------------------------------------
!
! IO  ds_thm           : datastructure for THM
! In  angl_naut        : nautical angles
!                        (1) Alpha - clockwise around Z0
!                        (2) Beta  - counterclockwise around Y1
!                        (1) Gamma - clockwise around X
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: repere(7)
    real(kind=8) :: h(6), hi(6), g, e, nu
    real(kind=8) :: e1i, e2i, e3i, gi
    real(kind=8) :: nu12i, nu13i, nu23i, nui
!
! --------------------------------------------------------------------------------------------------
!
    repere(:) = 0.d0
    repere(1) = 1.d0
    repere(2) = angl_naut(1)
    repere(3) = angl_naut(2)
    repere(4) = angl_naut(3)
!
! - Prepare Hook matrix coefficient
!
    if (ds_thm%ds_material%elas%id .eq. 1) then
        e = ds_thm%ds_material%elas%e
        nu = ds_thm%ds_material%elas%nu
        g = 2*e*(1.d0 + nu)
    else
        g = ds_thm%ds_material%elas%g
    endif
    call separ_RI_elas_3D(ds_thm%ds_material%elas%id ,&
                          ds_thm%ds_material%elas%nu ,&
                          g, nui ,gi, &
                          ds_thm%ds_material%elas%e_l,&
                          ds_thm%ds_material%elas%e_t,&
                          ds_thm%ds_material%elas%e_n,&
                          ds_thm%ds_material%elas%nu_lt,&
                          ds_thm%ds_material%elas%nu_ln,&
                          ds_thm%ds_material%elas%nu_tn,&
                          e1i     , e2i  , e3i  ,&
                          nu12i   , nu13i, nu23i,&
                          h, hi)
!
! - Compute matrix
!
    call matrHooke3d(ds_thm%ds_material%elas%id, repere,&
                     h, g,&
                     g1 = ds_thm%ds_material%elas%g_lt,&
                     g2 = ds_thm%ds_material%elas%g_ln,&
                     g3 = ds_thm%ds_material%elas%g_tn,&
                     matr_elas = ds_thm%ds_material%elas%d)
!
end subroutine
