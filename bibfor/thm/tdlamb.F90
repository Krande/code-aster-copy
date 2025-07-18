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
subroutine tdlamb(ds_thm, angl_naut, ndim, tdlamt)
!
    use THM_type
!
    implicit none
!
#include "asterfort/matrot.h"
#include "asterfort/utbtab.h"
#include "asterfort/THM_type.h"
!
    type(THM_DS), intent(in) :: ds_thm
    real(kind=8), intent(in) :: angl_naut(3)
    integer(kind=8), intent(in) :: ndim
    real(kind=8), intent(out) :: tdlamt(ndim, ndim)
!
! --------------------------------------------------------------------------------------------------
!
! THM
!
! Compute tensor of derivatives (by temperature) for thermal conductivity
!
! --------------------------------------------------------------------------------------------------
!
! In  ds_thm           : datastructure for THM
! In  angl_naut        : nautical angles
!                        (1) Alpha - clockwise around Z0
!                        (2) Beta  - counterclockwise around Y1
!                        (1) Gamma - clockwise around X
! In  ndim             : dimension of space
! Out tdlamt           : tensor of derivatives for thermal conductivity
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: tdlamti(3, 3), passag(3, 3), work(3, 3), tk2(3, 3)
!
! --------------------------------------------------------------------------------------------------
!
    work(:, :) = 0.d0
    passag(:, :) = 0.d0
    tdlamti(:, :) = 0.d0
    tk2(:, :) = 0.d0
    tdlamt(:, :) = 0.d0
!
    if (ds_thm%ds_material%ther%cond_type .eq. THER_COND_ISOT) then
        tdlamt(1, 1) = ds_thm%ds_material%ther%dlambda
        tdlamt(2, 2) = ds_thm%ds_material%ther%dlambda
        if (ndim .eq. 3) then
            tdlamt(3, 3) = ds_thm%ds_material%ther%dlambda
        end if
    else if (ds_thm%ds_material%ther%cond_type .eq. THER_COND_ISTR) then
        if (ndim .eq. 3) then
            tdlamti(1, 1) = ds_thm%ds_material%ther%dlambda_tl
            tdlamti(2, 2) = ds_thm%ds_material%ther%dlambda_tl
            tdlamti(3, 3) = ds_thm%ds_material%ther%dlambda_tn
            call matrot(angl_naut, passag)
            call utbtab('ZERO', 3, 3, tdlamti, passag, work, tk2)
            tdlamt = tk2
        end if
    else if (ds_thm%ds_material%ther%cond_type .eq. THER_COND_ORTH) then
        tdlamti(1, 1) = ds_thm%ds_material%ther%dlambda_tl
        tdlamti(2, 2) = ds_thm%ds_material%ther%dlambda_tt
        call matrot(angl_naut, passag)
        call utbtab('ZERO', 3, 3, tdlamti, passag, work, tk2)
        tdlamt(1, 1) = tk2(1, 1)
        tdlamt(2, 2) = tk2(2, 2)
        tdlamt(1, 2) = tk2(1, 2)
        tdlamt(2, 1) = tk2(2, 1)
    else
! ----- No thermic
    end if
!
end subroutine
