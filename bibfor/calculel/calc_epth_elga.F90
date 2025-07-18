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

subroutine calc_epth_elga(fami, ndim, poum, kpg, ksp, &
                          j_mater, angl_naut, epsi_ther)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/matrot.h"
#include "asterfort/get_elas_id.h"
#include "asterfort/utpslg.h"
#include "asterfort/verift.h"
!
!
    character(len=*), intent(in) :: fami
    integer(kind=8), intent(in) :: ndim
    character(len=*), intent(in) :: poum
    integer(kind=8), intent(in) :: kpg
    integer(kind=8), intent(in) :: ksp
    integer(kind=8), intent(in) :: j_mater
    real(kind=8), intent(in) :: angl_naut(3)
    real(kind=8), intent(out) :: epsi_ther(6)
!
! --------------------------------------------------------------------------------------------------
!
! Compute thermic strains
!
! For isoparametric elements
!
! --------------------------------------------------------------------------------------------------
!
! In  fami         : Gauss family for integration point rule
! In  ndim         : dimension of space
! In  poum         : parameters evaluation
!                     '-' for previous temperature
!                     '+' for current temperature
!                     'T' for current and previous temperature
! In  kpg          : current point gauss
! In  ksp          : current "sous-point" gauss
! In  j_mater      : coded material address
! In  angl_naut    : nautical angles (for non-isotropic materials)
! Out epsi_ther    : thermal strain tensor
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16) :: elas_keyword
    integer(kind=8) :: elas_id
    real(kind=8) :: p_glob_loca(3, 3), epsi_ther_vect(6)
    real(kind=8) :: epsth, epsth_anis(3)
    real(kind=8) :: vepst1(6), vepst2(6)
!
! --------------------------------------------------------------------------------------------------
!
    epsi_ther(1:6) = 0.d0
!
! - Get elasticity type
!
    call get_elas_id(j_mater, elas_id, elas_keyword)
    ASSERT(elas_id .le. 3)
!
! - Non-isotropic elasticity: prepare basis
!
    if (elas_id .gt. 1) then
        call matrot(angl_naut, p_glob_loca)
    end if
!
! - Compute (local) thermic strains
!
    if (elas_id .eq. 1) then
        if (elas_keyword .eq. 'ELAS_META') then
            call verift(fami, kpg, ksp, poum, j_mater, &
                        epsth_meta_=epsth)
        else
            call verift(fami, kpg, ksp, poum, j_mater, &
                        epsth_=epsth)
        end if
        epsi_ther(1) = epsth
        epsi_ther(2) = epsth
        epsi_ther(3) = epsth
    else if (elas_id .eq. 2) then
        call verift(fami, kpg, ksp, poum, j_mater, &
                    epsth_anis_=epsth_anis)
        epsi_ther_vect(1) = epsth_anis(1)
        epsi_ther_vect(2) = epsth_anis(2)
        epsi_ther_vect(3) = epsth_anis(3)
        epsi_ther_vect(4) = 0.d0
        epsi_ther_vect(5) = 0.d0
        epsi_ther_vect(6) = 0.d0
    else if (elas_id .eq. 3) then
        call verift(fami, kpg, ksp, poum, j_mater, &
                    epsth_anis_=epsth_anis)
        epsi_ther_vect(1) = epsth_anis(1)
        epsi_ther_vect(2) = epsth_anis(1)
        epsi_ther_vect(3) = epsth_anis(2)
        epsi_ther_vect(4) = 0.d0
        epsi_ther_vect(5) = 0.d0
        epsi_ther_vect(6) = 0.d0
    else
        ASSERT(.false.)
    end if
!
! - Non-isotropic elasticity: rotate strains
!
    if (elas_id .gt. 1) then
        vepst1(1) = epsi_ther_vect(1)
        vepst1(2) = epsi_ther_vect(4)
        vepst1(3) = epsi_ther_vect(2)
        vepst1(4) = epsi_ther_vect(5)
        vepst1(5) = epsi_ther_vect(6)
        vepst1(6) = epsi_ther_vect(3)
        call utpslg(1, 3, p_glob_loca, vepst1, vepst2)
        epsi_ther(1) = vepst2(1)
        epsi_ther(2) = vepst2(3)
        epsi_ther(3) = vepst2(6)
        epsi_ther(4) = vepst2(2)
        epsi_ther(5) = vepst2(4)
        epsi_ther(6) = vepst2(5)
        if (ndim .eq. 2) epsi_ther(3) = epsi_ther_vect(3)
    end if
!
end subroutine
