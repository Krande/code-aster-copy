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
subroutine calc_epth_elga(fami, ndim, poum, kpg, ksp, &
                          materCodeJv, anglNaut, epsiTher)
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/ElasticityMaterial_type.h"
#include "asterfort/get_elas_id.h"
#include "asterfort/matrot.h"
#include "asterfort/utmess.h"
#include "asterfort/utpslg.h"
#include "asterfort/verift.h"
#include "jeveux.h"
!
    character(len=*), intent(in) :: fami
    integer(kind=8), intent(in) :: ndim
    character(len=*), intent(in) :: poum
    integer(kind=8), intent(in) :: kpg, ksp
    integer(kind=8), intent(in) :: materCodeJv
    real(kind=8), intent(in) :: anglNaut(3)
    real(kind=8), intent(out) :: epsiTher(6)
!
! --------------------------------------------------------------------------------------------------
!
! Compute thermic strains
!
! For isoparametric elements
!
! --------------------------------------------------------------------------------------------------
!
! In  fami             : Gauss family for integration point rule
! In  ndim             : dimension of space
! In  poum             : parameters evaluation
!                     '-' for previous temperature
!                     '+' for current temperature
!                     'T' for current and previous temperature
! In  kpg              : current point gauss
! In  ksp              : current "sous-point" gauss
! In  materCodeJv      : coded material address
! In  anglNaut         : nautical angles for definition of basis for non-isotropic elasticity
! Out epsiTher         : thermal strain tensor
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: pgl(3, 3), epsiTherLoca(6)
    real(kind=8) :: epsiTherIsot, epsiTherAnis(3)
    real(kind=8) :: vepst1(6), vepst2(6)
    integer(kind=8) :: elasID
    character(len=32) :: elasKeyword
!
! --------------------------------------------------------------------------------------------------
!
    epsiTher = 0.d0
    epsiTherLoca = 0.d0

! - Get elasticity type
    call get_elas_id(materCodeJv, elasID, elasKeyword)
    if (elasID .gt. ELAS_ISTR) then
        call utmess("F", "COMPOR5_16", sk=elasKeyword)
    end if

! - Non-isotropic elasticity: prepare basis
    if (elasID .gt. ELAS_ISOT) then
        call matrot(anglNaut, pgl)
    end if

! - Compute (local) thermic strains
    if (elasID .eq. ELAS_ISOT) then
        if (elasKeyword .eq. 'ELAS_META') then
            call verift(fami, kpg, ksp, poum, materCodeJv, &
                        epsth_meta_=epsiTherIsot)
        else
            call verift(fami, kpg, ksp, poum, materCodeJv, &
                        epsth_=epsiTherIsot)
        end if
        epsiTher(1) = epsiTherIsot
        epsiTher(2) = epsiTherIsot
        epsiTher(3) = epsiTherIsot

    else if (elasID .eq. ELAS_ORTH) then
        call verift(fami, kpg, ksp, poum, materCodeJv, &
                    epsth_anis_=epsiTherAnis)
        epsiTherLoca(1) = epsiTherAnis(1)
        epsiTherLoca(2) = epsiTherAnis(2)
        epsiTherLoca(3) = epsiTherAnis(3)
        epsiTherLoca(4) = 0.d0
        epsiTherLoca(5) = 0.d0
        epsiTherLoca(6) = 0.d0

    else if (elasID .eq. ELAS_ISTR) then
        call verift(fami, kpg, ksp, poum, materCodeJv, &
                    epsth_anis_=epsiTherAnis)
        epsiTherLoca(1) = epsiTherAnis(1)
        epsiTherLoca(2) = epsiTherAnis(1)
        epsiTherLoca(3) = epsiTherAnis(2)
        epsiTherLoca(4) = 0.d0
        epsiTherLoca(5) = 0.d0
        epsiTherLoca(6) = 0.d0

    else
        call utmess("F", "COMPOR5_16", sk=elasKeyword)

    end if

! - Non-isotropic elasticity: rotate strains
    if (elasID .gt. ELAS_ISOT) then
        vepst1(1) = epsiTherLoca(1)
        vepst1(2) = epsiTherLoca(4)
        vepst1(3) = epsiTherLoca(2)
        vepst1(4) = epsiTherLoca(5)
        vepst1(5) = epsiTherLoca(6)
        vepst1(6) = epsiTherLoca(3)
        call utpslg(1, 3, pgl, vepst1, vepst2)
        epsiTher(1) = vepst2(1)
        epsiTher(2) = vepst2(3)
        epsiTher(3) = vepst2(6)
        epsiTher(4) = vepst2(2)
        epsiTher(5) = vepst2(4)
        epsiTher(6) = vepst2(5)
        if (ndim .eq. 2) then
            epsiTher(3) = epsiTherLoca(3)
        end if
    end if
!
end subroutine
