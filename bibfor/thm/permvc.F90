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
subroutine permvc(ds_thm, satur, &
                  krl, dkrl_dsatur, krg, dkrg_dsatur)
!
    use THM_type
!
    implicit none
!
#include "asterfort/kfomvc.h"
#include "asterfort/regup1.h"
#include "asterfort/regup2.h"
!
    type(THM_DS), intent(in) :: ds_thm
    real(kind=8), intent(in) :: satur
    real(kind=8), intent(out) :: krl, dkrl_dsatur
    real(kind=8), intent(out) :: krg, dkrg_dsatur
!
! --------------------------------------------------------------------------------------------------
!
! THM - Permeability (HYDR_VGC)
!
! Evaluate permeability for liquid and gaz
!
! --------------------------------------------------------------------------------------------------
!
! In  ds_thm           : datastructure for THM
! In  satur            : saturation
! Out krl              : value of kr(liquid)
! Out dkrl_dsatur      : value of d(kr(liquid))/dSatur
! Out krg              : value of kr(gaz)
! Out dkrg_dsatur      : value of d(kr(gaz))/dSatur
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: s1, sr, satur_min, satur_max
    real(kind=8) :: n, m, usm
    real(kind=8) :: y0w, y0wp, y1, a1, b1, c1
    real(kind=8) :: ar, br
!
! --------------------------------------------------------------------------------------------------
!
    krl = 0.d0
    dkrl_dsatur = 0.d0
    krg = 0.d0
    dkrg_dsatur = 0.d0
!
! - Get parameters
!
    n = ds_thm%ds_material%hydr%n
    sr = ds_thm%ds_material%hydr%sr
    satur_max = ds_thm%ds_material%hydr%smax
    m = 1.d0-1.d0/n
    usm = 1.d0/m
    s1 = (satur-sr)/(1.d0-sr)
    satur_min = sr+(1.d0-sr)*(1.d0-satur_max)
!
    if ((satur .lt. satur_max) .and. (satur .gt. satur_min)) then
        call kfomvc(sr, m, usm, satur, &
                    krl, krg, dkrl_dsatur, dkrg_dsatur)
    else if (satur .ge. satur_max) then
        call kfomvc(sr, m, usm, satur_max, &
                    y0w, krg, y0wp, dkrg_dsatur)
        y1 = 1.d0
        call regup2(satur_max, y0w, y0wp, y1, a1, &
                    b1, c1)
        krl = a1*satur*satur+b1*satur+c1
        dkrl_dsatur = 2.d0*a1*satur+b1
    else if (satur .le. satur_min) then
        call kfomvc(sr, m, usm, satur_min, &
                    y0w, krg, y0wp, dkrg_dsatur)
        call regup1(satur_min, y0w, y0wp, ar, br)
        krl = ar*satur+br
        dkrl_dsatur = ar
    end if
!
end subroutine
