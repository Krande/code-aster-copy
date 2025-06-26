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

subroutine epthmc(fami, nno, ndim, nbsig, npg, &
                  shape_func, angl_naut, time, j_mater, &
                  option, epsi_varc)
!
    implicit none
!
#include "asterfort/epstmc.h"
#include "asterfort/lteatt.h"
!
!
    character(len=*), intent(in) :: fami
    integer(kind=8), intent(in) :: nno
    integer(kind=8), intent(in) :: ndim
    integer(kind=8), intent(in) :: nbsig
    integer(kind=8), intent(in) :: npg
    real(kind=8), intent(in) :: shape_func(1)
    real(kind=8), intent(in) :: angl_naut(3)
    real(kind=8), intent(in) :: time
    integer(kind=8), intent(in) :: j_mater
    character(len=16), intent(in) :: option
    real(kind=8), intent(out) :: epsi_varc(1)
!
! --------------------------------------------------------------------------------------------------
!
! Compute variable commands strains (thermics, drying, etc.)
!
! --------------------------------------------------------------------------------------------------
!
! In  fami         : Gauss family for integration point rule
! In  nno          : number of nodes
! In  ndim         : dimension of space
! In  nbsig        : number of stress tensor components
! In  npg          : number of Gauss points
! In  shape_func   : shape function
! In  j_mater      : coded material address
! In  repere       : definition of basis (for non-isotropic materials)
! In  time         : current time
! In  option       : name of option to compute
! Out epsi_varc    : command variables strains
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: epsi_ther(6), epsi_hydr(6), epsi_sech(6), epsi_anel(6), epsi_pres(6)
    character(len=16) :: optio2
    integer(kind=8) :: i, kpg, ndim2
    real(kind=8) :: zero
!
! --------------------------------------------------------------------------------------------------
!
    zero = 0.d0
    epsi_varc(1:nbsig*npg) = zero
    ndim2 = ndim
    if (lteatt('FOURIER', 'OUI')) then
        ndim2 = 2
    end if
!
! - Loop on Gauss points
!
    do kpg = 1, npg
!
! ----- Thermic strains
!
        optio2 = ' '
        call epstmc(fami, ndim, time, '+', kpg, &
                    1, angl_naut, j_mater, optio2, &
                    epsi_ther)
!
! ----- Hydric strains
!
        optio2 = option(1:9)//'_HYDR'
        call epstmc(fami, ndim, time, '+', kpg, &
                    1, angl_naut, j_mater, optio2, &
                    epsi_hydr)
!
! ----- Drying strains
!
        optio2 = option(1:9)//'_SECH'
        call epstmc(fami, ndim, time, '+', kpg, &
                    1, angl_naut, j_mater, optio2, &
                    epsi_sech)
!
! ----- Anelastic strains (given by user)
!
        optio2 = option(1:9)//'_EPSA'
        call epstmc(fami, ndim, time, '+', kpg, &
                    1, angl_naut, j_mater, optio2, &
                    epsi_anel)
!
! ----- Pressure strains
!
        optio2 = option(1:9)//'_PTOT'
        call epstmc(fami, ndim, time, '+', kpg, &
                    1, angl_naut, j_mater, optio2, &
                    epsi_pres)
!
! ----- Total command variables strains
!
        do i = 1, nbsig
            epsi_varc(i+nbsig*(kpg-1)) = epsi_varc(i+nbsig*(kpg-1))+ &
                                         epsi_ther(i)+ &
                                         epsi_hydr(i)+ &
                                         epsi_sech(i)+ &
                                         epsi_anel(i)+ &
                                         epsi_pres(i)
        end do
    end do
!
end subroutine
