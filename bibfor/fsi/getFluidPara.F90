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
! person_in_charge: mickael.abbas at edf.fr
!
subroutine getFluidPara(j_mater, &
                        rho_, cele_r_, pesa_, alpha_, cele_i_, r_)
!
    implicit none
!
#include "asterfort/rcvalb.h"
!
! --------------------------------------------------------------------------------------------------
!
! Utilities for FSI
!
! Get material properties for fluid
!
! --------------------------------------------------------------------------------------------------
!
! In  j_mater          : coded material address
! Out rho              : density of fluid
! Out alpha            : reflection coefficient
! Out cele_i           : imaginary part of sound speed in fluid
! Out r                : geometric lenght for first order Sommerfeld condition
! Out cele_r           : real part of sound speed in fluid
! Out pesa             : gravity
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), intent(in) :: j_mater
    real(kind=8), optional, intent(out) :: rho_, cele_r_, pesa_, alpha_, cele_i_, r_
    integer(kind=8), parameter :: nb_resu = 6
    integer(kind=8) :: icodre(nb_resu)
    character(len=16), parameter :: resu_name(nb_resu) = &
                                    (/'RHO      ', 'COEF_AMOR', 'CELE_I   ', 'LONG_CARA', &
                                      'CELE_R   ', 'PESA_Z   '/)
    real(kind=8) :: resu_vale(nb_resu)
    real(kind=8) :: rho, cele_r, pesa, alpha, cele_i, r_impe
    character(len=8) :: fami
    character(len=1) :: poum
    integer(kind=8) :: ipg, ispg
!
! --------------------------------------------------------------------------------------------------
!
    fami = 'FPG1'
    ipg = 1
    ispg = 1
    poum = '+'
!
    if (present(pesa_)) then
        call rcvalb(fami, ipg, ispg, poum, j_mater, &
                    ' ', 'FLUIDE', 0, ' ', [0.d0], &
                    nb_resu, resu_name, resu_vale, icodre, 1)
        rho = resu_vale(1)
        alpha = resu_vale(2)
        cele_i = resu_vale(3)
        r_impe = resu_vale(4)
        cele_r = resu_vale(5)
        pesa = resu_vale(6)
    else
        if (present(cele_r_)) then
            call rcvalb(fami, ipg, ispg, poum, j_mater, &
                        ' ', 'FLUIDE', 0, ' ', [0.d0], &
                        5, resu_name, resu_vale, icodre, 1)
            rho = resu_vale(1)
            alpha = resu_vale(2)
            cele_i = resu_vale(3)
            r_impe = resu_vale(4)
            cele_r = resu_vale(5)
        else
            call rcvalb(fami, ipg, ispg, poum, j_mater, &
                        ' ', 'FLUIDE', 0, ' ', [0.d0], &
                        4, resu_name, resu_vale, icodre, 1)
            rho = resu_vale(1)
            alpha = resu_vale(2)
            cele_i = resu_vale(3)
            r_impe = resu_vale(4)
        end if
    end if
!
    if (present(rho_)) then
        rho_ = rho
    end if
    if (present(cele_r_)) then
        cele_r_ = cele_r
    end if
    if (present(pesa_)) then
        pesa_ = pesa
    end if
    if (present(alpha_)) then
        alpha_ = alpha
    end if
    if (present(cele_i_)) then
        cele_i_ = cele_i
    end if
    if (present(r_)) then
        r_ = r_impe
    end if

!
end subroutine
