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
subroutine leverettIsotTher(c, temp, imate, hygr, dpc_, poro_, t0_C_, beta_, pc_)
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/rcvala.h"
#include "asterfort/leverettIsot.h"
!
!   .................................................................
!   evaluation of hygrometry with leverett isotherm for thermic
!
!   c (in) : concentration en eau
!   temp (in) : temperature
!
!   hygr (out) : hygrometry
!   dpc_ (out optional) : capillar pressure derivative
!   poro_ (out optional) : porosity
!   t0_C_ (out optional) : parameter TEMP_0_C
!   beta_ (out optional) : parameter VG_N
!   pc_ (out optional)  : capillar pressure
!
!   .................................................................
    integer(kind=8), intent(in) :: imate
    real(kind=8), intent(in) :: c, temp
    real(kind=8), intent(out) :: hygr
    real(kind=8), intent(out), optional :: dpc_, beta_, poro_, t0_C_, pc_
!
    integer(kind=8), parameter :: nbres = 6
    integer(kind=8)           :: codret(nbres)
    real(kind=8)      :: valres(nbres)
    character(len=16) :: nomres(nbres)
    real(kind=8)      :: alpha, ad, satu, dpc, beta, poro, t0_C, pc, coef_p
!
    nomres(1) = 'VG_PR'
    nomres(2) = 'VG_N'
    nomres(3) = 'ATH'
    nomres(4) = 'TEMP_0_C'
    nomres(5) = 'PORO'
    nomres(6) = "COEF_UNITE_P"
!
    call rcvala(imate, ' ', 'BETON_DESORP', 0, ' ', [0.d0], &
                nbres, nomres, valres, codret, 0)
    ASSERT(codret(1) .eq. 0)
    alpha = valres(1)
    beta = valres(2)
    ad = valres(3)
    t0_C = valres(4)
    poro = valres(5)
    coef_p = valres(6)
    satu = c/poro/1.d3
!
    call leverettIsot(temp, satu, alpha, beta, ad, t0_C, hygr, dpc, pc_=pc)

    if (present(dpc_)) then
        dpc_ = dpc
        beta_ = beta
        poro_ = poro
        t0_C_ = t0_C
    end if

    if (present(pc_)) then
        pc_ = pc/coef_p
    end if

end subroutine leverettIsotTher
