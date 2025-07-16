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
subroutine leverettIsotMeca(fami, kpg, ksp, imate, hygr_prev, hygr_curr)
    implicit none
#include "asterfort/assert.h"
#include "asterfort/leverettIsot.h"
#include "asterfort/rcvala.h"
#include "asterfort/rcvarc.h"
#include "asterfort/utmess.h"
!
!.................................................................
!   evaluation of hygrometry with leverett isotherm
!.................................................................
    character(len=*), intent(in) :: fami
    integer(kind=8), intent(in) :: kpg, ksp, imate
    real(kind=8), intent(out) :: hygr_prev, hygr_curr
!
    integer(kind=8)           :: codret(5), iret
    real(kind=8)      :: valres(5)
    character(len=16) :: nomres(5)
    real(kind=8)      :: sech, temp, alpha, beta, ad, t0_C, poro
    real(kind=8)      :: satu
!
    nomres(1) = 'VG_PR'
    nomres(2) = 'VG_N'
    nomres(3) = 'ATH'
    nomres(4) = 'TEMP_0_C'
    nomres(5) = 'PORO'
    call rcvala(imate, ' ', 'BETON_DESORP', 0, ' ', [0.d0], &
                5, nomres, valres, codret, 0)
    ASSERT(codret(1) .eq. 0)
    alpha = valres(1)
    beta = valres(2)
    ad = valres(3)
    t0_C = valres(4)
    poro = valres(5)
!
    call rcvarc(' ', 'SECH', '-', fami, kpg, ksp, sech, iret)
    if (iret .ne. 0) call utmess('F', 'COMPOR2_95', sk='SECH')
!
    call rcvarc(' ', 'TEMP', '-', fami, kpg, ksp, temp, iret)
    if (iret .ne. 0) call utmess('F', 'COMPOR2_95', sk='TEMP')

    satu = sech/poro/1.d3
    call leverettIsot(temp, satu, alpha, beta, ad, t0_C, hygr_prev)
!
    call rcvarc(' ', 'SECH', '+', fami, kpg, ksp, sech, iret)
    if (iret .ne. 0) call utmess('F', 'COMPOR2_95', sk='SECH')
!
    call rcvarc(' ', 'TEMP', '+', fami, kpg, ksp, temp, iret)
    if (iret .ne. 0) call utmess('F', 'COMPOR2_95', sk='TEMP')

    satu = sech/poro/1.d3
    call leverettIsot(temp, satu, alpha, beta, ad, t0_C, hygr_curr)

end subroutine leverettIsotMeca
