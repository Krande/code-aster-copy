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
subroutine leverettIsot(temp, satuIn, alpha, beta, ad, t0_C, hygr, dpc, pc_)
    implicit none
#include "rgi_module.h"
#include "asterfort/utmess.h"
#include "asterc/r8t0.h"
    real(kind=8), intent(in) :: temp, satuIn, alpha, beta, ad, t0_C
    real(kind=8), intent(out) :: hygr
    real(kind=8), intent(out), optional :: dpc, pc_
!
    real(kind=8) :: gamma0, gamma, tempK, t0_K, dtemp, KT_K0, a, pc
    real(kind=8) :: satu, K0_KT, tz0
!       to Kelvin
    tz0 = r8t0()
    tempK = temp+tz0
    t0_K = t0_C+tz0
    dtemp = tempK-t0_K
    satu = satuIn

    if (satu .lt. 0.d0 .or. satu .gt. 1.d0) call utmess('F', 'COMPOR2_96', sr=satu)
    if (temp .gt. 300.d0) call utmess('F', 'COMPOR2_97')
!
    if (satu .lt. EPSIL) then
        satu = EPSIL
        call utmess('A', 'COMPOR2_98')
    end if
    if (satu .gt. 1.d0-EPSIL) then
        satu = 1.d0-EPSIL
        call utmess('A', 'COMPOR2_98')
    end if
!
    gamma0 = surfaceTension(t0_K)
    gamma = surfaceTension(tempK)

    a = waterDensity(t0_K)*KGAZP*t0_K/(alpha*WATERMOLARMASS)
    KT_K0 = 10.d0**(ad*(2.d-3*dtemp-1d-6*dtemp**2))
    K0_KT = 1.d0/KT_K0
    pc = -a*((satu**(-1.d0/beta)-1.d0)**(1.d0-beta))*(gamma0/gamma)*sqrt(K0_KT)
    hygr = hrKelvinLaw(pc, tempK)

    if (present(dpc)) then
!       Isothermal desorption derivative of Van-guenuchten
        dpc = ((gamma/gamma0)*((KT_K0)**0.5)*a*(-beta+1.d0)/(beta)) &
              *(satu**(-(1.d0+beta)/beta)) &
              *(-1.d0+satu**(-1.d0/beta))**(-beta)
    end if

    if (present(pc_)) then
!       Isothermal desorption
        pc_ = pc
    end if

contains
!   --------------------------------------------------------------------
    function surfaceTension(tempK)
        real(kind=8) :: tempK
        real(kind=8) :: surfaceTension
!
!   Surface tension of water
!
        surfaceTension = 0.1558d0*(1.d0-(tempK/647.1d0))**1.26d0
!
    end function surfaceTension
!   --------------------------------------------------------------------
    function waterDensity(tempK)
        real(kind=8) :: tempK
        real(kind=8) :: waterDensity
!
!   Density of liquid water
!
        waterDensity = 314.4d0+685.6d0*(1.d0 &
                                        -((tempK-tz0)/374.14d0)**(1.d0/0.55d0))**0.55d0
!
    end function waterDensity
!   --------------------------------------------------------------------
    function hrKelvinLaw(pc, tempK)
        real(kind=8) :: pc, tempK
        real(kind=8) :: hrKelvinLaw
!
        hrKelvinLaw = exp(pc*WATERMOLARMASS/(waterDensity(tempK)*KGAZP*tempK))
!
    end function hrKelvinLaw
!
end subroutine leverettIsot
