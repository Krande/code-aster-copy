! --------------------------------------------------------------------
! Copyright (C) 1991 - 2024 - EDF R&D - www.code-aster.org
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
subroutine leverettIsot(temp, satuIn, alpha, beta, ad, t0_C, hygr, dpc)
    implicit none
#include "asterfort/utmess.h"
    real(kind=8), intent(in) :: temp, satuIn, alpha, beta, ad, t0_C
    real(kind=8), intent(out) :: hygr
    real(kind=8), intent(out), optional :: dpc
!
    real(kind=8) :: gamma0, gamma, tempK, t0_K, dtemp, KT_K0, a, pc
    real(kind=8) :: satu, K0_KT
    real(kind=8), parameter :: epsil = 1d-6
!   waterMolarMass (g/mol)
    real(kind=8), parameter :: waterMolarMass = 18.01528d-3
!   idealGasConstant (J/K/mol)
    real(kind=8), parameter :: idealGasConstant = 8.314d0
!       to Kelvin
    tempK = temp+273.15d0
    t0_K = t0_C+273.15d0
    dtemp = tempK-t0_K
    satu = satuIn

    if (satu .lt. 0.d0 .or. satu .gt. 1.d0) call utmess('F', 'COMPOR2_96', sr=satu)
    if (temp .gt. 300.d0) call utmess('F', 'COMPOR2_97')
!
    if (satu .lt. epsil) then
        satu = epsil
        call utmess('A', 'COMPOR2_98')
    end if
    if (satu .gt. 1.d0-epsil) then
        satu = 1.d0-epsil
        call utmess('A', 'COMPOR2_98')
    end if
!
    gamma0 = surfaceTension(t0_K)
    gamma = surfaceTension(tempK)

    a = waterDensity(t0_K)*idealGasConstant*t0_K/(alpha*waterMolarMass)
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
                                        -((tempK-273.15d0)/374.14d0)**(1.d0/0.55d0))**0.55d0
!
    end function waterDensity
!   --------------------------------------------------------------------
    function hrKelvinLaw(pc, tempK)
        real(kind=8) :: pc, tempK
        real(kind=8) :: hrKelvinLaw
!
        hrKelvinLaw = exp(pc*waterMolarMass/(waterDensity(tempK)*idealGasConstant*tempK))
!
    end function hrKelvinLaw
!
end subroutine leverettIsot
