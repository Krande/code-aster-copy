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
!
subroutine rftDiffusion(fami, kpg, ksp, poum, imate, c, &
                        temp, diff)
    implicit none
#include "asterfort/assert.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcvalb.h"
#include "asterfort/rcvarc.h"
#include "asterfort/utmess.h"
    character(len=*), intent(in) :: fami, poum
    integer, intent(in) :: kpg, ksp, imate
    real(kind=8), intent(in) :: c, temp
    real(kind=8), intent(out) :: diff
! ......................................................................
!   RFT law : diffusion coefficient calculation (SECH_RFT)
!
!
!   c (in) : concentration en eau
!   temp (in) : temperature
!
!   diff (out) : coefficient de diffusion
! ......................................................................
    integer           :: codret(5), nbpar
    real(kind=8)      :: valres(5), hygr, valpar(2), dpc
    real(kind=8)      :: richardsDiffusionCoef, vapourDiffusionCoef
    real(kind=8)      :: perm_in, qsr_k, poro, a_mil, b_mil, t0_C, vg_m_p
    real(kind=8)      :: t0_K, tempK, beta, satu
    character(len=16) :: nomres(5)
    character(len=8) :: nompar(2)
!   waterMolarMass (g/mol)
    real(kind=8), parameter :: waterMolarMass = 18.01528d-3
!   idealGasConstant (J/K/mol)
    real(kind=8), parameter :: idealGasConstant = 8.314d0
!
!   --------------------------------------------------------------------
!
!
    !   rft parameters
    nomres(1) = 'PERM_IN'
    nomres(2) = 'QSR_K'
    nomres(3) = 'A_MIL'
    nomres(4) = 'B_MIL'
    nomres(5) = 'VG_M_P'

    call rcvalb(fami, kpg, ksp, poum, imate, &
                ' ', 'SECH_RFT', 0, ' ', [0.d0], &
                5, nomres, valres, codret, 1)

    perm_in = valres(1)
    qsr_k = valres(2)
    a_mil = valres(3)
    b_mil = valres(4)
    vg_m_p = valres(5)

    nomres(1) = 'FONC_DESORP'
    nbpar = 2
    nompar(1) = 'TEMP'
    valpar(1) = c
    nompar(2) = 'TSEC'
    valpar(2) = temp
    call rcvalb(fami, kpg, ksp, poum, imate, &
                ' ', 'BETON_DESORP', nbpar, nompar, valpar, &
                1, nomres, valres, codret, 0)
!
    if (codret(1) .eq. 0) then
        call utmess('F', 'ALGORITH10_22')
    else
!       leverett isotherm
        call leverettIsotTher(c, temp, imate, hygr, dpc, poro, t0_C, beta)
    end if
!
    satu = c/poro/1.d3
    t0_K = t0_C+273.15d0
    tempK = temp+273.15d0

    call vapourDiffusion(satu, tempK, dpc, hygr, poro, a_mil, b_mil, &
                         vapourDiffusionCoef)

    call richardsDiffusion(satu, tempK, perm_in, t0_K, qsr_k, &
                           vg_m_p, poro, beta, dpc, richardsDiffusionCoef)

    diff = vapourDiffusionCoef+richardsDiffusionCoef
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
contains
!
!   --------------------------------------------------------------------
    subroutine leverettIsotTher(c, temp, imate, hygr, dpc, poro, t0_C, beta)
!
#include "asterfort/rcvala.h"
#include "asterfort/leverettIsot.h"
!
!      .................................................................
!      evaluation of hygrometry with leverett isotherm for thermic
!
!      c (in) : concentration en eau
!      temp (in) : temperature
!
!      hygr (out) : hygrometry
!      dpc (out) : capillar pressure derivative
!
!      .................................................................
        integer, intent(in) :: imate
        real(kind=8), intent(in) :: c, temp
        real(kind=8), intent(out) :: hygr, dpc, beta, poro, t0_C
!
        integer           :: codret(5)
        real(kind=8)      :: valres(5)
        character(len=16) :: nomres(5)
        real(kind=8)      :: alpha, ad, satu
!
        nomres(1) = 'VG_PR'
        nomres(2) = 'VG_N'
        nomres(3) = 'ATH'
        nomres(4) = 'TEMP_0_C'
        nomres(5) = 'PORO'
!
        call rcvala(imate, ' ', 'BETON_DESORP', 0, ' ', [0.d0], &
                    5, nomres, valres, codret, 0)
        ASSERT(codret(1) .eq. 0)
        alpha = valres(1)
        beta = valres(2)
        ad = valres(3)
        t0_C = valres(4)
        poro = valres(5)
        satu = c/poro/1.d3
!
        call leverettIsot(temp, satu, alpha, beta, ad, t0_C, hygr, dpc)

    end subroutine leverettIsotTher
!   --------------------------------------------------------------------
!   --------------------------------------------------------------------
    subroutine vapourDiffusion(satu, tempK, dpc, hygr, poro, a_mil, b_mil, &
                               vapDiff)
!
!      .................................................................
!      vapour
!      .................................................................
        real(kind=8), intent(in) :: satu, tempK, dpc, hygr
        real(kind=8), intent(in) :: poro, a_mil, b_mil
        real(kind=8), intent(out) :: vapDiff
!
        real(kind=8) :: fickDiff, Pv, pvs
        real(kind=8), parameter :: alpha_rankine = 13.7d0
        real(kind=8), parameter :: beta_rankine = 5120.d0
        real(kind=8), parameter :: Pa = 101325.d0
!
        fickDiff = fickDiffusion(satu, poro, a_mil, b_mil, tempk)
        ! Rankine equation for calculating saturated vapor pressure
        pvs = Pa*exp(alpha_rankine-(beta_rankine/tempK))
        !
        Pv = hygr*pvs
        vapDiff = (fickDiff*Pv*dpc*(waterMolarMass/(idealGasConstant*TempK))**2) &
                  /(poro*waterDensity(tempK)**2)

    end subroutine vapourDiffusion
!   --------------------------------------------------------------------
!   --------------------------------------------------------------------
    subroutine richardsDiffusion(satu, tempK, perm_in, t0_K, qsr_k, &
                                 vg_m_p, poro, beta, dpc, richDiff)
!
!      .................................................................
!      vapour
!      .................................................................
        real(kind=8), intent(in) :: satu, tempK, perm_in, t0_K, qsr_k
        real(kind=8), intent(in) :: vg_m_p, poro, beta, dpc
        real(kind=8), intent(out) :: richDiff
!
        real(kind=8) :: liquPerm, liquVisc, vgmRelativePerm
!
        ! liquid permeability
        liquPerm = perm_in*exp(exp((tempK-t0_K)/qsr_k)-1.0d0)
        ! liquid viscosity
        liquVisc = 0.6612*(tempk-229d0)**(-1.562d0)
        ! Mualem-Van-guenuchten relative permeability
        vgmRelativePerm = (satu**vg_m_p)*(1.d0-(1.d0-satu**(1.d0/beta))**beta)**2
        !
        richDiff = (liquPerm/(poro*liquVisc))*vgmRelativePerm*dpc

    end subroutine richardsDiffusion
!   --------------------------------------------------------------------
    function fickDiffusion(satu, poro, a_mil, b_mil, tempk)
        real(kind=8) :: satu, poro, a_mil, b_mil, tempk
        real(kind=8) :: fickDiffusion
!
        real(kind=8) :: D0, resitanceFactor
!
!   Effective diffusion coefficient of concrete
!
        ! Diffusion of water vapour in air
        D0 = 0.217d0*1d-4*((tempk/273.15))**1.88
        ! Millington resistance factor
        resitanceFactor = ((poro)**a_mil)*(1.d0-satu)**b_mil
        fickDiffusion = D0*resitanceFactor
!
    end function fickDiffusion
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
!
end subroutine rftDiffusion
