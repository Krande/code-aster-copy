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
! aslint: disable=W1504
! person_in_charge: sylvie.granet at edf.fr
!
subroutine thmCpl001(ds_thm, &
                     lMatr, lSigm, lVari, angl_naut, &
                     ndim, nbvari, &
                     dimdef, dimcon, &
                     adcome, adcote, adcp11, &
                     addeme, addete, addep1, &
                     temp, &
                     dtemp, dp1, &
                     deps, epsv, depsv, &
                     tbiot, &
                     phi, rho11, satur, &
                     congem, congep, &
                     vintm, vintp, dsde, &
                     retcom)
!
    use THM_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/appmas.h"
#include "asterfort/calor.h"
#include "asterfort/capaca.h"
#include "asterfort/dhdt.h"
#include "asterfort/dhwdp1.h"
#include "asterfort/dilata.h"
#include "asterfort/dileau.h"
#include "asterfort/dmdepv.h"
#include "asterfort/dmwdp1.h"
#include "asterfort/dmwdt.h"
#include "asterfort/dqdeps.h"
#include "asterfort/dqdp.h"
#include "asterfort/dqdt.h"
#include "asterfort/dspdp1.h"
#include "asterfort/enteau.h"
#include "asterfort/inithm.h"
#include "asterfort/sigmap.h"
#include "asterfort/unsmfi.h"
#include "asterfort/viemma.h"
#include "asterfort/viporo.h"
#include "asterfort/virhol.h"
#include "asterfort/THM_type.h"
!
    type(THM_DS), intent(in) :: ds_thm
    aster_logical, intent(in) :: lMatr, lSigm, lVari
    real(kind=8), intent(in) :: angl_naut(3)
    integer(kind=8), intent(in) :: ndim, nbvari
    integer(kind=8), intent(in) :: dimdef, dimcon
    integer(kind=8), intent(in) :: adcome, adcote, adcp11
    integer(kind=8), intent(in) :: addeme, addete, addep1
    real(kind=8), intent(in) :: temp
    real(kind=8), intent(in) :: dtemp, dp1
    real(kind=8), intent(in) :: epsv, depsv, deps(6), tbiot(6)
    real(kind=8), intent(out) :: phi, rho11, satur
    real(kind=8), intent(in) :: congem(dimcon)
    real(kind=8), intent(inout) :: congep(dimcon)
    real(kind=8), intent(in) :: vintm(nbvari)
    real(kind=8), intent(inout) :: vintp(nbvari)
    real(kind=8), intent(inout) :: dsde(dimcon, dimdef)
    integer(kind=8), intent(out) :: retcom
!
! --------------------------------------------------------------------------------------------------
!
! THM
!
! Compute generalized stress and matrix for coupled quantities - 'LIQU_SATU'
!
! --------------------------------------------------------------------------------------------------
!
! In  ds_thm           : datastructure for THM
! In  perman           : .true. for no-transient problem
! In  option           : option to compute
! In  angl_naut        : nautical angles
!                        (1) Alpha - clockwise around Z0
!                        (2) Beta  - counterclockwise around Y1
!                        (1) Gamma - clockwise around X
! In  ndim             : dimension of space (2 or 3)
! In  nbvari           : total number of internal state variables
! In  dimdef           : dimension of generalized strains vector
! In  dimcon           : dimension of generalized stresses vector
! In  adcome           : adress of mechanic components in generalized stresses vector
! In  adcote           : adress of thermic components in generalized stresses vector
! In  adcp11           : adress of first component and first phase in generalized stresses vector
! In  addeme           : adress of mechanic components in generalized strains vector
! In  addete           : adress of thermic components in generalized strains vector
! In  addep1           : adress of capillary pressure in generalized strains vector
! In  temp             : temperature at end of current time step
! In  dtemp            : increment of temperature
! In  dp1              : increment of capillary pressure (here, only liquid)
! In  deps             : increment of mechanical strains vector
! In  epsv             : current volumic strain
! In  depsv            : increment of volumic strain
! In  tbiot            : Biot tensor
! Out phi              : porosity
! Out rho11            : volumic mass for liquid
! Out satur            : saturation
! In  congem           : generalized stresses - At begin of current step
! IO  congep           : generalized stresses - At end of current step
! In  vintm            : internal state variables - At begin of current step
! IO  vintp            : internal state variables - At end of current step
! IO  dsde             : derivative matrix
! Out retcom           : return code for error
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: i
    real(kind=8), parameter :: rac2 = sqrt(2.d0)
    real(kind=8) :: phim, phi0
    real(kind=8) :: alp11, alp12
    real(kind=8) :: cp11, cp12, cp21, cp22
    real(kind=8) :: rho11m, rho110, rho0
    real(kind=8) :: rho12, rho21, rho22
    real(kind=8) :: em, alpliq, cliq, csigm
    real(kind=8) :: coeps
    real(kind=8) :: m11m
    real(kind=8) :: epsvm, cs, mdal(6), dalal, alpha0, alphfi, cbiot, unsks
    aster_logical :: l_emmag
    real(kind=8) :: saturm, dsatur
    real(kind=8) :: dpad, dp2, signe
    real(kind=8) :: dmdeps(6), sigmp(6), dsdp1(6)
    real(kind=8) :: dqeps(6)
    integer(kind=8)      :: advico, advihy, vihrho, vicphi
    real(kind=8) :: ep, surf, shut, sbjh, wbjh, dpi
!
! --------------------------------------------------------------------------------------------------
!
    rho11 = 0.d0
    rho12 = 0.d0
    rho21 = 0.d0
    rho22 = 0.d0
    cp12 = 0.d0
    cp21 = 0.d0
    cp22 = 0.d0
    alp11 = 0.d0
    alp12 = 0.d0
    signe = -1.d0
    retcom = 0
    phi = 0.d0
    satur = 0.d0
    dp2 = 0.d0
    dpad = 0.d0
    dpi = 0.d0
    ep = 0.d0
    surf = 0.d0
    shut = 0.d0
    sbjh = 0.d0
    wbjh = 0.d0
!
! - Get storage parameters for behaviours
!
    advico = ds_thm%ds_behaviour%advico
    advihy = ds_thm%ds_behaviour%advihy
    vihrho = ds_thm%ds_behaviour%vihrho
    vicphi = ds_thm%ds_behaviour%vicphi
!
! - Get initial parameters
!
    phi0 = ds_thm%ds_parainit%poro_init
!
! - Get material parameters
!
    rho0 = ds_thm%ds_material%solid%rho
    csigm = ds_thm%ds_material%solid%cp
    rho110 = ds_thm%ds_material%liquid%rho
    cliq = ds_thm%ds_material%liquid%unsurk
    alpliq = ds_thm%ds_material%liquid%alpha
    cp11 = ds_thm%ds_material%liquid%cp
!
! - Storage coefficient
!
    l_emmag = ds_thm%ds_material%hydr%l_emmag
    em = ds_thm%ds_material%hydr%emmag
!
! - Evaluation of initial saturation
!
    ASSERT(ds_thm%ds_behaviour%satur_type .eq. SATURATED)
    saturm = 1.d0
    satur = 1.d0
    dsatur = 0.d0
!
! - Evaluation of initial porosity
!
    phi = vintm(advico+vicphi)+phi0
    phim = vintm(advico+vicphi)+phi0
!
! - Evaluation of initial mass/volumic mass
!
    m11m = congem(adcp11)
    rho11 = vintm(advihy+vihrho)+rho110
    rho11m = vintm(advihy+vihrho)+rho110
!
! - Prepare initial parameters for coupling law
!
    call inithm(ds_thm, &
                angl_naut, tbiot, phi0, &
                epsv, depsv, &
                epsvm, cs, mdal, dalal, &
                alpha0, alphfi, cbiot, unsks)
!
! ==================================================================================================
!
! Internal state variables
!
! ==================================================================================================
!
! - Evaluation of porosity and save it in internal variables
!
    if (lVari) then
! ----- Compute standard porosity
        if (ds_thm%ds_elem%l_dof_meca) then
            if (ds_thm%ds_elem%l_jhms) then
                phi = vintp(advico+vicphi)
            else
                call viporo(ds_thm, nbvari, &
                            advico, vicphi, &
                            dtemp, dp1, dp2, &
                            deps, depsv, &
                            signe, satur, unsks, phi0, &
                            cs, tbiot, cbiot, &
                            alpha0, alphfi, &
                            vintm, vintp, &
                            phi, phim, retcom)
            end if
        end if
! ----- Compute porosity with storage coefficient
        if (l_emmag) then
            call viemma(nbvari, vintm, vintp, &
                        advico, vicphi, &
                        phi0, dp1, dp2, signe, satur, &
                        em, phi, phim)
        end if
    end if
    if (retcom .ne. 0) then
        goto 30
    end if
!
! - Evaluation of volumic mass and save it in internal variables
!
    if (lVari) then
! ----- Compute volumic mass for water
        if (ds_thm%ds_elem%l_dof_ther) then
            call virhol(nbvari, vintm, vintp, &
                        advihy, vihrho, &
                        dtemp, dp1, dp2, dpad, &
                        cliq, alpliq, signe, &
                        rho110, rho11, rho11m, &
                        retcom)
        else
            call virhol(nbvari, vintm, vintp, &
                        advihy, vihrho, &
                        dtemp, dp1, dp2, dpad, &
                        cliq, 0.d0, signe, &
                        rho110, rho11, rho11m, &
                        retcom)
        end if
    end if
    if (retcom .ne. 0) then
        goto 30
    end if
!
! - Update differential thermal expansion ratio
!
    if (ds_thm%ds_elem%l_dof_meca .and. .not. ds_thm%ds_elem%l_jhms) then
        call dilata(ds_thm, angl_naut, phi, tbiot, alphfi)
    end if
!
! - Update Biot modulus
!
    if (ds_thm%ds_elem%l_dof_meca .and. .not. ds_thm%ds_elem%l_jhms) then
        call unsmfi(ds_thm, phi, tbiot, cs)
    end if
!
! ==================================================================================================
!
! Generalized stresses
!
! ==================================================================================================
!
    if (ds_thm%ds_elem%l_dof_ther) then
! ----- Compute thermal expansion of liquid
        alp11 = dileau(satur, phi, alphfi, alpliq)
! ----- Compute specific heat capacity
        call capaca(rho0, rho11, rho12, rho21, rho22, &
                    satur, phi, &
                    csigm, cp11, cp12, cp21, cp22, &
                    dalal, temp, coeps, retcom)
        if (retcom .ne. 0) then
            goto 30
        end if
        if (lSigm) then
! --------- Update enthalpy of liquid
            congep(adcp11+ndim+1) = congep(adcp11+ndim+1)+ &
                                    enteau(dtemp, alpliq, temp, &
                                           rho11, dp2, dp1, dpad, &
                                           signe, cp11)
! --------- Update "reduced" heat Q'
            congep(adcote) = congep(adcote)+ &
                             calor(mdal, temp, dtemp, deps, &
                                   dp1, dp2, signe, &
                                   alp11, alp12, coeps, ndim)
        end if
    end if
!
! - Update mechanical stresses from pressures
!
    if (lSigm) then
        if (ds_thm%ds_elem%l_dof_meca .and. .not. ds_thm%ds_elem%l_jhms) then
            call sigmap(ds_thm, &
                        satur, signe, tbiot, dp2, dp1, dpi, &
                        sigmp)
            do i = 1, 3
                congep(adcome+6+i-1) = congep(adcome+6+i-1)+sigmp(i)
            end do
            do i = 4, 6
                congep(adcome+6+i-1) = congep(adcome+6+i-1)+sigmp(i)*rac2
            end do
        end if
    end if
!
! - Compute quantity of mass from change of volume, porosity and saturation
!
    if (lSigm) then
! ----- Update quantity of mass of liquid
        congep(adcp11) = appmas(m11m, &
                                phi, phim, &
                                satur, saturm, &
                                rho11, rho11m, &
                                epsv, epsvm)
    end if
!
! ==================================================================================================
!
! Tangent matrix
!
! ==================================================================================================
!
    if (lMatr) then
!
! ----- Mechanic
!
        if (ds_thm%ds_elem%l_dof_meca .and. .not. ds_thm%ds_elem%l_jhms) then
! --------- Derivative of _pressure part_ of stresses by capillary pressure
            call dspdp1(ds_thm, signe, tbiot, satur, dsdp1, phi0, ep, surf, sbjh, wbjh)
            do i = 1, 3
                dsde(adcome+6+i-1, addep1) = dsde(adcome+6+i-1, addep1)+ &
                                             dsdp1(i)
            end do
            do i = 4, 6
                dsde(adcome+6+i-1, addep1) = dsde(adcome+6+i-1, addep1)+ &
                                             dsdp1(i)*rac2
            end do
! --------- Derivative of quantity of mass by volumic mass - Mechanical part (strains)
            call dmdepv(rho11, satur, tbiot, dmdeps)
            do i = 1, 6
                dsde(adcp11, addeme+ndim-1+i) = dsde(adcp11, addeme+ndim-1+i)+ &
                                                dmdeps(i)
            end do
        end if
!
! ----- Thermic
!
        if (ds_thm%ds_elem%l_dof_ther) then
! --------- Derivative of enthalpy of liquid by capillary pressure
            dsde(adcp11+ndim+1, addep1) = dsde(adcp11+ndim+1, addep1)+ &
                                          dhwdp1(signe, alpliq, temp, rho11)
! --------- Derivative of enthalpy of liquid by temperature
            dsde(adcp11+ndim+1, addete) = dsde(adcp11+ndim+1, addete)+ &
                                          dhdt(cp11)
! --------- Derivative of quantity of mass of liquid by temperature
            dsde(adcp11, addete) = dsde(adcp11, addete)+ &
                                   dmwdt(rho11, phi, satur, cliq, 0.d0, alp11)
! --------- Derivative of "reduced" heat Q' by temperature
            dsde(adcote, addete) = dsde(adcote, addete)+ &
                                   dqdt(coeps)
! --------- Derivative of "reduced" heat Q' by capillary pressure
            dsde(adcote, addep1) = dsde(adcote, addep1)+ &
                                   dqdp(signe, alp11, temp)
! --------- Derivative of "reduced" heat Q' by mechanical strains
            if (ds_thm%ds_elem%l_dof_meca .and. .not. ds_thm%ds_elem%l_jhms) then
                call dqdeps(mdal, temp, dqeps)
                do i = 1, 6
                    dsde(adcote, addeme+ndim-1+i) = dsde(adcote, addeme+ndim-1+i)+ &
                                                    dqeps(i)
                end do
            end if
        end if
! ----- Derivative of quantity of mass of liquid by capillary pressure
        dsde(adcp11, addep1) = dsde(adcp11, addep1)+ &
                               dmwdp1(rho11, signe, satur, dsatur, phi, &
                                      cs, cliq, 1.d0, l_emmag, em)
    end if
!
30  continue
!
end subroutine
