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
! aslint: disable=W1306

module lcgtn_module

    use scalar_newton_module, only: &
        newton_state, &
        utnewt

    use tenseur_dime_module, only: &
        proten, &
        kron, &
        voigt, &
        identity

    use visc_norton_module, only: &
        VISCO, &
        ViscoInit => Init, &
        f_dka, &
        dkv_dka, &
        f_vsc, &
        dkv_vsc, &
        ddka_vsc, &
        solve_slope_dka

    implicit none
    private
    public:: CONSTITUTIVE_LAW, Init, InitGradVari, InitViscoPlasticity, Integrate

#include "asterf_types.h"
#include "asterc/r8nnem.h"
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utmess.h"

    ! Material characteristics (without viscosity)
    type MATERIAL
        real(kind=8) :: lambda, deuxmu, troismu, troisk, young, nu
        real(kind=8) :: r0, rh, r1, g1, r2, g2, rk, p0, gk
        real(kind=8) :: q1, q2, f0, fc, fr, hc, dv
        real(kind=8) :: fn, pn, sn, c0, kf, ki, epc, b0
        real(kind=8) :: sig0
        real(kind=8) :: dam_bkn
        real(kind=8) :: c = 0.d0
        real(kind=8) :: r = 0.d0
        type(VISCO)  :: norton
    end type MATERIAL

    ! Parameters related to the coupling between Norton viscosity and GTN
    type VISC_GTN
        aster_logical :: is_visc = ASTER_FALSE
        type(VISCO) :: norton
        real(kind=8):: dka = 0.d0
        real(kind=8):: ts = 0.d0
        real(kind=8):: vsc = 0.d0
    end type VISC_GTN

    ! Porosity and damage
    type PORO_DAM
        real(kind=8):: poro, poro_nucl, poro_grow
        real(kind=8):: dam, dam_coal
        real(kind=8):: epcum
        real(kind=8):: dam_rate
    end type PORO_DAM

    ! GTN class
    type CONSTITUTIVE_LAW
        integer(kind=8)       :: exception = 0
        aster_logical :: elas, rigi, vari, pred
        aster_logical :: grvi = ASTER_FALSE
        aster_logical :: loaded = ASTER_FALSE
        integer(kind=8)       :: ndimsi, itemax
        real(kind=8)  :: theta
        real(kind=8)  :: cvuser
        real(kind=8)  :: dt
        real(kind=8)  :: phi = 0.d0
        real(kind=8)  :: telh
        real(kind=8)  :: telq
        real(kind=8)  :: dam
        real(kind=8)  :: jac
        real(kind=8)  :: kam
        type(VISC_GTN):: vgtn
        type(MATERIAL):: mat
    end type CONSTITUTIVE_LAW

    ! Post-treatment
    type POST_TREATMENT
        real(kind=8) :: sieq_erx = 0.d0
        real(kind=8) :: sieq_ecr = 0.d0
        real(kind=8) :: sieq_vsc = 0.d0
        real(kind=8) :: sieq_nlc = 0.d0
        integer(kind=8)      :: arret = 0
    end type POST_TREATMENT

    integer(kind=8), parameter:: ELASTIC_STATE = 0
    integer(kind=8), parameter:: PLASTIC_STATE = 1
    integer(kind=8), parameter:: SINGULAR_STATE = 2
    integer(kind=8), parameter:: BROKEN_STATE = 3

    aster_logical:: DBG = ASTER_FALSE
    real(kind=8) :: dbg_tmp1, dbg_tmp2, dbg_tmp3, dbg_tmp4

contains

! ==================================================================================================
!  OBJECT CREATION AND INITIALISATION
! ==================================================================================================

    function Init(ndimsi, option, fami, kpg, ksp, imate, itemax, precvg, parm_theta, deltat) &
        result(self)

        implicit none

        integer(kind=8), intent(in)          :: kpg, ksp, imate, itemax, ndimsi
        real(kind=8), intent(in)    :: precvg, deltat, parm_theta
        character(len=16), intent(in):: option
        character(len=*), intent(in) :: fami
        type(CONSTITUTIVE_LAW)      :: self
! --------------------------------------------------------------------------------------------------
! ndimsi        symmetric tensor dimension (2*ndim)
! option        computation option
! fami          Gauss point set
! kpg           Gauss point number
! ksp           Layer number (for structure elements)
! imate         material pointer
! itemax        max number of iterations for the solver
! precvg        required accuracy (with respect to stress level))
! parm_theta    theta value for the porosity prediction (theta-predictor)
! deltat    time increment
! --------------------------------------------------------------------------------------------------
        integer(kind=8), parameter   :: nbel = 2, nbec = 10, nben = 16
! --------------------------------------------------------------------------------------------------
        integer(kind=8)             :: iokel(nbel), iokec(nbec), ioken(nben)
        real(kind=8)        :: valel(nbel), valec(nbec), valen(nben)
        real(kind=8)        :: r8nan
        character(len=16)   :: nomel(nbel), nomec(nbec), nomen(nben)
! --------------------------------------------------------------------------------------------------
        data nomel/'E', 'NU'/
        data nomec/'R0', 'RH', 'R1', 'GAMMA_1', 'R2', 'GAMMA_2', 'RK', 'P0', 'GAMMA_M', &
            'EPSP_LUDERS'/
        data nomen/'Q1', 'Q2', 'PORO_INIT', 'COAL_PORO', 'COAL_ACCE', 'PORO_RUPT', &
            'NUCL_GAUSS_PORO', 'NUCL_GAUSS_PLAS', 'NUCL_GAUSS_DEV', &
            'NUCL_CRAN_PORO', 'NUCL_CRAN_INIT', 'NUCL_CRAN_FIN', &
            'NUCL_EPSI_PENTE', 'NUCL_EPSI_INIT', 'ENDO_CRIT_VISC', 'ENDO_CRIT_RUPT'/
! --------------------------------------------------------------------------------------------------

        ! Variables non initialisees
        r8nan = r8nnem()
        self%telh = r8nan
        self%telq = r8nan
        self%dam = r8nan
        self%jac = r8nan
        self%kam = r8nan

        ! Parametres generaux
        self%ndimsi = ndimsi
        self%itemax = itemax
        self%cvuser = precvg
        self%theta = parm_theta
        self%dt = deltat

        ! Options de calcul
        self%elas = option .eq. 'RIGI_MECA_ELAS' .or. option .eq. 'FULL_MECA_ELAS'
        self%rigi = option .eq. 'RIGI_MECA_TANG' .or. option .eq. 'RIGI_MECA_ELAS' &
                    .or. option .eq. 'FULL_MECA' .or. option .eq. 'FULL_MECA_ELAS'
        self%vari = option .eq. 'FULL_MECA_ELAS' .or. option .eq. 'FULL_MECA' &
                    .or. option .eq. 'RAPH_MECA'
        self%pred = option .eq. 'RIGI_MECA_ELAS' .or. option .eq. 'RIGI_MECA_TANG'

        ! Elasticity material parameters
        call rcvalb(fami, kpg, ksp, '+', imate, ' ', 'ELAS', 0, ' ', [0.d0], &
                    nbel, nomel, valel, iokel, 2)
        self%mat%young = valel(1)
        self%mat%nu = valel(2)

        ! Hardening material parameters
        call rcvalb(fami, kpg, ksp, '+', imate, ' ', 'ECRO_NL', 0, ' ', [0.d0], &
                    nbec, nomec, valec, iokec, 0)
        if (iokec(1) .ne. 0) call utmess('F', 'COMPOR1_52')
        self%mat%r0 = valec(1)
        self%mat%rh = merge(valec(2), 0.d0, iokec(2) .eq. 0)
        self%mat%r1 = merge(valec(3), 0.d0, iokec(3) .eq. 0)
        self%mat%g1 = merge(valec(4), 0.d0, iokec(4) .eq. 0)
        self%mat%r2 = merge(valec(5), 0.d0, iokec(5) .eq. 0)
        self%mat%g2 = merge(valec(6), 0.d0, iokec(6) .eq. 0)
        self%mat%rk = merge(valec(7), 0.d0, iokec(7) .eq. 0)
        self%mat%p0 = merge(valec(8), 0.d0, iokec(8) .eq. 0)
        self%mat%gk = merge(valec(9), 1.d0, iokec(9) .eq. 0)

        ! Luders hardening is not implemented in the GTN model
        ASSERT(iokec(10) .ne. 0)

        !  Damage parameters
        call rcvalb(fami, kpg, ksp, '+', imate, ' ', 'GTN', 0, ' ', [0.d0], &
                    nben, nomen, valen, ioken, 0)
        self%mat%q1 = valen(1)
        self%mat%q2 = valen(2)
        self%mat%f0 = valen(3)
        self%mat%fc = merge(valen(4), r8nan, ioken(4) .eq. 0)
        self%mat%hc = merge(valen(5)-1, r8nan, ioken(5) .eq. 0)
        self%mat%fr = merge(valen(6), r8nan, ioken(6) .eq. 0)
        self%mat%fn = merge(valen(7), 0.00d0, ioken(7) .eq. 0)
        self%mat%pn = merge(valen(8), 0.10d0, ioken(8) .eq. 0)
        self%mat%sn = merge(valen(9), 0.05d0, ioken(8) .eq. 0)
        self%mat%c0 = merge(valen(10), 0.00d0, ioken(10) .eq. 0)
        self%mat%ki = merge(valen(11), 0.05d0, ioken(11) .eq. 0)
        self%mat%kf = merge(valen(12), 0.15d0, ioken(12) .eq. 0)
        self%mat%b0 = merge(valen(13), 0.00d0, ioken(13) .eq. 0)
        self%mat%epc = merge(valen(14), 0.80d0, ioken(14) .eq. 0)
        self%mat%dv = merge(valen(15), 1.d0, ioken(15) .eq. 0)
        self%mat%dam_bkn = merge(valen(16), 1.d0, ioken(16) .eq. 0)

        ASSERT(ioken(1) .eq. 0)
        ASSERT(ioken(4) .eq. 0 .eqv. (ioken(5) .eq. 0 .or. ioken(6) .eq. 0))

! --------------------------------------------------------------------------------------------------
        ! Tests and complements
! --------------------------------------------------------------------------------------------------

        ! Elasticity
        ASSERT(self%mat%young .gt. 0.d0)
        ASSERT(self%mat%nu .gt. -1.d0 .and. self%mat%nu .lt. 0.5d0)

        self%mat%lambda = self%mat%young*self%mat%nu/((1+self%mat%nu)*(1-2*self%mat%nu))
        self%mat%deuxmu = self%mat%young/(1+self%mat%nu)
        self%mat%troismu = 1.5d0*self%mat%deuxmu
        self%mat%troisk = self%mat%young/(1.d0-2.d0*self%mat%nu)

        ! Damage parameters
        ASSERT(self%theta .ge. 0.d0 .and. self%theta .le. 1.d0)
        ASSERT(self%mat%q1 .gt. 0)
        ASSERT(self%mat%q2 .gt. 0)
        ASSERT(self%mat%f0 .gt. 0)

        ! No coalescence
        if (ioken(4) .ne. 0) then
            self%mat%fc = 0.d0
            self%mat%hc = 0.d0
            self%mat%fr = 1.d0/self%mat%q1

            ! Coalescence slope is given
        else if (ioken(5) .eq. 0) then
            ASSERT(self%mat%fc .lt. 1.d0/self%mat%q1)
            ASSERT(self%mat%hc .ge. 0)
            ASSERT(1.d0/self%mat%q1 .le. self%mat%fc+(1+self%mat%hc)*(1-self%mat%fc))

            self%mat%fr = self%mat%fc+(1.d0/self%mat%q1-self%mat%fc)/(1+self%mat%hc)

            ! Fracture porosity is given
        else if (ioken(6) .eq. 0) then
            ASSERT(self%mat%fc .lt. self%mat%fr)
            ASSERT(self%mat%fr .le. 1)
            ASSERT(self%mat%fr .le. 1.d0/self%mat%q1)

            self%mat%hc = max(0.d0, (1.d0/self%mat%q1-self%mat%fc)/(self%mat%fr-self%mat%fc)-1)

        end if
        ASSERT(self%mat%f0 .le. self%mat%fr)

        ! Coupling between damage and viscosity
        ASSERT(self%mat%dv .ge. 0.d0)
        ASSERT(self%mat%dv .le. 1.d0)

        ! Strees reference as the inital yield threshold
        self%mat%sig0 = f_ecro(self, 0.d0)

        ! Damage threshold above which the point is considered broken (numerical motivations)
        self%mat%dam_bkn = min(self%mat%dam_bkn, 1.d0-sqrt(self%cvuser))
        ASSERT(self%mat%dam_bkn .ge. 0.d0)

    end function Init

! ==================================================================================================
!  COMPLEMENTARY INITIALISATION FOR GRAD_VARI
! ==================================================================================================

    subroutine InitGradVari(self, fami, kpg, ksp, imate, lag, apg)

        implicit none
        integer(kind=8), intent(in)          :: kpg, ksp, imate
        real(kind=8), intent(in)     :: lag, apg
        character(len=*), intent(in) :: fami
        type(CONSTITUTIVE_LAW), intent(inout)      :: self
! --------------------------------------------------------------------------------------------------
! fami      Gauss point set
! kpg       Gauss point number
! ksp       Layer number (for structure elements)
! imate     material pointer
! lag       Lagrangian value
! apg       nonlocal hardening variable
! --------------------------------------------------------------------------------------------------
        integer(kind=8), parameter:: nb = 2
! --------------------------------------------------------------------------------------------------
        integer(kind=8)             :: iok(nb)
        real(kind=8)        :: vale(nb)
        character(len=16)   :: nom(nb)
! --------------------------------------------------------------------------------------------------
        data nom/'C_GRAD_VARI', 'PENA_LAGR'/
! --------------------------------------------------------------------------------------------------

        self%grvi = ASTER_TRUE

        call rcvalb(fami, kpg, ksp, '+', imate, ' ', 'NON_LOCAL', 0, ' ', [0.d0], &
                    nb, nom, vale, iok, 2)
        self%mat%c = vale(1)
        self%mat%r = vale(2)
        self%phi = lag+self%mat%r*apg

    end subroutine InitGradVari

! ==================================================================================================
!  COMPLEMENTARY INITIALISATION FOR VISCOPLASTICITY
! ==================================================================================================

    subroutine InitViscoPlasticity(self, visc, fami, kpg, ksp, imate, deltat)

        implicit none
        aster_logical, intent(in)            :: visc
        integer(kind=8), intent(in)                  :: kpg, ksp, imate
        real(kind=8), intent(in)             :: deltat
        character(len=*), intent(in)         :: fami
        type(CONSTITUTIVE_LAW), intent(inout):: self
! --------------------------------------------------------------------------------------------------
! visc      True if viscosity is present
! fami      Gauss point set
! kpg       Gauss point number
! ksp       Layer number (for structure elements)
! imate     material pointer
! deltat    time increment (instap - instam)
! --------------------------------------------------------------------------------------------------

        self%mat%norton = ViscoInit(visc, fami, kpg, ksp, imate, deltat)

    end subroutine InitViscoPlasticity

! ==================================================================================================
!  INTEGRATION OF THE CONSTITUTIVE LAW (MAIN ROUTINE)
! ==================================================================================================

    subroutine Integrate(self, eps, vim, sig, vip, deps_sig, dphi_sig, deps_vi, dphi_vi)

        implicit none
        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8), intent(in)           :: eps(:), vim(:)
        real(kind=8), intent(out)          :: sig(:), vip(:), deps_sig(:, :)
        real(kind=8), intent(out), optional:: dphi_sig(:), deps_vi(:), dphi_vi
! --------------------------------------------------------------------------------------------------
! eps       strain at the end of the time step
! vim       internal variables at the beginning of the time step
! sig       stress at the end of the time step (resi) or the beginning of the time step (not resi)
! vip       internal variables at the end of the time step
! deps_sig  derivee dsig / deps
! dphi_sig  derivee dsig / dphi   (grad_vari)
! deps_vi   derivee dka  / deps  (grad_vari)
! dphi_vi   derivee dka  / dphi  (grad_vari)
! --------------------------------------------------------------------------------------------------
        integer(kind=8)         :: state
        real(kind=8)    :: kam, dka, epm(self%ndimsi), ep(self%ndimsi), rac2(self%ndimsi)
        real(kind=8)    :: vdum1(self%ndimsi), vdum2(self%ndimsi), rdum
        real(kind=8)    :: sigm(self%ndimsi)
        real(kind=8)    :: epcum_m, poro_m, poro_nucl_m, dam_m, dam_rate_m, dam_extr
        type(PORO_DAM):: pd_m, pd
        type(POST_TREATMENT):: post
! --------------------------------------------------------------------------------------------------

        if (DBG) write (6, *) 'START Integrate'
        if (DBG) write (6, *) 'eps = ', eps
        if (DBG) write (6, *) 'vim = ', vim

        ASSERT(present(dphi_sig) .eqv. self%grvi)
        ASSERT(present(deps_vi) .eqv. self%grvi)
        ASSERT(present(dphi_vi) .eqv. self%grvi)

        ! unpack internal variables
        rac2 = voigt(self%ndimsi)
        kam = vim(1)
        poro_m = vim(2)
        epm = vim(4:3+self%ndimsi)*rac2
        epcum_m = vim(10)
        poro_nucl_m = vim(11)
        dam_m = vim(12)
        dam_rate_m = vim(13)
        sigm = vim(14:13+self%ndimsi)*rac2

        ! Initial porosity and damage
        pd_m = Set_Porosity_Damage(self, poro_m, poro_nucl_m, dam_m, epcum_m, dam_rate_m)

        ! Extrapolated damage for plasticity computation
        dam_extr = Extrapolate_Damage(self, pd_m)

        ! Special behaviour in case of a broken point
        if (is_broken(self, pd_m%dam) .or. is_broken(self, dam_extr)) then
            if (self%grvi) then
                call Broken_Point(self, eps, sigm, state, dka, ep, &
                                  sig, deps_sig, dphi_sig, deps_vi, dphi_vi)
            else
                call Broken_Point(self, eps, sigm, state, dka, ep, &
                                  sig, deps_sig, vdum1, vdum2, rdum)
            end if

            ! (Visco) Plasticity behaviour integration
        else
            if (self%grvi) then
                call ComputePlasticity(self, eps, kam, epm, dam_extr, state, dka, ep, &
                                       sig, deps_sig, dphi_sig, deps_vi, dphi_vi)
            else
                call ComputePlasticity(self, eps, kam, epm, dam_extr, state, dka, ep, &
                                       sig, deps_sig, vdum1, vdum2, rdum)
            end if
            if (self%exception .ne. 0) goto 999
        end if

        ! Porosity update
        if (self%vari) then
            pd = Update_Porosity_Damage(self, state, epm, ep, kam+dka, pd_m)

        end if

! Post-treatment
        if (self%vari) then
            post = Compute_Post(self, pd_m, pd, dam_extr, sig, dka)
        end if

! pack internal variables
        if (self%vari) then
            vip(1) = kam+dka
            vip(2) = pd%poro
            vip(3) = state
            vip(8:9) = 0
            vip(4:3+self%ndimsi) = ep/rac2
            vip(10) = pd%epcum
            vip(11) = pd%poro_nucl
            vip(12) = pd%dam
            vip(13) = pd%dam_rate
            vip(18:19) = 0
            vip(14:13+self%ndimsi) = sig/rac2
            vip(20) = dam_extr
            vip(21) = post%sieq_erx
            vip(22) = post%sieq_ecr
            vip(23) = post%sieq_vsc
            vip(24) = post%sieq_nlc
            vip(25) = post%arret
        end if

999     continue
        if (DBG) write (6, *) 'END Integrate'
    end subroutine Integrate

! ==================================================================================================
!  PLASTICITY AND DAMAGE COMPUTATION
! ==================================================================================================

    subroutine ComputePlasticity(self, eps, kam, epm, dam, state, dka, ep, &
                                 t, deps_t, dphi_t, deps_ka, dphi_ka)

        implicit none

        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8), intent(in)             :: eps(:), kam, epm(:), dam
        integer(kind=8), intent(out)         :: state
        real(kind=8), intent(out)            :: dka, ep(:), t(:)
        real(kind=8), intent(out)            :: deps_t(:, :), dphi_t(:), deps_ka(:), dphi_ka
! --------------------------------------------------------------------------------------------------
! eps       final strain
! kam       initial hardening variable
! epm       initial plastic strain
! dam       current damage
! state     final regime state
! dka       increment of hardening variable
! ep        final plastic strain
! t         final stress
! deps_t    derivate dt / deps
! dphi_t    derivate dt / dphi   (grad_vari)
! deps_ka   derivate dka / deps  (grad_vari)
! dphi_ka   derivate dka / dphi  (grad_vari)
! --------------------------------------------------------------------------------------------------
        integer(kind=8)     :: iteint, iteext, i, typmin, typmax, flow_type
        real(kind=8):: kr(size(eps)), cvsigm, cv_fine, cveps, cv_g, jac
        real(kind=8):: tel(size(eps)), telh, telq, teld(size(eps))
        real(kind=8):: tels, deph, depd(size(eps))
        real(kind=8):: gamm1, desh, desq, pin, chim1, cs
        real(kind=8):: kv, kvmin, kvmax
        real(kind=8):: dkas, dkaz, dkamin, dkamax, dkaini, coef
        real(kind=8):: tsm, tss, tshmin, tshmax
        real(kind=8):: kml1, kml2, lbd, lbd1, lbd2
        real(kind=8):: equint, equext, d_equint, d_equext, p, q, ts
        real(kind=8):: p1, p2, fg1, fg2, sgn, dts_p, dkv_ts, dka_ts
        real(kind=8):: lambda_bar, deuxmu_bar
        real(kind=8):: mat(2, 2), vec(2)
        real(kind=8):: djac_p, dteh_p, dteq_p, dphi_p, deps_p(size(eps))
        real(kind=8):: deps_jac(size(eps)), deps_teh(size(eps)), deps_teq(size(eps))
        real(kind=8):: dp_q, dts_q, djac_q, dteh_q, dteq_q, dphi_q, deps_q(size(eps))
        real(kind=8):: djac_ts, dteh_ts, dteq_ts, dphi_ts
        real(kind=8):: djac_ka, dteh_ka, dteq_ka
        type(newton_state):: memint, memext
! ----------------------------------------------------------------------

        if (DBG) write (6, *) 'START Compute_Plasticity'

! --------------------------------------------------------------------------------------------------
!   Initialisation
! --------------------------------------------------------------------------------------------------

        kr = kron(self%ndimsi)
        ep = epm
        dka = 0.d0

        cvsigm = self%mat%sig0*self%cvuser
        cv_fine = 1.d-1*cvsigm
        cveps = cv_fine/self%mat%young

        ! Adjust viscosity with respect to damage
        call Set_Visc_Damage(self, dam)

        ! Caracteristiques des deformations
        jac = exp(sum(eps(1:3)))

        ! Contrainte elastique
        tel = self%mat%lambda*sum(eps(1:3)-epm(1:3))*kr+self%mat%deuxmu*(eps-epm)
        telh = sum(tel(1:3))/3.d0
        teld = tel-telh*kr
        telq = sqrt(1.5d0*dot_product(teld, teld))

        ! Perturbation en cvsigm*1e-3 pour eviter la singularite TELH=0
        telh = sign(max(cvsigm*1.d-3, abs(telh)), telh)

        ! Copie des informations dans self pour utilisation des fonctions du module
        self%kam = kam
        self%telh = telh
        self%telq = telq
        self%dam = dam
        self%jac = jac
        self%loaded = ASTER_TRUE

        ! Algorithmic parameters for transition between weak and strong viscosity
        call Set_Visc_Transition(self)

        ! Precision sur la fonction G
        cv_g = (1-dam)*self%cvuser

        ! Frontieres des regimes elastique et singulier
        gamm1 = (1-dam)**2/(2*dam)
        desh = telh/self%mat%troisk
        desq = telq/self%mat%deuxmu
        if (desh**2+desq**2 .eq. 0) then
            pin = 0
        else
            chim1 = (self%mat%q2*desq)**2/(9*desh**2/dam+(self%mat%q2*desq)**2)
            cs = 2*gamm1*(1-chim1)/(1+sqrt(1+2*gamm1*chim1*(1-chim1)))
            pin = 2/self%mat%q2*acosh(1+cs)*abs(desh)+2.d0/3.d0*(1-dam)*sqrt(1-cs/gamm1)*desq
        end if
        dkas = jac*pin
        tsm = f_ts_hat(self, dka=0.d0)
        tss = f_ts_hat(self, dka=dkas)

        if (DBG) write (6, *) 'telh  = ', telh
        if (DBG) write (6, *) 'telq  = ', telq
        if (DBG) write (6, *) 'jac   = ', jac
        if (DBG) write (6, *) 'dkas  = ', dkas
        if (DBG) write (6, *) 'tsm   = ', tsm
        if (DBG) write (6, *) 'tss   = ', tss

        if (DBG) write (6, *) 'INTEGRATION'
        if (DBG) write (6, *) 'rigi = ', self%rigi
        if (DBG) write (6, *) 'vari = ', self%vari

! --------------------------------------------------------------------------------------------------
!  INTEGRATION DE LA LOI DE COMPORTEMENT
! --------------------------------------------------------------------------------------------------

! ----------------------------------------------------------------------
!  REGIME SINGULIER SANS PLASTICITE
! ----------------------------------------------------------------------

        ! Elasticite (avec marge) si tels < T(kam) + cvsigm
        if (abs(tsm) .le. cvsigm) then
            state = SINGULAR_STATE
            ep = eps
            t = 0
            goto 500
        end if

! ----------------------------------------------------------------------
!  REGIME ELASTIQUE
! ----------------------------------------------------------------------

        ! Elasticite (avec marge) si tels < T(kam) + cvsigm
        if (is_tels_less_than(self, tsm+cvsigm)) then
            state = ELASTIC_STATE
            t = tel
            goto 500
        end if

! ----------------------------------------------------------------------
!  REGIME PLASTIQUE SINGULIER
! ----------------------------------------------------------------------

        ! Cas singulier standard (tss < 0) ou si T* petit (tels<cvsigm)
        if (tss .le. cvsigm .or. is_tels_less_than(self, cvsigm)) then

            ! Calcul de dkaz tq T_hat(dkaz) = 0
            if (abs(tss) .le. cvsigm) then
                dkaz = dkas
            else
                dkamin = merge(dkas, 0.d0, tss .lt. 0)
                call Solve_ts_hat(self, 0.d0, cvsigm, dkamin, dkaz)
                if (self%exception .eq. 1) goto 999
            end if

            ! Mise a jour de la plasticite et des contraintes
            dka = dkaz
            ep = eps
            t = 0
            state = SINGULAR_STATE
            goto 500
        end if

! ----------------------------------------------------------------------
!  REGIME PLASTIQUE REGULIER
! ----------------------------------------------------------------------

        ! 1. Calcul d'une borne min dkamin = max(0.d0 , dkaz) avec tolerance

        if (tsm .ge. cv_fine) then
            dkamin = 0.d0
            typmin = 1
            tshmin = tsm
        else
            call Solve_ts_hat(self, 2*cv_fine, cv_fine, 0.d0, dkaz)
            if (self%exception .eq. 1) goto 999
            dkamin = dkaz
            typmin = 2
            tshmin = f_ts_hat(self, dka=dkamin)
        end if

        ! Test si la borne min conduit a une solution facile
        ! precision numerique pour garantir dkamin < lbd
        if (typmin .eq. 2) then
            p = Solve_g_hat(self, tshmin, cv_g, 0.d0)
            if (self%exception .eq. 1) goto 999
            lbd = f_lbd_hat(self, p, tshmin)
            if (lbd .le. dkamin) then
                p = 0.d0
                q = 0.d0
                ts = 0.d0
                dka = dkamin
                state = SINGULAR_STATE
                goto 400
            end if
        end if

        ! 2. Calcul d'une borne max dkamax = min(dkas, dkael) ou T_hat(dkael)=tels

        tels = Compute_Star(self, 1.d0, 1.d0, cv_fine)
        if (self%exception .eq. 1) goto 999

        if (.not. is_tels_less_than(self, tss+cv_fine)) then
            dkamax = dkas
            typmax = 1
            tshmax = tss
        else
            call Solve_ts_hat(self, tels-1.5d0*cv_fine, 0.5d0*cv_fine, 0.d0, dkamax)
            if (self%exception .eq. 1) goto 999
            typmax = 2
            tshmax = f_ts_hat(self, dka=dkamax)
        end if
        if (DBG) write (6, *) 'typmax, dkamax, tshmax', typmax, dkamax, tshmax

        ASSERT(dkamin .le. dkamax)
        ASSERT(tshmin .le. tshmax)
        ASSERT(tshmin .le. tels)
        ASSERT(tshmin .gt. 0)
        ASSERT(.not. is_tels_less_than(self, tshmax))

        ! Test si la borne max conduit a une solution facile
        ! precision numerique pour garantir lbd < dkamax
        if (typmax .eq. 2) then
            p = Solve_g_hat(self, tshmax, cv_g, 0.d0)
            lbd = f_lbd_hat(self, p, tshmax)
            if (lbd .ge. dkamax) then
                p = 1.d0
                q = 1.d0
                ts = tels
                dka = dkamax
                state = ELASTIC_STATE
                goto 400
            end if
        end if

        ! 3. Estimation initiale par resolution du probleme de von Mises

        call Solve_Mises(self, cv_fine, dkaz, flow_type, dkaini)
        if (self%exception .eq. 1) goto 999

        ! Si Mises ne fournit pas une initialisation pertinente
        !   -> initialisation par une methode de corde
        if (flow_type .eq. 0 .or. flow_type .eq. 2 .or. &
            (flow_type .eq. 1 .and. (dkaini .le. dkamin .or. dkaini .ge. dkamax))) then
            if (DBG) write (6, *) 'tels     =', tels
            coef = (tels-tshmin)/(tels-tshmin+tshmax)
            dkaini = dkamin+coef*(dkamax-dkamin)
        end if

        ASSERT(dkaini .ge. dkamin .and. dkaini .le. dkamax)
        if (DBG) write (6, *) 'flow_type=', flow_type
        if (DBG) write (6, *) 'kam      = ', kam
        if (DBG) write (6, *) 'dkaini   =', dkaini
        if (DBG) write (6, *) 'dkamin   =', dkamin, '  tshmin = ', tshmin
        if (DBG) write (6, *) 'dkamax   =', dkamax, '  tshmax = ', tshmax

        ! 4. Preparation du changement de variable en presence de viscoplasticite

        ! La solution est-elle dans le domaine ou la viscoplasticite est dominante
        if (self%vgtn%is_visc) then
            if (self%vgtn%dka .le. dkamin) then
                self%vgtn%norton%active = ASTER_FALSE
            else if (self%vgtn%dka .ge. dkamax) then
                self%vgtn%norton%active = ASTER_TRUE
            else
                p = Solve_g_hat(self, self%vgtn%ts, cv_g, 0.d0)
                self%vgtn%norton%active = f_lbd_hat(self, p, self%vgtn%ts) .lt. self%vgtn%dka
            end if
        else
            self%vgtn%norton%active = ASTER_FALSE
        end if

        ! Correction des bornes si viscoplasticite dominante ou non
        if (self%vgtn%is_visc) then
            if (self%vgtn%norton%active) then
                kvmin = f_vsc(self%vgtn%norton, dka=dkamin)
                kvmax = min(self%vgtn%vsc, f_vsc(self%vgtn%norton, dka=dkamax))
                kv = min(self%vgtn%vsc, f_vsc(self%vgtn%norton, dka=dkaini))
            else
                kvmin = max(self%vgtn%dka, dkamin)
                kvmax = dkamax
                kv = max(self%vgtn%dka, dkaini)
            end if
        else
            kvmin = dkamin
            kvmax = dkamax
            kv = dkaini
        end if
        ASSERT(kvmin .le. kv .and. kv .le. kvmax)

        ! 5. Resolution du systeme non lineaire de deux equations
        ! Resolution par rapport a kv de : kml := dka - lbd(p,Ts) == 0

        do iteext = 1, self%itemax

            ! 3.1 Calcul de Ts(ka)
            if (DBG) write (6, *) '3.1'
            dka = f_dka(self%vgtn%norton, kv)
            ts = f_ts_hat(self, kv=kv)
            if (DBG) write (6, *) 'kv = ', kv, '  dka = ', dka, '  ts = ', ts

            ! 3.2 Resolution de G_hat(p,ts) == -prec/2 par rapport a p
            !     --> borne min de l'encadrement de la solution exacte
            !     (Rappel: G croissant par rapport a p)
            p1 = Solve_g_hat(self, ts, 0.5d0*cv_g, -0.5d0*cv_g)
            if (self%exception .eq. 1) goto 999
            fg1 = f_g_hat(self, p1, ts)
            lbd1 = f_lbd_hat(self, p1, ts)
            kml1 = dka-lbd1
            if (DBG) write (6, *) 'p1=', p1, '  lbd1=', lbd1, ' kml1=', kml1, '  fg1=', fg1

            ! Test si convergence de la boucle exterieure
            if (abs(ts-f_ts_hat(self, dka=lbd1)) .le. cvsigm) then
                if (DBG) write (6, *) 'CVG-1.SIGM'
                p = p1
                dka = lbd1
                exit
            end if
            if (abs(dka-lbd1) .le. cveps) then
                if (DBG) write (6, *) 'CVG-1.EPSI'
                p = p1
                exit
            end if

            ! 3.3 Resolution de G_hat(p,ts)=prec/2
            !     --> borne max de l'encadrement de la solution exacte
            if (DBG) write (6, *) '3.3'

            p2 = Solve_g_hat(self, ts, 0.5d0*cv_g, 0.5d0*cv_g, p_ini=p1)
            if (self%exception .eq. 1) goto 999
            fg2 = f_g_hat(self, p2, ts)

            lbd2 = f_lbd_hat(self, p2, ts)
            kml2 = dka-lbd2
            if (DBG) then
                dbg_tmp1 = f_lbd_hat(self, p2, ts)
                dbg_tmp2 = f_g_hat(self, p2, ts)
                write (6, *) 'p2=', p2, '  lbd2=', lbd2, '=', dbg_tmp1, ' kml2=', kml2, &
                    '  fg2=', fg2, '=', dbg_tmp2
            end if

            ! Test de convergence de la boucle exterieure sur la seconde borne
            if (abs(ts-f_ts_hat(self, dka=lbd2)) .le. cvsigm) then
                if (DBG) write (6, *) 'CVG-2.SIGM'
                p = p2
                dka = lbd2
                exit
            end if
            if (abs(dka-lbd2) .le. cveps) then
                if (DBG) write (6, *) 'CVG-2.EPSI'
                p = p2
                exit
            end if

            ! Valeur de p optimale -> celle qui conduit au meilleur lambda
            ASSERT(fg1*fg2 .le. 0)
            ASSERT(p2 .ge. p1)
            p = merge(p1, p2, abs(dka-lbd1) .le. abs(dka-lbd2))

            ! 3.4 Si l'intervalle (p1,p2) permet de trouver une solution KML(p,ts,ka)==0
            !   -> recherche de p
            if (DBG) write (6, *) '3.4'
            if (kml1*kml2 .lt. 0) then
                if (DBG) then
                    write (6, *) '3.4 -> GO'
                    dbg_tmp1 = f_g_hat(self, p1, ts)
                    dbg_tmp2 = f_g_hat(self, p2, ts)
                    dbg_tmp3 = f_lbd_hat(self, p1, ts)
                    dbg_tmp4 = f_lbd_hat(self, p2, ts)
                    write (6, *) 'dka = ', dka, '  ts = ', ts
                    write (6, *) 'lbd1 = ', lbd1, ' = ', dbg_tmp3
                    write (6, *) 'lbd2 = ', lbd2, ' = ', dbg_tmp4
                    write (6, *) 'p1=', p1, ' kml1=', kml1, '=', dka-lbd1, '  fg1=', &
                        fg1, '=', dbg_tmp1
                    write (6, *) 'p2=', p2, ' kml2=', kml2, '=', dka-lbd2, '  fg2=', &
                        fg2, '=', dbg_tmp2
                end if

                ! croissance ou decroissance de KML(p,ts,ka) par rapport a p dans [p1,p2]
                sgn = sign(1.d0, kml2-kml1)

                ! resolution dans l'intervalle
                if (DBG) write (6, *) 'utnewt 3B'
                do iteint = 1, self%itemax
                    lbd = f_lbd_hat(self, p, ts)
                    equint = sgn*(dka-lbd)
                    d_equint = -sgn*dp_lbd_hat(self, p, ts)
                    if (DBG) write (6, *) 'ite = ', iteint, '  p = ', p, &
                        '  lbd = ', lbd, '  equ = ', equint, '  d_equ = ', d_equint
                    if (DBG) then
                        dbg_tmp1 = f_ts_hat(self, dka=lbd)
                        write (6, *) 'cvg = ', ts-dbg_tmp1
                    end if
                    if (abs(ts-f_ts_hat(self, dka=lbd)) .le. cvsigm) then
                        dka = lbd
                        exit
                    end if
                    if (abs(dka-lbd) .le. cveps) exit
                    p = utnewt(p, equint, d_equint, iteint, memint, xmin=p1, xmax=p2)
                end do
                if (DBG) write (6, *) 'utnewt 3E'
                if (iteint .gt. self%itemax) then
                    self%exception = 1
                    goto 999
                end if
                exit
            end if

            ! 3.5 Sinon KML(p,ts,ka) <> 0 dans l'intervalle (p1,p2) -> nouvel itere pour dka
            equext = dka-f_lbd_hat(self, p, ts)

            dts_p = -dts_g_hat(self, p, ts)/dp_g_hat(self, p, ts)
            dkv_ts = dkv_ts_hat(self, kv)
            d_equext = dkv_dka(self%vgtn%norton, kv)-(dp_lbd_hat(self, p, ts)*dts_p &
                                                      +dts_lbd_hat(self, p, ts))*dkv_ts
            if (DBG) write (6, *) 'iteext=', iteext, '  kv=', kv, '  equext=', equext, &
                '  d_equext=', d_equext
            if (DBG) then
                dbg_tmp1 = f_lbd_hat(self, p, ts)
                write (6, *) 'dka=', dka, '  lambda=', dbg_tmp1, '  p=', p, '  ts=', ts
            end if
            if (DBG) write (6, *) 'utnewt 4B'
            kv = utnewt(kv, equext, d_equext, iteext, memext, xmin=kvmin, xmax=kvmax)
            if (DBG) write (6, *) 'utnewt 4E'

        end do
        if (iteext .gt. self%itemax) then
            self%exception = 1
            goto 999
        end if
        if (DBG) write (6, *) 'Calcul de q'
        q = f_q_hat(self, p, ts)
        state = PLASTIC_STATE

400     continue
        if (DBG) write (6, *) '6.1'
        deph = (1-p)*telh/self%mat%troisk
        depd = (1-q)*teld/self%mat%deuxmu
        ep = ep+deph*kr+depd
        t = p*telh*kr+q*teld

500     continue

        if (DBG) write (6, *) 'dka=', dka
        if (DBG) write (6, *) 'state=', state
        if (DBG) write (6, *) 'tau=', t
        if (DBG) write (6, *) 'deph=', deph
        if (DBG) write (6, *) 'depd=', depd
        if (DBG) write (6, *) 'ep=', ep
        if (DBG) write (6, *) 'tsm=', tsm

! ======================================================================
!                           MATRICES TANGENTES
! ======================================================================

        if (.not. self%rigi) goto 999

        if (DBG) write (6, *) 'MATRICE TANGENTE'
        deps_t = 0
        dphi_t = 0
        deps_ka = 0
        dphi_ka = 0

        ! Correction eventuelle pour choisir l'operateur tangent en phase de prediction
        ! Si elastique mais presque plastique -> plastique

        if (self%pred .and. state .eq. ELASTIC_STATE .and. .not. self%elas) then
            if (.not. is_tels_less_than(self, tsm-cvsigm)) then
                if (DBG) write (6, *) 'Matrice tangente plastique forcee'
                state = PLASTIC_STATE
                p = 1.d0
                q = 1.d0
                ts = f_ts_hat(self, dka=dka)
                t = tel
            end if
        end if

        ! Regime elastique (seul deps_t est non nulle, egale a la matrice d'elasticite)
        if (state .eq. ELASTIC_STATE .or. self%elas) then
            deps_t(1:3, 1:3) = self%mat%lambda
            do i = 1, size(eps)
                deps_t(i, i) = deps_t(i, i)+self%mat%deuxmu
            end do

            ! Regime plastique regulier
        else if (state .eq. PLASTIC_STATE) then

            ! Variations des inconnues principales (p,ka)
            mat(1, 1) = dts_m_hat(self, p, ts, dka)
            mat(1, 2) = -dts_g_hat(self, p, ts)
            mat(2, 1) = -dp_m_hat(self, p, ts, dka)
            mat(2, 2) = dp_g_hat(self, p, ts)
            mat = mat/(mat(1, 1)*mat(2, 2)-mat(1, 2)*mat(2, 1))

            vec = -matmul(mat, [dteh_g_hat(self, p, ts), dteh_m_hat(self, p, ts, dka)])
            dteh_p = vec(1)
            dteh_ts = vec(2)

            vec = -matmul(mat, [dteq_g_hat(self, p, ts), dteq_m_hat(self, p, ts, dka)])
            dteq_p = vec(1)
            dteq_ts = vec(2)

            vec = -matmul(mat, [0.d0, djac_m_hat(self, p, ts, dka)])
            djac_p = vec(1)
            djac_ts = vec(2)

            vec = -matmul(mat, [0.d0, dphi_m_hat(self, p, ts, dka)])
            dphi_p = vec(1)
            dphi_ts = vec(2)

            ! Variations de q
            dp_q = dp_q_hat(self, p, ts)
            dts_q = dts_q_hat(self, p, ts)
            dteh_q = dp_q*dteh_p+dts_q*dteh_ts+dteh_q_hat(self, p, ts)
            dteq_q = dp_q*dteq_p+dts_q*dteq_ts
            djac_q = dp_q*djac_p+dts_q*djac_ts
            dphi_q = dp_q*dphi_p+dts_q*dphi_ts

            ! Variations de ka
            dka_ts = dka_ts_hat(self, dka)
            dteh_ka = (dteh_ts)/dka_ts
            dteq_ka = (dteq_ts)/dka_ts
            djac_ka = (djac_ts-djac_ts_hat(self, dka))/dka_ts
            dphi_ka = (dphi_ts-dphi_ts_hat(self, dka))/dka_ts

            ! Variations des invariants par rapport a epsilon
            deps_teh = self%mat%troisk/3.d0*kr
            deps_teq = self%mat%troismu*teld/telq
            deps_jac = jac*kr

            ! dt/deps
            lambda_bar = (p*self%mat%troisk-q*self%mat%deuxmu)/3.d0
            deuxmu_bar = q*self%mat%deuxmu
            deps_t(1:3, 1:3) = lambda_bar
            do i = 1, size(eps)
                deps_t(i, i) = deps_t(i, i)+deuxmu_bar
            end do
            deps_p = djac_p*deps_jac+dteh_p*deps_teh+dteq_p*deps_teq
            deps_q = djac_q*deps_jac+dteh_q*deps_teh+dteq_q*deps_teq
            deps_t = deps_t+proten(telh*kr, deps_p)+proten(teld, deps_q)

            ! dka/deps
            deps_ka = djac_ka*deps_jac+dteh_ka*deps_teh+dteq_ka*deps_teq

            ! dt/dphi
            dphi_t = dphi_p*telh*kr+dphi_q*teld

            ! dka/dphi -> deja calcule ci-dessus

            ! Regime singulier (deps_t et dphi_t sont nulles)
        else if (state .eq. SINGULAR_STATE) then
            deps_jac = jac*kr
            dka_ts = dka_ts_hat(self, dka)
            dphi_ka = -dphi_ts_hat(self, dka)/dka_ts
            djac_ka = -djac_ts_hat(self, dka)/dka_ts
            deps_ka = djac_ka*deps_jac

        else
            ASSERT(.false.)
        end if

!    Fin de la routine
999     continue

        if (DBG) write (6, *) 'Code retour = ', self%exception
        ASSERT(.not. DBG .or. self%exception .eq. 0)
        if (DBG) write (6, *) 'END Compute_Plasticity'
    end subroutine ComputePlasticity

! ==================================================================================================
!  Is Tel_star less than a value ?
! ==================================================================================================

    aster_logical function is_tels_less_than(self, thre)

        implicit none

        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8), intent(in)::thre
! --------------------------------------------------------------------------------------------------
! thre      threshold: return True if tels < thre
! --------------------------------------------------------------------------------------------------

        if (DBG) write (6, *) 'BEGIN is_tels_less_than'
        is_tels_less_than = ASTER_FALSE
        ASSERT(self%loaded)

        ! Precaution pour le calcul numerique de G
        if (thre .le. 0) goto 100
        if (self%telq .gt. thre*sqrt(1+self%dam**2)) goto 100
        if (1.5d0*self%mat%q2*self%telh .gt. thre*acosh((1+self%dam**2)/(2*self%dam))) goto 100

        !   On exploite le caractere decroissant de G et le fait que G(Tel, Tels) = 0
        is_tels_less_than = to_aster_logical(f_g(self, 1.d0, 1.d0, thre) .lt. 0)

100     continue
        if (DBG) write (6, *) 'END is_tels_less_than'
    end function is_tels_less_than

! ==================================================================================================
!  Solve f_ts_hat(dka) = ts with dka > 0
! ==================================================================================================

    subroutine Solve_ts_hat(self, ts_target, cv_ts, dkamin, dka)

        implicit none

        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8), intent(in)  ::ts_target, cv_ts, dkamin
        real(kind=8), intent(out) ::dka
! --------------------------------------------------------------------------------------------------
! ts_target     target equivalent stress
! cv_ts         expected accuracy on the residual
! exi_ka        True if a positive solution exist, else False
! dka           solution hardening variable
! --------------------------------------------------------------------------------------------------
        aster_logical:: buffer
        integer(kind=8)     :: ite
        real(kind=8):: ts_min, equ, d_equ, kv, kvmin
        type(newton_state):: mem
! --------------------------------------------------------------------------------------------------
        if (DBG) write (6, *) 'BEGIN Solve_ts_hat'
        buffer = self%vgtn%norton%active

        ! Lower bound
        ts_min = f_ts_hat(self, dka=dkamin)
        if (abs(ts_min-ts_target) .le. cv_ts) then
            dka = dkamin
            goto 999
        end if

        ! Existence
        ASSERT(ts_min-ts_target .le. 0)

        ! Zone de viscosite dominante (chgt de variable) ou non
        if (.not. self%vgtn%is_visc) then
            self%vgtn%norton%active = ASTER_FALSE
            kv = dkamin
        else if (self%vgtn%ts .le. ts_target) then
            self%vgtn%norton%active = ASTER_FALSE
            kv = max(dkamin, self%vgtn%dka)
        else
            self%vgtn%norton%active = ASTER_TRUE
            kv = max(0.d0, f_vsc(self%vgtn%norton, dka=dkamin))
        end if
        kvmin = kv

        ! Newton method
        do ite = 1, self%itemax
            equ = f_ts_hat(self, kv=kv)-ts_target
            d_equ = dkv_ts_hat(self, kv)
            if (abs(equ) .le. cv_ts) exit
            if (DBG) write (6, *) 'utnewt 5B'
            kv = utnewt(kv, equ, d_equ, ite, mem, xmin=kvmin)
            if (DBG) write (6, *) 'utnewt 5E'
        end do
        if (ite .gt. self%itemax) then
            self%exception = 1
            goto 999
        end if
        dka = f_dka(self%vgtn%norton, kv)

999     continue
        self%vgtn%norton%active = buffer
        if (DBG) write (6, *) 'END Solve_ts_hat'
    end subroutine Solve_ts_hat

! ==================================================================================================
!  Solve f_g_hat(p,ts) = g0 with respect to p
! ==================================================================================================

    function Solve_g_hat(self, ts, cvg, g0, p_ini) result(p)

        implicit none

        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8), intent(in)  ::ts, g0, cvg
        real(kind=8), intent(in), optional::p_ini
        real(kind=8) ::p
! --------------------------------------------------------------------------------------------------
! ts     equivalent stress
! g0     right hand side term
! cvg    accuracy
! --------------------------------------------------------------------------------------------------
        integer(kind=8)     :: ite
        real(kind=8):: equ, d_equ, pmin, pmax, equ_min, equ_max, gm
        type(newton_state):: mem
! --------------------------------------------------------------------------------------------------
        if (DBG) write (6, *) 'BEGIN Solve_g_hat ', ts, cvg, g0

        ! Test des valeurs aux bornes
        pmax = bnd_pmax(self, ts, g0)
        p = pmax
        equ_max = f_g_hat(self, pmax, ts)-g0
        if (DBG) write (6, *) 'pmax, equmax = ', pmax, equ_max
        if (abs(equ_max) .le. cvg) goto 999
        ASSERT(equ_max .ge. 0)

        pmin = bnd_pmin(self, ts, g0)
        p = pmin
        equ_min = f_g_hat(self, pmin, ts)-g0
        if (abs(equ_min) .le. cvg) goto 999
        ASSERT(equ_min .le. 0)

        ASSERT(pmin .le. pmax)

        if (DBG) then
            write (6, *) 'pmin = ', pmin, '  pmax = ', pmax
            write (6, *) 'gmin-g0=', equ_min, 'gmax-g0=', equ_max
        end if

        ! Point initial
        if (present(p_ini)) then
            p = max(p_ini, pmin)
            p = min(p, pmax)
        else
            ! Cord method
            gm = -equ_min/(equ_max-equ_min)
            p = (1-gm)*pmin+gm*pmax
            p = min(max(p, pmin), pmax)
        end if

        if (DBG) write (6, *) 'utnewt 1B'
        do ite = 1, self%itemax
            equ = f_g_hat(self, p, ts)-g0
            d_equ = dp_g_hat(self, p, ts)
            if (abs(equ) .le. cvg) exit
            if (DBG) write (6, *) p, equ, d_equ, ite
            p = utnewt(p, equ, d_equ, ite, mem, xmin=pmin, xmax=pmax)
        end do
        if (DBG) write (6, *) 'utnewt 1E'
        if (ite .gt. self%itemax) then
            self%exception = 1
            goto 999
        end if
999     continue
        if (DBG) write (6, *) 'END Solve_g_hat'
    end function Solve_g_hat

! --------------------------------------------------------------------------------------------------
! Minimal bound pmin for equation G_hat(p,ts) = g0
! --------------------------------------------------------------------------------------------------

    function bnd_pmin(self, ts, g0) result(pmin)
        implicit none
        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8), intent(in)::ts, g0
        real(kind=8)::pmin
! --------------------------------------------------------------------------------------------------
        real(kind=8):: re, al, a, l0, b, c0, fmin, p0, dmin, coef, sh, gmin
! --------------------------------------------------------------------------------------------------

        if (DBG) write (6, *) 'BEGIN bnd_pmin'
        ! G(0,ts)-g0 < 0
        ASSERT(self%loaded)
        ASSERT(g0 .ge. -(1-self%dam)**2)

        ! Initialisation
        re = abs(self%telh/ts)
        al = self%mat%troismu/self%mat%troisk
        a = al*re
        l0 = 0.5d0*self%mat%q2*self%dam
        b = 1.5d0*self%mat%q2*re
        c0 = (4*al)/(3*self%mat%q2**2*self%dam)

        if (DBG) then
            write (6, *) 'ts   = ', ts
            write (6, *) 'g0   = ', g0
            write (6, *) 'telh = ', self%telh
            write (6, *) 're   = ', re
            write (6, *) 'al   = ', al
            write (6, *) 'l0   = ', l0
            write (6, *) 'b    = ', b
            write (6, *) 'c0   = ', c0
        end if

        ! Denominateur minimal
        if (c0 .le. 1) then
            fmin = 0.d0
        else if (acosh(c0) .ge. b) then
            fmin = l0*sinh(b)-a
        else
            p0 = acosh(c0)/b
            fmin = l0*sinh(b*p0)-a*p0
        end if
        if (DBG) write (6, *) 'fmin = ', fmin

        dmin = a+fmin

        coef = ((0.5d0*self%telq*self%mat%q2*self%dam)/(ts*dmin))**2+self%dam
        sh = sqrt((g0+(1-self%dam)**2)/coef)
        pmin = asinh(sh)/b
        if (DBG) write (6, *) 'pmin=', pmin

        if (DBG .and. pmin .ge. 1) then
            dbg_tmp1 = f_g_hat(self, 0.99999999d0, ts)
            write (6, *) 'g(1)=', dbg_tmp1
        end if

        ASSERT(pmin .lt. 1.d0)
        gmin = f_g_hat(self, pmin, ts)
        if (DBG) write (6, *) 'gmin=', gmin

        ASSERT(gmin .le. g0)

        if (DBG) write (6, *) 'END bnd_pmin'
    end function bnd_pmin

! --------------------------------------------------------------------------------------------------
! Maximal bound pmax for equation G_hat(p,ts) = g0
! --------------------------------------------------------------------------------------------------

    function bnd_pmax(self, ts, g0) result(pmax)
        implicit none
        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8), intent(in)::ts, g0
        real(kind=8)::pmax
! --------------------------------------------------------------------------------------------------
        real(kind=8):: chm, arg, p
! --------------------------------------------------------------------------------------------------

        if (DBG) write (6, *) 'BEGIN bnd_pmax'
        ! G(0,ts)-g0 < 0
        ASSERT(self%loaded)
        ASSERT(g0 .ge. -(1-self%dam)**2)

        chm = (g0+(1-self%dam)**2)/(2*self%dam)
        arg = acosh(1+chm)
        p = (2.d0*ts*arg)/(3*self%mat%q2*abs(self%telh))
        pmax = min(1.d0, p)

        if (DBG) then
            dbg_tmp1 = f_g(self, 0.d0, p, ts)
            dbg_tmp2 = f_g_hat(self, pmax, ts)
            write (6, *) 'p = ', p
            write (6, *) 'pmax = ', pmax
            write (6, *) 'gapp-g0 = ', dbg_tmp1-g0
            write (6, *) 'gmax-g0 = ', dbg_tmp2-g0
        end if

        if (DBG) write (6, *) 'END bnd_pmax'
    end function bnd_pmax

! ==================================================================================================
!  Solve von Mises problem (with plastic flow)
! ==================================================================================================

    subroutine Solve_Mises(self, cv_ts, dkaz, flow_type, dka)

        implicit none

        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8), intent(in)  :: cv_ts, dkaz
        integer(kind=8), intent(out)      :: flow_type
        real(kind=8), intent(out) :: dka
! --------------------------------------------------------------------------------------------------
! cv_ts         precision souhaitee sur T(ka) + c0 + c1*ka = 0
! dkaz          plastic increment such as T(dkaz)=0
! flow_type     0 = elastic, 1 = regular flow, 2 = singular flow
! dka           solution: increment of the hardening variable
! --------------------------------------------------------------------------------------------------
        aster_logical:: buffer
        integer(kind=8)     :: ite
        real(kind=8):: dkasvm, equ, d_equ, c0, c1, tsm, tss, kv, kvmin, kvmax
        type(newton_state):: mem
! --------------------------------------------------------------------------------------------------
        if (DBG) write (6, *) 'BEGIN Solve_Mises'
        ASSERT(self%loaded)

        ! Initialisation
        buffer = self%vgtn%norton%active

        ! Elastic solution
        tsm = f_ts_hat(self, dka=0.d0)
        if (self%telq/(1-self%dam) .le. tsm+cv_ts) then
            dka = 0.d0
            flow_type = 0
            goto 999
        end if

        ! Coefficient du terme affine
        c1 = self%mat%troismu/self%jac/(1-self%dam)**2
        c0 = -self%telq/(1-self%dam)

        ! Valeur de transition regulier / singulier pour von Mises
        dkasvm = -c0/c1
        tss = f_ts_hat(self, dka=dkasvm)

        ! Regime singulier pour von Mises (GTN regulier)
        if (tss .le. cv_ts) then
            flow_type = 2
            dka = dkaz
            goto 999
        end if

        ! Regime d'ecoulement regulier
        flow_type = 1

        ! Viscosite active (chgt de variable) ou non
        if (self%vgtn%is_visc .and. (self%vgtn%ts+c0+c1*self%vgtn%dka) .gt. 0) then
            self%vgtn%norton%active = ASTER_TRUE
            kv = 0.d0
            kvmin = 0.d0
            kvmax = min(self%vgtn%vsc, f_vsc(self%vgtn%norton, dka=dkasvm))
        else if (self%vgtn%is_visc) then
            self%vgtn%norton%active = ASTER_FALSE
            kv = self%vgtn%dka
            kvmin = self%vgtn%dka
            kvmax = dkasvm
        else
            self%vgtn%norton%active = ASTER_FALSE
            kv = -(tsm+c0)/(tss-tsm-c0)*dkasvm
            kvmin = 0.d0
            kvmax = dkasvm
        end if
        if (DBG) write (6, *) 'kvmin=', kvmin, '  kv=', kv, &
            'kvmax=', kvmax, '  v_active=', self%vgtn%norton%active
        ASSERT(kvmin .le. kv .and. kv .le. kvmax)

        do ite = 1, self%itemax
            dka = f_dka(self%vgtn%norton, kv)
            equ = f_ts_hat(self, kv=kv)+c0+c1*dka
            if (abs(equ) .le. cv_ts) exit
            d_equ = dkv_ts_hat(self, kv)+c1*dkv_dka(self%vgtn%norton, kv)
            if (DBG) write (6, *) 'utnewt 6B'
            kv = utnewt(kv, equ, d_equ, ite, mem, xmin=kvmin, xmax=kvmax)
            if (DBG) write (6, *) 'utnewt 6E'
        end do
        if (ite .gt. self%itemax) then
            self%exception = 1
            goto 999
        end if

999     continue
        self%vgtn%norton%active = buffer
        if (DBG) write (6, *) 'END Solve_Mises'
    end subroutine Solve_Mises

! ==================================================================================================
!  Compute the GTN equivalent norm of a stress tensor (by lower values)
! ==================================================================================================

    function Compute_Star(self, p, q, cv_ts) result(ts)

        implicit none

        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8), intent(in):: cv_ts, p, q
        real(kind=8)  ::ts
! --------------------------------------------------------------------------------------------------
! p             normalised hydrostatic stress TH = p*TelH
! q             normalised von Mises stress   TQ = q*TelQ
! cv_ts         expected accuracy on the residual
! return ts     T_star
! --------------------------------------------------------------------------------------------------
! Methode:
! On cherche la racine x --> G(p,q,1/x) fonction croissante, convexe et non bornee
! On assure in fine G(p,q,ts) > 0 et G(p,q,ts+cv_ts) < 0
! Une attention particuliere est portee aux questions d'arrondis numeriques
! --------------------------------------------------------------------------------------------------
        integer(kind=8)     :: ite
        real(kind=8):: tsmin, tsmax, xmin, xmax, gmin, gmax
        real(kind=8):: x, equ, d_equ, ts_m, ts_p, equ2, dts
        real(kind=8):: a, b, d, c1, c2
        type(newton_state):: mem
! --------------------------------------------------------------------------------------------------
        if (DBG) write (6, *) 'BEGIN Compute_Star'
        ASSERT(self%loaded)
        if (DBG) write (6, *) 'th=', p*self%telh, '  tq=', q*self%telq, '  cv_ts=', cv_ts

        ! Initialisation
        a = q*self%telq
        b = 1.5d0*self%mat%q2*abs(p*self%telh)
        d = self%dam
        ASSERT(1-d .gt. 0)
        if (DBG) write (6, *) 'a,b,d = ', a, b, d

        ! Cas d'un tenseur nul (a=b=0)
        if (b .eq. 0 .and. a .eq. 0) then
            ts = 0
            goto 999
        end if

        ! Cas d'un tenseur purement deviatorique (b=0) ou d'un endommagement nul (d=0)
        if (b .eq. 0 .or. d .eq. 0) then
            ts = a/(1-d)
            goto 999
        end if

        ! Cas d'un tenseur purement hydrostatique (a=0)
        if (a .eq. 0) then
            ts = b/acosh(1+(1-d)**2/(2*d))
            goto 999
        end if

        ! Borne min sur ts (borne max sur x = 1/ts) ---> G(tsmin) > 0  (rappel: G(0)=+inf)
        c1 = (1-d)/sqrt(a**2+d*b**2)
        c2 = acosh(1+(1-d)**2/(2*d))/b
        xmax = min(c1, c2)
        if (xmax .eq. 0.d0) then
            ASSERT(c1 .ne. 0.d0)
            c2 = (1-d)/(b*sqrt(d))
            xmax = min(c1, c2)
        end if
        ASSERT(xmax .gt. 0.d0)
        tsmin = 1.d0/xmax
        gmax = f_g(self, p, q, tsmin)
        if (DBG) write (6, *) 'xmax,tsmin,gmax : ', xmax, tsmin, gmax

        ! La borne min sur Ts est-elle une solution
        ts_m = tsmin
        ts_p = ts_m+cv_ts
        if (gmax .ge. 0) then
            if (f_g(self, p, q, ts_p) .le. 0) goto 500
        end if

        ! Correction de la borne tsmin pour des raisons d'arrondis si gmax=G(tsmin) < 0
        if (gmax .le. 0) then

            ! Recherche de la borne min sur ts par dichotomie et decalage (car G(0)=+inf)
            dts = cv_ts
            do ite = 1, 10*self%itemax
                ts_p = tsmin
                tsmin = max(0.5d0*tsmin, tsmin-dts)
                dts = dts*1.d1
                equ = f_g(self, p, q, tsmin)
                if (equ .ge. 0) exit
            end do
            if (ite .gt. self%itemax) then
                self%exception = 1
                goto 999
            else if (ts_p-tsmin .le. cv_ts) then
                ts_m = tsmin
                goto 500
            end if
            xmax = 1.d0/tsmin
            gmax = equ
        end if
        ASSERT(tsmin .gt. 0.d0)

        ! Borne max sur ts (borne min sur x=1/ts) --> G(tsmax) < 0
        xmin = acosh(1+(1-d)**2/(2*(d+(a/b)**2)))/b
        if (xmin .eq. 0.d0) xmin = (1-d)/(2*b*sqrt(d+(a/b)**2))
        ASSERT(xmin .gt. 0.d0)
        tsmax = 1.d0/xmin
        gmin = f_g(self, p, q, tsmax)
        if (DBG) write (6, *) 'xmin,tsmax,gmin : ', xmin, tsmax, gmin

        ! La borne max sur ts est-elle solution ?
        if (gmin .le. 0) then
            ts_p = tsmax
            if (tsmax-tsmin .le. cv_ts) then
                ts_m = tsmin
                goto 500
            end if
            ts_m = ts_p-cv_ts
            if (f_g(self, p, q, ts_m) .ge. 0) goto 500
        end if

        ! Correction de la borne tsmax en cas d'erreur d'arrondis
        if (gmin .ge. 0) then

            ! Recherche de la borne max sur ts par amplification d'un facteur 2 et decalage
            dts = cv_ts
            do ite = 1, 10*self%itemax
                ts_m = tsmax
                tsmax = min(tsmax+dts, 2*tsmax)
                dts = dts*1.d1
                equ = f_g(self, p, q, tsmax)
                if (equ .le. 0) exit
            end do
            if (ite .gt. self%itemax) then
                self%exception = 1
                goto 999
            else if (tsmax-ts_m .le. cv_ts) then
                ts_p = tsmax
                goto 500
            end if
            xmin = 1.d0/tsmax
            gmin = equ
        end if

        if (DBG) write (6, *) 'xmax, tsmin, gmax = ', xmax, tsmin, gmax
        if (DBG) write (6, *) 'xmin, tsmax, gmin = ', xmin, tsmax, gmin

        ASSERT(xmin .le. xmax)
        ASSERT(gmin .le. gmax)
        ASSERT(tsmin .le. tsmax)

        ! Resolution G(p,q,1/x)=0 pour tirer profit de la convexite de G par rapport a x
        ts = tsmin
        x = xmax
        do ite = 1, self%itemax
            equ = f_g(self, p, q, ts)
            if (DBG) write (6, *) 'ite = ', ite, '   x = ', x, '   ts = ', ts, '   equ = ', equ

            ! Construction de l'encadrement et test (monotonie de la fonction)
            if (equ .ge. 0) then
                ts_m = ts
                ts_p = ts+cv_ts
                if (ts_p .ge. tsmax) goto 500
                equ2 = f_g(self, p, q, ts_p)
                if (DBG) write (6, *) 'ts_p=', ts_p, '  equ2=', equ2
                if (equ2 .le. 0) goto 500
            else
                ts_p = ts
                ts_m = ts-cv_ts
                if (ts_m .le. tsmin) then
                    ts_m = tsmin
                    goto 500
                end if
                equ2 = f_g(self, p, q, ts_m)
                if (DBG) write (6, *) 'ts_m=', ts_m, '  equ2=', equ2
                if (equ2 .ge. 0) goto 500
            end if

            ! Itere de Newton
            d_equ = -dts_g(self, p, q, ts)*ts**2
            if (DBG) write (6, *) 'utnewt 7B'
            x = utnewt(x, equ, d_equ, ite, mem, xmin=xmin, xmax=xmax)
            ts = 1.d0/x
            if (DBG) write (6, *) 'utnewt 7E'
        end do
        if (ite .gt. self%itemax) then
            self%exception = 1
            goto 999
        end if

        ! On fournit la borne inferieure de l'encadrement (inferieure a la valeur reelle)
500     continue
        ts = ts_m

999     continue
        if (DBG) write (6, *) 'END Compute_Star ', self%exception
    end function Compute_Star

! ==================================================================================================
!  Viscous parameters required by the coupling between damage and viscosity
! ==================================================================================================

    subroutine Set_Visc_Damage(self, dam)

        implicit none

        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8), intent(in):: dam
! --------------------------------------------------------------------------------------------------
! dam   damage level
! --------------------------------------------------------------------------------------------------
        real(kind=8):: gam
! --------------------------------------------------------------------------------------------------

        if (DBG) write (6, *) 'BEGIN Set_Visc_Damage'

        self%vgtn%is_visc = self%mat%norton%visc

        if (self%vgtn%is_visc) then
            self%vgtn%norton = self%mat%norton
            if (dam .le. self%mat%dv) then
                gam = 1
            else
                gam = (1-dam)/(1-self%mat%dv)
            end if
            self%vgtn%norton%v0 = self%mat%norton%v0/gam
            self%vgtn%norton%k = self%mat%norton%k/gam
        end if

        if (DBG) write (6, *) 'END Set_Visc_Damage'
    end subroutine Set_Visc_Damage

! ==================================================================================================
!  Viscous parameters for the algorithmic transition between weak and strong viscosity
! ==================================================================================================

    subroutine Set_Visc_Transition(self)

        implicit none

        type(CONSTITUTIVE_LAW), intent(inout):: self
! --------------------------------------------------------------------------------------------------
        real(kind=8):: dka_r0
! --------------------------------------------------------------------------------------------------

        if (DBG) write (6, *) 'BEGIN Set_Visc_Transition'

        if (self%vgtn%is_visc) then
            dka_r0 = dka_ecro(self, 0.d0)
            self%vgtn%dka = solve_slope_dka(self%vgtn%norton, dka_r0)
            self%vgtn%vsc = f_vsc(self%vgtn%norton, dka=self%vgtn%dka)
            self%vgtn%ts = f_ts_hat(self, dka=self%vgtn%dka)
        end if

        if (DBG) write (6, *) 'END Set_Visc_Transition'
    end subroutine Set_Visc_Transition

! ==================================================================================================
!  SET POROSITY AND DAMAGE VALUES
! ==================================================================================================

    function Set_Porosity_Damage(self, poro_in, poro_nucl, dam_in, epcum, dam_rate) result(pd)

        implicit none

        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8), intent(in):: poro_in, poro_nucl, dam_in, epcum, dam_rate
        type(PORO_DAM):: pd
! --------------------------------------------------------------------------------------------------
! poro_in   porosity
! poro_nucl nucleation porosity
! dam_in    damage
! epcum     cumulated plastic strain
! dam_rate  damage rate
! --------------------------------------------------------------------------------------------------
! Invariants to fullfil:
!       poro = poro_nucl + poro_grow
!       dam  = dam_coal + q1*poro
! --------------------------------------------------------------------------------------------------
        real(kind=8):: poro, poro_grow, dam, dam_coal
! --------------------------------------------------------------------------------------------------

        if (DBG) write (6, *) 'BEGIN Set_Porosity_Damage'

        ! Local variables
        poro = poro_in
        dam = dam_in

        ! Find a compatible state for porosity even though initial poro is not acceptable
        poro_grow = poro-poro_nucl
        if (poro_grow .lt. self%mat%f0) then
            poro_grow = self%mat%f0
            poro = poro_grow+poro_nucl
        end if

        ! Find a compatible state for damage even though initial dam is not acceptable
        dam_coal = dam-self%mat%q1*poro
        if (dam_coal .lt. 0.d0) then
            dam_coal = 0.d0
            dam = self%mat%q1*poro
        end if

        ! Store the variables
        pd%poro = poro
        pd%poro_nucl = poro_nucl
        pd%poro_grow = poro_grow
        pd%dam = dam
        pd%dam_coal = dam_coal
        pd%epcum = epcum
        pd%dam_rate = dam_rate

        if (DBG) write (6, *) 'END Set_Porosity_Damage'
    end function Set_Porosity_Damage

! ==================================================================================================
!  COMPUTATION OF THE STAR-POROSITY (INCLUDING COALESCENCE)
! ==================================================================================================

    function Star_Porosity(self, poro) result(star_poro)

        implicit none

        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8), intent(in):: poro
        real(kind=8)           :: star_poro
! --------------------------------------------------------------------------------------------------
! poro      porosity
! star_poro augmented porosity with coalescence
! --------------------------------------------------------------------------------------------------
        if (DBG) write (6, *) 'BEGIN Star_Porosity'
        star_poro = poro+self%mat%hc*max(0.d0, poro-self%mat%fc)
        if (DBG) write (6, *) 'END Star_Porosity'
    end function Star_Porosity

! ==================================================================================================
!  EXTRAPOLATION OF DAMAGE WITH RESPECT OF BOUNDS
! ==================================================================================================

    function Extrapolate_Damage(self, pd_m) result(dam_extr)

        implicit none

        type(CONSTITUTIVE_LAW), intent(inout):: self
        type(PORO_DAM), intent(in):: pd_m
        real(kind=8):: dam_extr
! --------------------------------------------------------------------------------------------------
! pd_m      initial porosity and damage
! --------------------------------------------------------------------------------------------------
        real(kind=8):: dam_min, dam_max
! --------------------------------------------------------------------------------------------------
        if (DBG) write (6, *) 'BEGIN Extrapolate_Damage'
        dam_min = self%mat%q1*(self%mat%f0+pd_m%poro_nucl)+pd_m%dam_coal
        dam_max = 1.d0
        dam_extr = pd_m%dam+self%theta*self%dt*pd_m%dam_rate
        dam_extr = min(max(dam_min, dam_extr), dam_max)
        if (DBG) write (6, *) 'END Extrapolate_Damage'
    end function Extrapolate_Damage

! ==================================================================================================
!  COMPUTATION OF NUCLEATION POROSITY
! ==================================================================================================

    real(kind=8) function Nucleation(self, ka, ep)

        implicit none

        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8), intent(in)::ka, ep
! --------------------------------------------------------------------------------------------------
! ka        hardening variable
! ep        cumulated plastic strain
! --------------------------------------------------------------------------------------------------
        real(kind=8):: fgauss, fcran, feps
! --------------------------------------------------------------------------------------------------
        if (DBG) write (6, *) 'BEGIN Nucleation'

        fgauss = 0.5d0*self%mat%fn*(erf((ka-self%mat%pn)/sqrt(2.d0)/self%mat%sn)+ &
                                    erf(self%mat%pn/sqrt(2.d0)/self%mat%sn))
        fcran = max(0.d0, min(self%mat%c0/(self%mat%kf-self%mat%ki)*(ka-self%mat%ki), self%mat%c0))
        feps = self%mat%b0*max(ep-self%mat%epc, 0.d0)
        Nucleation = fgauss+fcran+feps

        if (DBG) write (6, *) 'END Nucleation'
    end function Nucleation

! ==================================================================================================
!  UPDATE POROSITY AND DAMAGE (NUCLEATION AND GROWTH)
! ==================================================================================================
    function Update_Porosity_Damage(self, state, epm, ep, ka, pd_m) result(pd)
        implicit none
        type(CONSTITUTIVE_LAW), intent(inout):: self
        integer(kind=8), intent(in):: state
        real(kind=8), intent(in):: epm(:), ep(:), ka
        type(PORO_DAM), intent(in):: pd_m
        type(PORO_DAM):: pd
! --------------------------------------------------------------------------------------------------
! state     plastic state
! epm       initial plastic strain
! ep        final plastic strain
! ka        final hardening variable
! pd_m      initial porosity and damage
! pd        final (updated) porosity and damage
! --------------------------------------------------------------------------------------------------
        real(kind=8):: coef, dep(size(ep)), deptr
        real(kind=8):: poro, poro_nucl, poro_grow, dam, epcum, dam_rate, poro_nucl_th, dam_coal
! --------------------------------------------------------------------------------------------------
        if (state .eq. BROKEN_STATE) then
            epcum = pd_m%epcum
            poro_nucl = pd_m%poro_nucl
            poro = self%mat%fr
            dam = 1.d0
            goto 800
        end if

        dep = ep-epm
        deptr = sum(dep(1:3))

        ! Nucleation
        epcum = pd_m%epcum+sqrt(2.d0/3.d0*dot_product(dep, dep))
        poro_nucl = Nucleation(self, ka, epcum)

        ! Growth law (theta scheme)
        poro_nucl_th = self%theta*poro_nucl+(1-self%theta)*pd_m%poro_nucl
        poro_grow = (pd_m%poro_grow+deptr*(1-poro_nucl_th-(1-self%theta)*pd_m%poro_grow)) &
                    /(1+self%theta*deptr)

        ! Constraint poro_grow >= f0
        poro_grow = max(self%mat%f0, poro_grow)

        ! Constraint poro <= fr
        if (poro_nucl+poro_grow .gt. self%mat%fr) then
            coef = (self%mat%fr-pd_m%poro)/(poro_nucl+poro_grow-pd_m%poro)
            poro_nucl = coef*poro_nucl+(1-coef)*pd_m%poro_nucl
            poro_grow = coef*poro_grow+(1-coef)*pd_m%poro_grow
        end if

        ! Porosity
        poro = poro_nucl+poro_grow
        poro = min(self%mat%fr, poro)

        ! Coalescence Damage evolution (prop: dam<=1)
        dam_coal = max(pd_m%dam_coal, self%mat%q1*(Star_Porosity(self, poro)-poro))
        dam = dam_coal+self%mat%q1*poro
        dam = min(1.d0, dam)

800     continue
        dam_rate = (dam-pd_m%dam)/self%dt
        pd = Set_Porosity_Damage(self, poro, poro_nucl, dam, epcum, dam_rate)

    end function Update_Porosity_Damage

! ==================================================================================================
!  IS IT BROKEN ?
! ==================================================================================================

    aster_logical function is_broken(self, dam)

        implicit none

        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8), intent(in)::dam
! --------------------------------------------------------------------------------------------------
        if (DBG) write (6, *) 'Start is_broken'
        is_broken = dam .ge. self%mat%dam_bkn
        if (DBG) write (6, *) 'End is_broken'
    end function is_broken

! ==================================================================================================
!  RESPONSE OF A BROKEN POINT
! ==================================================================================================

    subroutine Broken_Point(self, eps, sigm, state, dka, ep, &
                            t, deps_t, dphi_t, deps_ka, dphi_ka)

        implicit none

        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8), intent(in)             :: eps(:), sigm(:)
        integer(kind=8), intent(out)                 :: state
        real(kind=8), intent(out)            :: dka, ep(:), t(:)
        real(kind=8), intent(out)            :: deps_t(:, :), dphi_t(:), deps_ka(:), dphi_ka
! --------------------------------------------------------------------------------------------------
! eps       final strain
! sigm      initial stress
! state     regime during time-step
! dka       increment of hardening variable
! ep        final plastic strain
! t         final stress
! deps_t    derivate dt / deps
! dphi_t    derivate dt / dphi   (grad_vari)
! deps_ka   derivate dka / deps  (grad_vari)
! dphi_ka   derivate dka / dphi  (grad_vari)
! --------------------------------------------------------------------------------------------------
        real(kind=8):: n, cv, a, dt0, um, red
! --------------------------------------------------------------------------------------------------
        if (DBG) write (6, *) 'Start Broken_Point'
        state = BROKEN_STATE
        dka = 0.d0
        ep = eps

        ! Instantaneous stress vanishing
        if (.not. self%mat%norton%visc) then
            red = 0
            goto 500
        end if

        ! The stress vaninshing is damped
        n = self%mat%norton%n
        cv = self%mat%norton%k*(1.d0-self%mat%dv)

        ! quasi zero damping
        if (cv .eq. 0) then
            red = 0
            goto 500
        end if

        a = self%mat%young/cv
        um = sqrt(dot_product(sigm, sigm))/cv

        ! Accelerated regime
        if (um .le. 1 .or. n .eq. 1) then
            red = exp(-a*self%dt)
            goto 500
        end if

        dt0 = (1-um**(1-n))/((n-1)*a)

        ! Mixed regime
        if (dt0 .le. self%dt) then
            red = exp(-a*(self%dt-dt0))/um
            goto 500
        end if

        ! Norton regime
        red = (1+a*(n-1)*self%dt*um**(n-1))**(1/(1-n))

500     continue
        if (DBG) write (6, *) 'red = ', red
        ASSERT(red .ge. 0)
        ASSERT(red .lt. 1)
        t = red*sigm

        if (self%rigi) then
            deps_t = 0.d0
            dphi_t = 0.d0
            deps_ka = 0.d0
            dphi_ka = 0.d0
        end if

        if (DBG) write (6, *) 'End Broken_Point'
    end subroutine Broken_Point

! ==================================================================================================
!  POST-TREATMENT
! ==================================================================================================

    function Compute_Post(self, pd_m, pd, dam_extr, t, dka) &
        result(post)

        implicit none
        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8), intent(in):: t(:), dam_extr, dka
        type(PORO_DAM), intent(in):: pd_m, pd
        type(POST_TREATMENT):: post
! --------------------------------------------------------------------------------------------------
! porom     initial porosity
! dam_extr  extrapolated damage
! poro      final porosity
! t         final stress
! dka       increment of internal variable
! --------------------------------------------------------------------------------------------------
        real(kind=8):: cv_sigm, cv_ts, att, tn, sig_star_extr, sig_star_updt
        type(CONSTITUTIVE_LAW)    :: temp
! --------------------------------------------------------------------------------------------------
        if (DBG) write (6, *) 'BEGIN Compute_Post'

        ! Stress split
        if (.not. is_broken(self, dam_extr)) then
            ASSERT(self%loaded)
            att = self%jac*(1-dam_extr)
            post%sieq_ecr = att*f_ecro(self, self%kam+dka)
            if (self%grvi) post%sieq_nlc = att*(self%mat%r*(self%kam+dka)-self%phi)
            if (self%vgtn%is_visc) post%sieq_vsc = att*f_vsc(self%vgtn%norton, dka=dka)
        end if

        ! Error due to explicit scheme measured in terms of stress value
        tn = sqrt(dot_product(t, t))
        cv_sigm = self%mat%sig0*self%cvuser

        ! Broken law during time step (no stress error)
        if (is_broken(self, dam_extr) .or. is_broken(self, pd%dam)) then
            post%sieq_erx = 0.d0

            ! Small stress implies small stress error
        else if (tn .le. cv_sigm) then
            post%sieq_erx = 0.d0

            ! Predictor vs. update error
        else
            ASSERT(self%loaded)
            cv_ts = 1.d-1*cv_sigm/(1-dam_extr)

            temp = self
            temp%telh = sum(t(1:3))/3.d0
            temp%telq = sqrt(1.5d0*abs(tn**2-3*temp%telh**2))

            temp%dam = dam_extr
            sig_star_extr = Compute_Star(temp, 1.d0, 1.d0, cv_ts)
            ASSERT(self%exception .eq. 0)
            if (DBG) write (6, *) 'sig_star_extr = ', sig_star_extr

            temp%dam = pd%dam
            sig_star_updt = Compute_Star(temp, 1.d0, 1.d0, cv_ts)
            ASSERT(self%exception .eq. 0)
            if (DBG) write (6, *) 'sig_star_updt = ', sig_star_updt

            post%sieq_erx = abs(1-sig_star_extr/sig_star_updt)*tn
        end if

        ! User stopping criteria
        post%arret = merge(1, 0, is_broken(self, pd_m%dam))

        if (DBG) write (6, *) 'END Compute_Post'
    end function Compute_Post

! ==================================================================================================
!  Liste des fonctions intermediaires et leurs derivees
!   --> Derivees: djac_*, dteh_*, dteq_*, dphi_*
! ==================================================================================================

    real(kind=8) function f_g(self, q, p, ts)
        implicit none
        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8)::q, p, ts
        ASSERT(self%loaded)
        f_g = (self%telq*q/ts)**2+2*self%dam*cosh(1.5d0*self%mat%q2*self%telh*p/ts)-1-self%dam**2
    end function f_g

    real(kind=8) function dq_g(self, q, p, ts)
        implicit none
        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8)::q, p, ts
        ASSERT(self%loaded)
        dq_g = 2*self%telq**2*q/ts**2
    end function dq_g

    real(kind=8) function dp_g(self, q, p, ts)
        implicit none
        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8)::q, p, ts
        ASSERT(self%loaded)
        if (DBG) write (6, *) 'dp_g: ', q, p, ts
        dp_g = 3*self%dam*self%mat%q2*self%telh/ts*sinh(1.5d0*self%mat%q2*self%telh*p/ts)
        if (DBG) write (6, *) 'dp_g ok'
    end function dp_g

    real(kind=8) function dts_g(self, q, p, ts)
        implicit none
        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8)::q, p, ts
        ASSERT(self%loaded)
        if (DBG) write (6, *) 'dts_g: ', q, p, ts
        dts_g = -2*(self%telq*q)**2/ts**3-3*self%dam*self%mat%q2*self%telh*p/ts**2 &
                *sinh(1.5d0*self%mat%q2*self%telh*p/ts)
        if (DBG) write (6, *) 'dts_g ok'
    end function dts_g

    real(kind=8) function dteh_g(self, q, p, ts)
        implicit none
        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8)::q, p, ts
        ASSERT(self%loaded)
        if (DBG) write (6, *) 'dteh_g: ', q, p, ts
        dteh_g = 3*self%dam*self%mat%q2*p/ts*sinh(1.5d0*self%mat%q2*self%telh*p/ts)
        if (DBG) write (6, *) 'dteh_g ok'
    end function dteh_g

    real(kind=8) function dteq_g(self, q, p, ts)
        implicit none
        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8)::q, p, ts
        ASSERT(self%loaded)
        dteq_g = 2*self%telq*(q/ts)**2
    end function dteq_g

    real(kind=8) function ddam_g(self, q, p, ts)
        implicit none
        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8)::q, p, ts
        ASSERT(self%loaded)
        ddam_g = 2*cosh(1.5d0*self%mat%q2*self%telh*p/ts)-2*self%dam

    end function ddam_g

    real(kind=8) function f_ecro(self, ka)
        implicit none
        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8)::ka
        f_ecro = self%mat%r0+self%mat%rh*ka+self%mat%r1*(1-exp(-self%mat%g1*ka)) &
                 +self%mat%r2*(1-exp(-self%mat%g2*ka))+self%mat%rk*(ka+self%mat%p0)**self%mat%gk
    end function f_ecro

    real(kind=8) function dka_ecro(self, ka)
        implicit none
        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8)::ka
        dka_ecro = self%mat%rh+self%mat%r1*self%mat%g1*exp(-self%mat%g1*ka) &
                   +self%mat%r2*self%mat%g2*exp(-self%mat%g2*ka) &
                   +self%mat%rk*self%mat%gk*(ka+self%mat%p0)**(self%mat%gk-1)
    end function dka_ecro

    real(kind=8) function f_ts_hat(self, dka, kv)
        implicit none
        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8), intent(in), optional::dka, kv
        real(kind=8):: vsc, ka

        ! one and one only among dka and kv is given
        ASSERT(present(dka) .eqv. .not. present(kv))
        ASSERT(self%loaded)

        if (present(dka)) then
            vsc = f_vsc(self%vgtn%norton, dka=dka)
            ka = self%kam+dka
        else
            vsc = f_vsc(self%vgtn%norton, kv=kv)
            ka = self%kam+f_dka(self%vgtn%norton, kv)
        end if

        f_ts_hat = self%jac*(f_ecro(self, ka)+self%mat%r*ka-self%phi+vsc)
    end function f_ts_hat

    real(kind=8) function dkv_ts_hat(self, kv)
        implicit none
        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8)::kv
        real(kind=8):: ka
        ASSERT(self%loaded)

        ka = self%kam+f_dka(self%vgtn%norton, kv)
        dkv_ts_hat = self%jac*((dka_ecro(self, ka)+self%mat%r)*dkv_dka(self%vgtn%norton, kv) &
                               +dkv_vsc(self%vgtn%norton, kv))
    end function dkv_ts_hat

    real(kind=8) function dka_ts_hat(self, dka)
        implicit none
        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8)::dka
        ASSERT(self%loaded)

        dka_ts_hat = self%jac*(dka_ecro(self, self%kam+dka)+ &
                               self%mat%r+ddka_vsc(self%vgtn%norton, dka))
    end function dka_ts_hat

    real(kind=8) function djac_ts_hat(self, dka)
        implicit none
        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8)::dka
        real(kind=8)::ka
        ASSERT(self%loaded)
        ka = self%kam+dka
        djac_ts_hat = f_ecro(self, ka)+self%mat%r*ka-self%phi+f_vsc(self%vgtn%norton, dka)
    end function djac_ts_hat

    real(kind=8) function dphi_ts_hat(self, dka)
        implicit none
        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8)::dka
        ASSERT(self%loaded)
        dphi_ts_hat = -self%jac
    end function dphi_ts_hat

    real(kind=8) function f_lambda(self, x)
        implicit none
        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8)::x
        ASSERT(self%loaded)
        if (DBG) write (6, *) 'f_lambda: ', x
        f_lambda = 0.5d0*self%dam*self%mat%q2*sinh(1.5d0*self%mat%q2*x)
        if (DBG) write (6, *) 'f_lambda ok'
    end function f_lambda

    real(kind=8) function dx_lambda(self, x)
        implicit none
        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8)::x
        ASSERT(self%loaded)
        dx_lambda = 0.75d0*self%dam*self%mat%q2**2*cosh(1.5d0*self%mat%q2*x)
    end function dx_lambda

    real(kind=8) function f_theta(self, x, y)
        implicit none
        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8)::x, y
        if (DBG) write (6, *) 'f_theta'
        f_theta = 1.d0/(y**2+3*x*f_lambda(self, x))
    end function f_theta

    real(kind=8) function dx_theta(self, x, y)
        implicit none
        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8)::x, y
        real(kind=8)::lbd
        if (DBG) write (6, *) 'dx_theta'
        lbd = f_lambda(self, x)
        dx_theta = -3*(lbd+x*dx_lambda(self, x))/(y**2+3*x*lbd)**2
    end function dx_theta

    real(kind=8) function dy_theta(self, x, y)
        implicit none
        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8)::x, y
        if (DBG) write (6, *) 'dy_theta'
        dy_theta = -(2*y)/(y**2+3*x*f_lambda(self, x))**2
    end function dy_theta

    real(kind=8) function f_q_hat(self, p, ts)
        implicit none
        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8)::p, ts
        real(kind=8)::muk, num, den
        ASSERT(self%loaded)
        if (DBG) write (6, *) 'f_q_hat: ', p, ts, self%telh
        muk = self%mat%troismu/self%mat%troisk
        num = ts*f_lambda(self, self%telh*p/ts)
        den = num+muk*self%telh*(1-p)
        f_q_hat = num/den
    end function f_q_hat

    real(kind=8) function dp_q_hat(self, p, ts)
        implicit none
        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8)::p, ts
        real(kind=8)::muk, num, den, d_num, d_den
        ASSERT(self%loaded)
        if (DBG) write (6, *) 'dp_q_hat'
        muk = self%mat%troismu/self%mat%troisk
        num = ts*f_lambda(self, self%telh*p/ts)
        den = num+muk*self%telh*(1-p)
        d_num = dx_lambda(self, self%telh*p/ts)*self%telh
        d_den = d_num-muk*self%telh
        dp_q_hat = (d_num*den-num*d_den)/den**2
    end function dp_q_hat

    real(kind=8) function dts_q_hat(self, p, ts)
        implicit none
        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8)::p, ts
        real(kind=8)::muk, num, den, lbd, d_num, d_den
        ASSERT(self%loaded)
        if (DBG) write (6, *) 'dts_q_hat'
        muk = self%mat%troismu/self%mat%troisk
        lbd = f_lambda(self, self%telh*p/ts)
        num = ts*lbd
        den = num+muk*self%telh*(1-p)
        d_num = lbd-dx_lambda(self, self%telh*p/ts)*(self%telh*p/ts)
        d_den = d_num
        dts_q_hat = (d_num*den-num*d_den)/den**2
    end function dts_q_hat

    real(kind=8) function dteh_q_hat(self, p, ts)
        implicit none
        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8)::p, ts
        real(kind=8)::muk, w, num, d_num, den, d_den, lbd0, a0p, seuil
        ASSERT(self%loaded)
        if (DBG) write (6, *) 'dteh_q_hat'
        muk = self%mat%troismu/self%mat%troisk
        w = 1.5d0*self%mat%q2*self%telh*p/ts
        seuil = (r8prem()*1200)**0.25d0
        if (abs(w) .le. seuil) then
            ! Precision numerique quand telh petit (dvp limite de lambda ordre 4)
            lbd0 = 0.5d0*self%dam*self%mat%q2
            a0p = 1.5d0*self%mat%q2*p
            dteh_q_hat = (12*(a0p)**3*lbd0*muk*(1-p)*ts**2*self%telh) &
                         /(lbd0*a0p*(6*ts**2+(a0p*self%telh)**2)+6*muk*ts**2*(1-p))**2
        else
            ! Derivee normale de lambda avec le sinh
            num = ts*f_lambda(self, self%telh*p/ts)
            d_num = dx_lambda(self, self%telh*p/ts)*p
            den = num+muk*self%telh*(1-p)
            d_den = d_num+muk*(1-p)
            dteh_q_hat = (d_num*den-num*d_den)/den**2
        end if
    end function dteh_q_hat

    real(kind=8) function f_g_hat(self, p, ts)
        implicit none
        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8)::p, ts
        real(kind=8)::q
        if (DBG) write (6, *) 'f_g_hat'
        q = f_q_hat(self, p, ts)
        f_g_hat = f_g(self, q, p, ts)
    end function f_g_hat

    real(kind=8) function dp_g_hat(self, p, ts)
        implicit none
        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8)::p, ts
        real(kind=8)::q
        q = f_q_hat(self, p, ts)
        dp_g_hat = dq_g(self, q, p, ts)*dp_q_hat(self, p, ts)+dp_g(self, q, p, ts)
    end function dp_g_hat

    real(kind=8) function dts_g_hat(self, p, ts)
        implicit none
        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8)::p, ts
        real(kind=8)::q
        q = f_q_hat(self, p, ts)
        dts_g_hat = dq_g(self, q, p, ts)*dts_q_hat(self, p, ts)+dts_g(self, q, p, ts)
    end function dts_g_hat

    real(kind=8) function dteh_g_hat(self, p, ts)
        implicit none
        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8)::p, ts
        real(kind=8)::q, d_q
        q = f_q_hat(self, p, ts)
        d_q = dteh_q_hat(self, p, ts)
        dteh_g_hat = dq_g(self, q, p, ts)*d_q+dteh_g(self, q, p, ts)
    end function dteh_g_hat

    real(kind=8) function dteq_g_hat(self, p, ts)
        implicit none
        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8)::p, ts
        real(kind=8)::q
        q = f_q_hat(self, p, ts)
        dteq_g_hat = dteq_g(self, q, p, ts)
    end function dteq_g_hat

    real(kind=8) function f_lbd_hat(self, p, ts)
        implicit none
        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8)::p, ts
        real(kind=8)::jk, q, num, den, lbd, the
        ASSERT(self%loaded)
        if (DBG) write (6, *) 'f_lbd_hat'
        jk = self%jac/self%mat%troisk
        q = f_q_hat(self, p, ts)
        lbd = f_lambda(self, self%telh*p/ts)
        the = f_theta(self, self%telh*p/ts, self%telq*q/ts)
        num = jk*self%telh*(1-p)
        den = the*lbd
        f_lbd_hat = num/den
    end function f_lbd_hat

    real(kind=8) function dp_lbd_hat(self, p, ts)
        implicit none
        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8)::p, ts
        real(kind=8)::jk, q, d_q, lbd, d_lbd, the, d_the, num, d_num, den, d_den
        ASSERT(self%loaded)
        if (DBG) write (6, *) 'dp_lbd_hat'
        jk = self%jac/self%mat%troisk
        q = f_q_hat(self, p, ts)
        d_q = dp_q_hat(self, p, ts)
        lbd = f_lambda(self, self%telh*p/ts)
        d_lbd = dx_lambda(self, self%telh*p/ts)*self%telh/ts
        the = f_theta(self, self%telh*p/ts, self%telq*q/ts)
        d_the = dx_theta(self, self%telh*p/ts, self%telq*q/ts)*self%telh/ts &
                +dy_theta(self, self%telh*p/ts, self%telq*q/ts)*self%telq/ts*d_q
        num = jk*self%telh*(1-p)
        d_num = -jk*self%telh
        den = the*lbd
        d_den = d_the*lbd+the*d_lbd
        dp_lbd_hat = (d_num*den-num*d_den)/den**2
    end function dp_lbd_hat

    real(kind=8) function dts_lbd_hat(self, p, ts)
        implicit none
        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8)::p, ts
        real(kind=8)::jk, q, d_q, lbd, d_lbd, the, d_the, num, den, d_den
        ASSERT(self%loaded)
        if (DBG) write (6, *) 'dts_lbd_hat'
        jk = self%jac/self%mat%troisk
        q = f_q_hat(self, p, ts)
        d_q = dts_q_hat(self, p, ts)
        lbd = f_lambda(self, self%telh*p/ts)
        d_lbd = dx_lambda(self, self%telh*p/ts)*(-self%telh*p/ts**2)
        the = f_theta(self, self%telh*p/ts, self%telq*q/ts)
        d_the = dx_theta(self, self%telh*p/ts, self%telq*q/ts)*(-self%telh*p/ts**2) &
                +dy_theta(self, self%telh*p/ts, self%telq*q/ts)*self%telq*(d_q*ts-q)/ts**2
        num = jk*self%telh*(1-p)
        den = the*lbd
        d_den = d_the*lbd+the*d_lbd
        dts_lbd_hat = -num*d_den/den**2
    end function dts_lbd_hat

    real(kind=8) function djac_lbd_hat(self, p, ts)
        implicit none
        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8)::p, ts
        real(kind=8)::jk, d_jk, q, num, d_num, den, lbd, the
        if (DBG) write (6, *) 'djac_lbd_hat'
        ASSERT(self%loaded)
        jk = self%jac/self%mat%troisk
        d_jk = 1/self%mat%troisk
        q = f_q_hat(self, p, ts)
        lbd = f_lambda(self, self%telh*p/ts)
        the = f_theta(self, self%telh*p/ts, self%telq*q/ts)
        num = jk*self%telh*(1-p)
        d_num = d_jk*self%telh*(1-p)
        den = the*lbd
        djac_lbd_hat = d_num/den
    end function djac_lbd_hat

    real(kind=8) function dteh_lbd_hat(self, p, ts)
        implicit none
        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8)::p, ts
        real(kind=8)::jk, q, d_q, num, d_num, den, d_den, lbd, d_lbd, the, d_the
        ASSERT(self%loaded)
        if (DBG) write (6, *) 'dteh_lbd_hat'
        jk = self%jac/self%mat%troisk
        q = f_q_hat(self, p, ts)
        d_q = dteh_q_hat(self, p, ts)
        lbd = f_lambda(self, self%telh*p/ts)
        d_lbd = dx_lambda(self, self%telh*p/ts)*p/ts
        the = f_theta(self, self%telh*p/ts, self%telq*q/ts)
        d_the = dx_theta(self, self%telh*p/ts, self%telq*q/ts)*p/ts &
                +dy_theta(self, self%telh*p/ts, self%telq*q/ts)*self%telq*d_q/ts
        num = jk*self%telh*(1-p)
        d_num = jk*(1-p)
        den = the*lbd
        d_den = d_the*lbd+the*d_lbd
        dteh_lbd_hat = (d_num*den-num*d_den)/den**2
    end function dteh_lbd_hat

    real(kind=8) function dteq_lbd_hat(self, p, ts)
        implicit none
        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8)::p, ts
        real(kind=8)::jk, q, num, den, d_den, lbd, the, d_the
        if (DBG) write (6, *) 'dteq_lbd_hat'
        ASSERT(self%loaded)
        jk = self%jac/self%mat%troisk
        q = f_q_hat(self, p, ts)
        lbd = f_lambda(self, self%telh*p/ts)
        the = f_theta(self, self%telh*p/ts, self%telq*q/ts)
        d_the = dy_theta(self, self%telh*p/ts, self%telq*q/ts)*q/ts
        num = jk*self%telh*(1-p)
        den = the*lbd
        d_den = d_the*lbd
        dteq_lbd_hat = -num*d_den/den**2
    end function dteq_lbd_hat

    real(kind=8) function dp_m_hat(self, p, ts, dka)
        implicit none
        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8)::p, ts, dka
        dp_m_hat = -dp_lbd_hat(self, p, ts)
    end function dp_m_hat

    real(kind=8) function dts_m_hat(self, p, ts, dka)
        implicit none
        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8)::p, ts, dka
        dts_m_hat = 1.d0/dka_ts_hat(self, dka)-dts_lbd_hat(self, p, ts)
    end function dts_m_hat

    real(kind=8) function djac_m_hat(self, p, ts, dka)
        implicit none
        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8)::p, ts, dka
        djac_m_hat = -djac_ts_hat(self, dka)/dka_ts_hat(self, dka)-djac_lbd_hat(self, p, ts)
    end function djac_m_hat

    real(kind=8) function dphi_m_hat(self, p, ts, dka)
        implicit none
        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8)::p, ts, dka
        dphi_m_hat = -dphi_ts_hat(self, dka)/dka_ts_hat(self, dka)
    end function dphi_m_hat

    real(kind=8) function dteh_m_hat(self, p, ts, dka)
        implicit none
        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8)::p, ts, dka
        dteh_m_hat = -dteh_lbd_hat(self, p, ts)
    end function dteh_m_hat

    real(kind=8) function dteq_m_hat(self, p, ts, dka)
        implicit none
        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8)::p, ts, dka
        dteq_m_hat = -dteq_lbd_hat(self, p, ts)
    end function dteq_m_hat

end module lcgtn_module
