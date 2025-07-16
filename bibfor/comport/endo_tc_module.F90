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

module endo_tc_module

    use scalar_newton_module, only: &
        newton_state, &
        utnewt

    use tenseur_dime_module, only: &
        kron, proten

    use endo_crit_tc_module, only: &
        CRITERION, &
        InitCrit => Init, &
        DerivativeCrit => Derivative

    use endo_rigi_unil_module, only: &
        UNILATERAL, &
        UNIL_MAT => MATERIAL, &
        InitUnilateral => Init, &
        ComputeStress, &
        ComputeStress_eig, &
        ComputeStiffness, &
        ComputeEnergy

    implicit none
    private
    public:: CONSTITUTIVE_LAW, Init, Integrate

#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utmess.h"
#include "asterc/r8pi.h"
#include "asterc/r8prem.h"
#include "asterc/r8vide.h"

! --------------------------------------------------------------------------------------------------

    ! Variables internes
    integer(kind=8), parameter:: ENDO = 1
    integer(kind=8), parameter:: ENDOTRAC = 2
    integer(kind=8), parameter:: HISTTRAC = 3
    integer(kind=8), parameter:: ENERTRAC = 4
    integer(kind=8), parameter:: ENDOCOMP = 5
    integer(kind=8), parameter:: HISTCOMP = 6
    integer(kind=8), parameter:: ENERCOMP = 7
    integer(kind=8), parameter:: SIGMVISC = 8
    integer(kind=8), parameter:: ENDOTOT = 9

    ! Material parameters
    type MATERIAL
        real(kind=8)   :: lambda
        real(kind=8)   :: deuxmu
        real(kind=8)   :: omega_bar
        real(kind=8)   :: m0
        real(kind=8)   :: d1
        real(kind=8)   :: r
        real(kind=8)   :: ft
        real(kind=8)   :: sig0
        real(kind=8)   :: fc
        real(kind=8)   :: beta
        real(kind=8)   :: crit_p
        real(kind=8)   :: coef_v
        type(UNIL_MAT) :: unil
    end type MATERIAL

    ! Required precision and impact on the different quantities
    type PRECISION
        real(kind=8):: cvuser
        real(kind=8):: cvsig
        real(kind=8):: cveps
        real(kind=8):: cvdam
    end type PRECISION

    ! Shared attibutes through the global variable self
    type CONSTITUTIVE_LAW
        integer(kind=8)         :: exception = 0
        aster_logical   :: elas, rigi, vari, pilo, pred, visc
        integer(kind=8)         :: ndimsi, itemax
        real(kind=8)    :: deltat
        type(PRECISION) :: prec
        type(MATERIAL)  :: mat
        type(UNILATERAL):: unil
        type(CRITERION) :: critTens, critComp
    end type CONSTITUTIVE_LAW

    ! Post-treatment results
    type POST_TREATMENT
        real(kind=8)  :: tens_stf
        real(kind=8)  :: comp_stf
        real(kind=8)  :: wpos
        real(kind=8)  :: wneg
        real(kind=8)  :: deq
        real(kind=8)  :: sigmvisc
    end type POST_TREATMENT

contains

! ==================================================================================================
!  OBJECT CREATION AND INITIALISATION
! ==================================================================================================

    function Init(ndimsi, option, fami, kpg, ksp, imate, deltat, itemax, precvg) result(self)

        implicit none

        integer(kind=8), intent(in)          :: kpg, ksp, imate, itemax, ndimsi
        real(kind=8), intent(in)    :: deltat, precvg
        character(len=16), intent(in):: option
        character(len=*), intent(in) :: fami
        type(CONSTITUTIVE_LAW)      :: self
! --------------------------------------------------------------------------------------------------
! ndimsi    symmetric tensor dimension (2*ndim)
! option    computation option
! fami      Gauss point set
! kpg       Gauss point number
! ksp       Layer number (for structure elements)
! imate     material pointer
! deltat    time increment
! itemax    max number of iterations for the solver
! precvg    required accuracy (with respect to stress level))
! --------------------------------------------------------------------------------------------------

        self%rigi = option .eq. 'RIGI_MECA_TANG' .or. option .eq. 'RIGI_MECA_ELAS' &
                    .or. option .eq. 'FULL_MECA' .or. option .eq. 'FULL_MECA_ELAS'
        self%vari = option .eq. 'FULL_MECA_ELAS' .or. option .eq. 'FULL_MECA' &
                    .or. option .eq. 'RAPH_MECA'
        self%pilo = option .eq. 'PILO_PRED_ELAS'
        self%elas = option .eq. 'RIGI_MECA_ELAS' .or. option .eq. 'FULL_MECA_ELAS'
        self%pred = .not. (self%vari .or. self%pilo)

        ASSERT(merge(1, 0, self%pred)+merge(1, 0, self%vari)+merge(1, 0, self%pilo) .eq. 1)
        ASSERT(.not. self%pred .or. self%rigi)

        self%ndimsi = ndimsi
        self%deltat = deltat
        self%itemax = itemax
        self%mat = GetMaterial(self, fami, kpg, ksp, imate, precvg)

        ! Expected accuracy for the stress, the strain and the damage
        self%prec%cvuser = precvg
        self%prec%cvsig = self%mat%ft*precvg
        self%prec%cveps = self%prec%cvsig/(self%mat%unil%lambda+self%mat%unil%deuxmu)
        self%prec%cvdam = precvg

    end function Init

! ==================================================================================================
!  INTEGRATION OF ENDO_LOCA_TC (MAIN ROUTINE)
! ==================================================================================================

    subroutine Integrate(self, eps, vim, sig, vip, dsde)

        implicit none

        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8), intent(in)    :: eps(:), vim(:)
        real(kind=8), intent(out)    :: sig(:), vip(:), dsde(:, :)
! ----------------------------------------------------------------------
! eps   current strain
! vim   internal variables at the beginning of the time step
! sig   stress at the end of the time step
! vip   internal variables at the end of the time step
! dsde  tangent matrix (ndimsi,ndimsi)
! --------------------------------------------------------------------------------------------------
        real(kind=8)        :: bem, cmaxm, be, cmax
        type(POST_TREATMENT):: post
! --------------------------------------------------------------------------------------------------

! unpack internal variables
        bem = vim(HISTTRAC)
        cmaxm = vim(HISTCOMP)

! Strain split tension/compression and strain measures in tension and compression
        call ComputeStrainQuantities(self, eps)

! damage behaviour integration
        call ComputeLaw(self, bem, cmaxm, be, cmax, sig, dsde)
        if (self%exception .ne. 0) goto 999

! Post-treatments
        if (self%vari) post = PostTreatment(self, be, cmax)

! pack internal variables
        if (self%vari) then
            vip(ENDO) = merge(post%deq, vim(1), post%deq .ne. r8vide())
            vip(ENDOTRAC) = 1.d0-post%tens_stf
            vip(HISTTRAC) = be
            vip(ENERTRAC) = post%wpos
            vip(ENDOCOMP) = 1.d0-post%comp_stf
            vip(HISTCOMP) = cmax
            vip(ENERCOMP) = post%wneg
            vip(SIGMVISC) = post%sigmvisc
            vip(ENDOTOT) = 1.d0-post%tens_stf*post%comp_stf
        end if

999     continue
    end subroutine Integrate

! ==================================================================================================
!  MATERIAL CHARACTERISTICS
! ==================================================================================================

    function GetMaterial(self, fami, kpg, ksp, imate, precvg) result(mat)

        implicit none

        type(CONSTITUTIVE_LAW), intent(inout):: self
        integer(kind=8), intent(in)                   :: kpg, ksp, imate
        character(len=*), intent(in)          :: fami
        real(kind=8), intent(in)              :: precvg
        type(MATERIAL)                       :: mat
! --------------------------------------------------------------------------------------------------
! fami      Gauss point set
! kpg       Gauss point number
! ksp       Layer number (for structure elements)
! imate     material pointer
! precvg    required precision
! --------------------------------------------------------------------------------------------------
        integer(kind=8), parameter:: nbel = 2, nben = 5, nbrg = 1
! --------------------------------------------------------------------------------------------------
        integer(kind=8) :: iokel(nbel), ioken(nben), iokrg(nbrg)
        real(kind=8) :: valel(nbel), valen(nben), valrg(nbrg)
        character(len=16) :: nomel(nbel), nomen(nben), nomrg(nbrg)
        real(kind=8):: pi
        real(kind=8):: cvsig, coef_trac, coef_comp
! --------------------------------------------------------------------------------------------------
        data nomel/'E', 'NU'/
        data nomen/'ENER_TRAC_RUPT_N', 'COEF_ECRO_TRAC', 'FT', 'SIGM_COMP_SEUIL', 'FC'/
        data nomrg/'TAU_REGU_VISC'/
! --------------------------------------------------------------------------------------------------

        pi = r8pi()

!  Elasticity
        call rcvalb(fami, kpg, ksp, '+', imate, ' ', 'ELAS', 0, ' ', [0.d0], &
                    nbel, nomel, valel, iokel, 2)

!   Damage parameters
        call rcvalb(fami, kpg, ksp, '+', imate, ' ', 'ENDO_LOCA_TC', 0, ' ', [0.d0], &
                    nben, nomen, valen, ioken, 2)

!   Regularisation parameters
        call rcvalb(fami, kpg, ksp, '+', imate, ' ', 'ENDO_LOCA_TC', 0, ' ', [0.d0], &
                    nbrg, nomrg, valrg, iokrg, 2)

        mat%lambda = valel(1)*valel(2)/((1+valel(2))*(1-2*valel(2)))
        mat%deuxmu = valel(1)/(1+valel(2))

        mat%omega_bar = valen(1)
        mat%m0 = 1.5d0*pi*(valen(2)+2)**(-1.5d0)
        mat%d1 = 0.75d0*pi*sqrt(1+valen(2))
        mat%r = (2*(mat%d1-1)-mat%m0)/(2-mat%m0)

        mat%ft = valen(3)
        mat%sig0 = valen(4)
        mat%fc = valen(5)
        mat%beta = mat%fc/mat%sig0

        mat%coef_v = 0.d0
        self%visc = valrg(1) .gt. 0.d0
        if (self%visc) then
            if (self%deltat < valrg(1)*r8prem()) then
                call utmess('F', 'COMPOR1_96', sr=self%deltat/valrg(1))
            end if
            mat%coef_v = valrg(1)/self%deltat
        end if

        mat%unil%lambda = mat%lambda
        mat%unil%deuxmu = mat%deuxmu

        ! Admissibilite du jeu de parametres
        if (mat%omega_bar*mat%m0 .le. 2) call utmess('F', 'COMPOR1_98', sr=mat%m0)
        ASSERT(mat%ft .gt. 0.d0)
        ASSERT(mat%omega_bar .gt. 0.d0)
        ASSERT(mat%m0 .gt. 0.d0)
        ASSERT(mat%m0 .lt. 2.d0)
        ASSERT(mat%d1 .gt. mat%m0)
        ASSERT(mat%r .gt. 1.d0)
        ASSERT(mat%sig0 .gt. 2*mat%ft)
        ASSERT(mat%beta .gt. 1.d0)
        ASSERT(mat%coef_v .ge. 0.d0)

        ! Calcul des parametres de regularisation en fonction de la precision souhaitee
        cvsig = mat%ft*precvg
        mat%unil%regbet = cvsig/(mat%lambda+mat%deuxmu)
        coef_trac = (mat%m0*mat%omega_bar-2)/(mat%m0*mat%omega_bar)*(cvsig/mat%ft)
        coef_comp = cvsig/mat%fc
        mat%crit_p = log(3.d0)*max((1+coef_trac)/coef_trac, (1+coef_comp)/coef_comp)

    end function GetMaterial

! ==================================================================================================
!  COMPUTE STRAIN QUANTITIES (STRAIN SPLIT AND STRAIN MEASURES)
! ==================================================================================================

    subroutine ComputeStrainQuantities(self, eps)

        implicit none

        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8), intent(in)             :: eps(:)

! --------------------------------------------------------------------------------------------------
! self  -> Compute self.unil (unilateral object), self.critTens and self.critComp (strain measures)
! eps   strain at the end of the time step
! --------------------------------------------------------------------------------------------------
        real(kind=8), dimension(self%ndimsi):: kr
        real(kind=8), dimension(self%ndimsi):: sigpos, signeg, tnsTens, tnsComp
! --------------------------------------------------------------------------------------------------

        ! Initialiation
        kr = kron(self%ndimsi)

        ! Split tension / compression
        self%unil = InitUnilateral(self%mat%unil, eps, self%prec%cveps)
        call ComputeStress(self%unil, sigpos, signeg)

        ! tensile damage measure
        tnsTens = (self%mat%lambda*sum(eps(1:3))*kr+self%mat%deuxmu*eps)/self%mat%ft
        self%critTens = InitCrit(self%mat%crit_p, tnsTens, self%prec%cvuser)

        ! Compressive strain measure
        tnsComp = -signeg/self%mat%sig0
        self%critComp = InitCrit(self%mat%crit_p, tnsComp, self%prec%cvuser)

    end subroutine ComputeStrainQuantities

! ==================================================================================================
!  DAMAGE AND STRESS COMPUTATION AND TANGENT OPERATOR (ACCORDING TO RIGI AND RESI)
! ==================================================================================================

    subroutine ComputeLaw(self, bem, cmaxm, be, cmax, sig, dsde)

        implicit none

        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8), intent(in)              :: bem, cmaxm
        real(kind=8), intent(out)             :: be, cmax
        real(kind=8), intent(out)             :: sig(:), dsde(:, :)
! --------------------------------------------------------------------------------------------------
! be    damage at the beginning (in) of the time-step then its end (out)
! cmax  maximal compressive driving variable
! sig   stress at the end of the time step
! dsde  tangent matrix (ndimsi,ndimsi)
! --------------------------------------------------------------------------------------------------
        integer(kind=8)     :: tens_state, comp_state
        real(kind=8), dimension(self%ndimsi):: kr
        real(kind=8), dimension(self%ndimsi):: sigpos, signeg
        real(kind=8), dimension(self%ndimsi, self%ndimsi):: deps_sigpos, deps_signeg
        real(kind=8):: quad, tens_stf, comp_stf, be_pert
        real(kind=8):: db_B, db_phi, dtns_chit(self%ndimsi), db_visc, dq_visc
        real(kind=8):: deps_chit(self%ndimsi), deps_b(self%ndimsi)
        real(kind=8):: deps_C(self%ndimsi), dtns_chic(self%ndimsi), deps_chic(self%ndimsi)
! --------------------------------------------------------------------------------------------------

        ! Initialiation
        kr = kron(self%ndimsi)

! --------------------------------------------------------------------------------------------------
!  Damage computation
! --------------------------------------------------------------------------------------------------

        ! Tensile damage
        call ComputeTensileDamage(self, bem, be, tens_state)
        if (self%exception .ne. 0) goto 999

        ! Compressive damage
        call ComputeCompressionDamage(self, cmaxm, cmax, comp_state)

! --------------------------------------------------------------------------------------------------
!  Stress computation
! --------------------------------------------------------------------------------------------------

        call ComputeStress(self%unil, sigpos, signeg)
        tens_stf = FB(self, be)
        comp_stf = FC(self, cmax)
        sig = comp_stf*(tens_stf*sigpos+signeg)

! --------------------------------------------------------------------------------------------------
!  Tangent matrix computation
! --------------------------------------------------------------------------------------------------

        if (self%rigi) then

            quad = self%critTens%chi**2

            ! Contribution elastique (non lineaire) a endommagement fixe
            call ComputeStiffness(self%unil, deps_sigpos, deps_signeg)
            dsde = comp_stf*(tens_stf*deps_sigpos+deps_signeg)

            ! Actualisation des regimes en phase de prediction (si valeurs proches des seuils)
            if (self%pred .and. .not. self%elas) then
                if (tens_state .eq. 0) then
                    be_pert = max(0.d0, be-self%prec%cvdam)
                    if (quad .ge. Fphi(self, be_pert)+Fvisc(self, be_pert-bem, quad)) tens_state = 1
                end if
                comp_state = merge(1, 0, self%critComp%chi .gt. max(1.d0, (cmaxm-self%prec%cvuser)))
            end if

            ! Contribution liee a la variation d'endommagement de traction
            if ((.not. self%elas) .and. (tens_state .eq. 1)) then
                dtns_chit = DerivativeCrit(self%critTens)
                deps_chit = (self%mat%lambda*sum(dtns_chit(1:3))*kr+self%mat%deuxmu*dtns_chit) &
                            /self%mat%ft
                db_phi = db_Fphi(self, be)
                db_visc = db_Fvisc(self, be-bem, quad)
                dq_visc = dquad_Fvisc(self, be-bem, quad)
                deps_b = (1-dq_visc)*2*self%critTens%chi*deps_chit/(db_phi+db_visc)
                db_B = db_FB(self, be)
                dsde = dsde+comp_stf*db_B*proten(sigpos, deps_b)
            end if

            ! Contribution liee a la variation d'endommagement de compression
            if ((.not. self%elas) .and. (comp_state .eq. 1)) then
                dtns_chic = DerivativeCrit(self%critComp)
                deps_chic = -matmul(deps_signeg, dtns_chic)/self%mat%sig0
                deps_C = Dx_FC(self, self%critComp%chi)*deps_chic
                dsde = dsde+proten(tens_stf*sigpos+signeg, deps_C)
            end if

        end if

999     continue
    end subroutine ComputeLaw

! ==================================================================================================
!  TENSILE DAMAGE COMPUTATION
! ==================================================================================================

    subroutine ComputeTensileDamage(self, bem, be, tens_state)

        implicit none

        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8), intent(in)              :: bem
        real(kind=8), intent(out)             :: be
        integer(kind=8), intent(out)                  :: tens_state
! --------------------------------------------------------------------------------------------------
! self
! bem   tensile damage at the beginning of the time step
! be    tensile damage at the end of the time step
! tens_state  state of damage
!               0 = elastic
!               1 = damage
!               2 = saturated
! --------------------------------------------------------------------------------------------------
        integer(kind=8)     :: iter
        real(kind=8), dimension(self%ndimsi):: sigpos, signeg
        real(kind=8):: quad, equ, db_equ, phi, visc
        real(kind=8):: bemin, equmin, bemax, equmax, signrm, dsig, dbemax, becpt, equcpt
        type(newton_state):: mem
! --------------------------------------------------------------------------------------------------

        ! Nonlinear elastic stresses and tangent (or secant) matrices
        call ComputeStress(self%unil, sigpos, signeg)

        ! tensile damage measure
        quad = self%critTens%chi**2

        ! Already saturated point
        if (bem .ge. 1.d0) then
            tens_state = 2
            be = 1.d0
            goto 999
        end if

        ! Elastic regime
        bemin = bem
        equmin = Fphi(self, bemin)+Fvisc(self, bemin-bem, quad)-quad
        if (equmin .ge. 0.d0) then
            tens_state = 0
            be = bem
            goto 999
        end if

        ! Saturated point
        bemax = 1.d0
        equmax = Fphi(self, bemax)+Fvisc(self, bemax-bem, quad)-quad
        if (equmax .le. 0.d0) then
            tens_state = 2
            be = 1.d0
            goto 999
        end if

        ! Damage regime
        tens_state = 1
        signrm = sqrt(dot_product(sigpos, sigpos))

        ! Secant estimation (lower bound thanks to the convexity of phi)
        be = (bemin*equmax-bemax*equmin)/(equmax-equmin)

        do iter = 1, self%itemax
            phi = Fphi(self, be)
            visc = Fvisc(self, be-bem, quad)
            equ = phi+visc-quad

            ! Convergence check: interval smaller than the accuracy with change of equation sign

            ! Required accuracy
            dsig = abs(db_FB(self, be))*signrm
            if (self%prec%cvdam*dsig .le. self%prec%cvsig) then
                dbemax = self%prec%cvdam
            else
                dbemax = self%prec%cvsig/dsig
            end if

            ! Search interval
            becpt = merge(min(bemax, be+dbemax), max(bemin, be-dbemax), equ .le. 0.d0)
            equcpt = Fphi(self, becpt)+Fvisc(self, becpt-bem, quad)-quad

            ! Sign change
            if (equ*equcpt .le. 0.d0) exit

            ! New iterate
            db_equ = Db_Fphi(self, be)+Db_Fvisc(self, be-bem, quad)
            be = utnewt(be, equ, db_equ, iter, mem, xmin=bemin, xmax=bemax)
        end do
        if (iter .gt. self%itemax) then
            self%exception = 1
            goto 999
        end if

999     continue
    end subroutine ComputeTensileDamage

! ==================================================================================================
!  COMPRESSION DAMAGE COMPUTATION
! ==================================================================================================

    subroutine ComputeCompressionDamage(self, cmaxm, cmax, comp_state)

        implicit none

        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8), intent(in)              :: cmaxm
        real(kind=8), intent(out)             :: cmax
        integer(kind=8), intent(out)                  :: comp_state
! --------------------------------------------------------------------------------------------------
! self
! bem   tensile damage at the beginning of the time step
! be    tensile damage at the end of the time step
! tens_state  state of damage
!               0 = elastic
!               1 = damage
! --------------------------------------------------------------------------------------------------

        comp_state = merge(1, 0, self%critComp%chi .gt. max(1.d0, cmaxm))
        cmax = max(cmaxm, self%critComp%chi)

    end subroutine ComputeCompressionDamage

! ==================================================================================================
!  POST-TREATMENTS
! ==================================================================================================

    function PostTreatment(self, be, cmax) result(post)

        implicit none
        type(CONSTITUTIVE_LAW)    :: self
        real(kind=8), intent(in)  :: be
        real(kind=8), intent(in)  :: cmax
        type(POST_TREATMENT)      :: post
! --------------------------------------------------------------------------------------------------
! be    damage at the end of the time-step
! cmax  maximal compressive driving variable
! --------------------------------------------------------------------------------------------------
        integer(kind=8):: state_free
        real(kind=8):: we, be_free, tens_stf_free, sigpos_eig(3), signeg_eig(3)
! --------------------------------------------------------------------------------------------------

        ! stiffness damage
        post%tens_stf = FB(self, be)
        post%comp_stf = FC(self, cmax)

        ! elastic energy
        call ComputeEnergy(self%unil, post%wpos, post%wneg)

        ! Directional damage
        we = post%wpos+post%wneg
        if (we .le. 0.d0) then
            post%deq = r8vide()
        else
            post%deq = 1.d0-post%comp_stf*(post%tens_stf*post%wpos+post%wneg)/we
            post%deq = max(0.d0, min(1.d0, post%deq))
        end if

        ! Viscous stress
        post%sigmvisc = 0
        if (self%visc) then
            self%visc = .FALSE.
            call ComputeStress_eig(self%unil, sigpos_eig, signeg_eig)
            call ComputeTensileDamage(self, be, be_free, state_free)
            tens_stf_free = FB(self, be_free)
            post%sigmvisc = post%comp_stf*(post%tens_stf-tens_stf_free)*sigpos_eig(1)
            self%visc = .TRUE.
        end if
    end function PostTreatment

! ==================================================================================================
!  Intermediate function h(b)
! ==================================================================================================

    function Fh(self, b) result(res)

        implicit none
        type(CONSTITUTIVE_LAW), intent(in):: self
        real(kind=8):: res
        real(kind=8):: b
! ---------------------------------------------------------------------
        real(kind=8):: c1, cr
! ---------------------------------------------------------------------

        c1 = 0.5d0*self%mat%omega_bar*self%mat%m0-1
        cr = 0.5d0*self%mat%omega_bar*(self%mat%d1-self%mat%m0)
        res = 1+c1*b+cr*b**self%mat%r

    end function Fh
! --------------------------------------------------------------------------------------------------

    function Db_Fh(self, b) result(res)

        implicit none
        type(CONSTITUTIVE_LAW), intent(in):: self
        real(kind=8):: res
        real(kind=8):: b
! --------------------------------------------------------------------------------------------------
        real(kind=8):: c1, cr
! --------------------------------------------------------------------------------------------------

        c1 = 0.5d0*self%mat%omega_bar*self%mat%m0-1
        cr = 0.5d0*self%mat%omega_bar*(self%mat%d1-self%mat%m0)
        res = c1+cr*self%mat%r*b**(self%mat%r-1)

    end function Db_Fh

! ==================================================================================================
!  Stiffness function FB(b)
! ==================================================================================================

    function FB(self, b) result(res)

        implicit none
        type(CONSTITUTIVE_LAW), intent(in):: self
        real(kind=8):: res
        real(kind=8):: b
! --------------------------------------------------------------------------------------------------
        res = (1-b)/Fh(self, b)
    end function FB
! --------------------------------------------------------------------------------------------------

    function Db_FB(self, b) result(res)

        implicit none
        type(CONSTITUTIVE_LAW), intent(in):: self
        real(kind=8):: res
        real(kind=8):: b
! --------------------------------------------------------------------------------------------------
        real(kind=8):: h, db_h
! --------------------------------------------------------------------------------------------------
        h = Fh(self, b)
        db_h = Db_Fh(self, b)
        res = -(h+(1-b)*db_h)/h**2

    end function Db_FB

! ==================================================================================================
!  Threshold function Phi(b)
! ==================================================================================================

    function Fphi(self, b) result(res)

        implicit none
        type(CONSTITUTIVE_LAW), intent(in):: self
        real(kind=8):: res
        real(kind=8):: b
! --------------------------------------------------------------------------------------------------
        res = Fh(self, b)**2

    end function Fphi
! --------------------------------------------------------------------------------------------------

    function Db_Fphi(self, b) result(res)

        implicit none
        type(CONSTITUTIVE_LAW), intent(in):: self
        real(kind=8):: res
        real(kind=8):: b
! --------------------------------------------------------------------------------------------------
        res = 2*Db_Fh(self, b)*Fh(self, b)

    end function Db_Fphi

! ==================================================================================================
!  Viscous function Fvisc(b)
! ==================================================================================================

    function Fvisc(self, delta_b, quad) result(res)

        implicit none
        type(CONSTITUTIVE_LAW), intent(in):: self
        real(kind=8):: res
        real(kind=8):: delta_b, quad
! --------------------------------------------------------------------------------------------------
        res = merge(self%mat%coef_v*quad*delta_b, 0.d0, self%visc)

    end function Fvisc
! --------------------------------------------------------------------------------------------------

    function Db_Fvisc(self, delta_b, quad) result(res)

        implicit none
        type(CONSTITUTIVE_LAW), intent(in):: self
        real(kind=8):: res
        real(kind=8):: delta_b, quad
! --------------------------------------------------------------------------------------------------
        res = merge(self%mat%coef_v*quad, 0.d0, self%visc)

    end function Db_Fvisc
! --------------------------------------------------------------------------------------------------

    function Dquad_Fvisc(self, delta_b, quad) result(res)

        implicit none
        type(CONSTITUTIVE_LAW), intent(in):: self
        real(kind=8):: res
        real(kind=8):: delta_b, quad
! --------------------------------------------------------------------------------------------------
        res = merge(self%mat%coef_v*delta_b, 0.d0, self%visc)

    end function Dquad_Fvisc
! --------------------------------------------------------------------------------------------------

! ==================================================================================================
!  Compression function C(x)
! ==================================================================================================

    function FC(self, x) result(res)

        implicit none
        type(CONSTITUTIVE_LAW), intent(in):: self
        real(kind=8):: res
        real(kind=8):: x
! --------------------------------------------------------------------------------------------------
        if (x .le. 1.d0) then
            res = 1.d0
        else
            res = (1.d0+(self%mat%beta-1.d0)*tanh((x-1.d0)/(self%mat%beta-1.d0)))/x
        end if

    end function FC
! --------------------------------------------------------------------------------------------------

    function Dx_FC(self, x) result(res)

        implicit none
        type(CONSTITUTIVE_LAW), intent(in):: self
        real(kind=8):: res
! --------------------------------------------------------------------------------------------------
        real(kind=8):: x, th, num, dx_num
! --------------------------------------------------------------------------------------------------

        if (x .le. 1.d0) then
            res = 0.d0
        else
            th = tanh((x-1.d0)/(self%mat%beta-1.d0))
            num = 1.d0+(self%mat%beta-1.d0)*th
            dx_num = 1-th**2
            res = (dx_num*x-num)/x**2
        end if

    end function Dx_FC
! --------------------------------------------------------------------------------------------------

end module endo_tc_module
