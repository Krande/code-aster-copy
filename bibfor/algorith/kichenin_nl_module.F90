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

module kichenin_nl_module

    use scalar_newton_module, only: &
        newton_state, &
        utnewt

    use tenseur_dime_module, only: &
        proten, &
        kron, &
        sph_norm, &
        deviator, &
        voigt, &
        identity

    implicit none
    private
    public:: CONSTITUTIVE_LAW, Init, Integrate

#include "asterf_types.h"
#include "asterc/r8gaem.h"
#include "asterfort/assert.h"
#include "asterfort/rcvalb.h"

! --------------------------------------------------------------------------------------------------

    ! Material characteristics

    type MATERIAL
        real(kind=8) :: troiskp, deuxmup, sc, prag
        real(kind=8) :: troiskv, deuxmuv, vsc, nud, ga
    end type MATERIAL

    ! VMIS_ISOT_NL class
    type CONSTITUTIVE_LAW
        integer(kind=8)       :: exception = 0
        aster_logical :: elas, rigi, resi, vari, pred, visc_lin
        integer(kind=8)       :: ndimsi, itemax
        real(kind=8)  :: cvuser, small
        type(MATERIAL):: mat
    end type CONSTITUTIVE_LAW

contains

! =====================================================================
!  OBJECT CREATION AND INITIALISATION
! =====================================================================

    function Init(ndimsi, option, fami, kpg, ksp, imate, itemax, precvg, deltat) &
        result(self)

        implicit none

        integer(kind=8), intent(in)          :: kpg, ksp, imate, itemax, ndimsi
        real(kind=8), intent(in)    :: precvg, deltat
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
! itemax    max number of iterations for the solver
! precvg    required accuracy (with respect to stress level)
! deltat    time-increment duration
! --------------------------------------------------------------------------------------------------
        integer(kind=8), parameter   :: nbel = 2, nbki = 7
! --------------------------------------------------------------------------------------------------
        integer(kind=8)             :: okel(nbel), okki(nbki)
        real(kind=8)        :: valel(nbel), valki(nbki)
        character(len=16)   :: nomel(nbel), nomki(nbki)
! --------------------------------------------------------------------------------------------------
        data nomel/'E', 'NU'/
        data nomki/'SIGC', 'PRAGER', 'E_VISC', 'NU_VISC', 'NU_AMOR', 'N_AMOR', 'ETA_AMOR'/
! --------------------------------------------------------------------------------------------------

        ! Parametres generaux
        self%ndimsi = ndimsi
        self%itemax = itemax
        self%cvuser = precvg
        self%small = 1.d0/r8gaem()

        ! Options de calcul
        self%elas = option .eq. 'RIGI_MECA_ELAS' .or. option .eq. 'FULL_MECA_ELAS'
        self%rigi = option .eq. 'RIGI_MECA_TANG' .or. option .eq. 'RIGI_MECA_ELAS' &
                    .or. option .eq. 'FULL_MECA' .or. option .eq. 'FULL_MECA_ELAS'
        self%resi = option .eq. 'FULL_MECA_ELAS' .or. option .eq. 'FULL_MECA' &
                    .or. option .eq. 'RAPH_MECA'
        self%vari = option .eq. 'FULL_MECA_ELAS' .or. option .eq. 'FULL_MECA' &
                    .or. option .eq. 'RAPH_MECA'
        self%pred = option .eq. 'RIGI_MECA_ELAS' .or. option .eq. 'RIGI_MECA_TANG'

        ASSERT(self%pred .or. self%resi)

        ! Elasticity material parameters for the plastic branch
        call rcvalb(fami, kpg, ksp, '+', imate, ' ', 'ELAS', 0, ' ', [0.d0], &
                    nbel, nomel, valel, okel, 2)
        self%mat%deuxmup = valel(1)/(1+valel(2))
        self%mat%troiskp = valel(1)/(1.d0-2.d0*valel(2))

        ! Other parameters
        call rcvalb(fami, kpg, ksp, '+', imate, ' ', 'KICHENIN_NL', 0, ' ', [0.d0], &
                    nbki, nomki, valki, okki, 2)
        self%mat%sc = valki(1)
        self%mat%prag = valki(2)
        self%mat%deuxmuv = valki(3)/(1+valki(4))
        self%mat%troiskv = valki(3)/(1.d0-2.d0*valki(4))
        self%mat%nud = valki(5)
        self%mat%ga = 1.d0/valki(6)
        self%mat%vsc = deltat/(valki(7)**self%mat%ga)

        ! Viscosite lineaire ou non
        self%visc_lin = (valki(6) .eq. 1.d0 .or. self%mat%ga .eq. 1.d0)

        ! Controle des parametres
        ASSERT(self%mat%sc .gt. 0)
        ASSERT(self%mat%sc .ge. 0)
        ASSERT(valki(3) .ge. 0)
        ASSERT(valki(4) .gt. -1.d0 .and. valki(4) .lt. 0.5d0)
        ASSERT(self%mat%nud .ge. -1.d0 .and. self%mat%nud .le. 0.5d0)
        ASSERT(self%mat%ga .ge. 1.d0)
        ASSERT(self%mat%vsc .ge. 0.d0)

    end function Init

! =====================================================================
!  INTEGRATION OF THE CONSTITUTIVE LAW (MAIN ROUTINE)
! =====================================================================

    subroutine Integrate(self, eps, vim, sig, vip, deps_sig)

        implicit none

        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8), intent(in)         :: eps(:), vim(:)
        real(kind=8), intent(out)         :: sig(:), vip(:), deps_sig(:, :)
! --------------------------------------------------------------------------------------------------
! eps       strain at the end of the time step
! vim       internal variables at the beginning of the time step
! sig       stress at the end of the time step (resi) or the beginning of the time step (not resi)
! vip       internal variables at the end of the time step
! deps_sig  derivee dsig / deps
! --------------------------------------------------------------------------------------------------
        real(kind=8)    :: kapl, kavs, ep(self%ndimsi), ev(self%ndimsi), rac2(self%ndimsi)
        real(kind=8)    :: sip(self%ndimsi), deps_sip(self%ndimsi, self%ndimsi)
        real(kind=8)    :: siv(self%ndimsi), deps_siv(self%ndimsi, self%ndimsi)
! --------------------------------------------------------------------------------------------------

! unpack internal variables
        rac2 = voigt(self%ndimsi)
        kapl = vim(1)
        ep = vim(2:1+self%ndimsi)*rac2
        kavs = vim(8)
        ev = vim(9:8+self%ndimsi)*rac2

! Plasticity integration
        call ComputePlasticity(self, eps, kapl, ep, sip, deps_sip)
        if (self%exception .ne. 0) goto 999

! Viscoelasticity integration
        call ComputeViscoElasticity(self, eps, kavs, ev, siv, deps_siv)
        if (self%exception .ne. 0) goto 999

! Combination
        sig = sip+siv
        if (self%rigi) deps_sig = deps_sip+deps_siv

! pack internal variables
        if (self%vari) then
            vip = 0.d0
            vip(1) = kapl
            vip(2:1+self%ndimsi) = ep/rac2
            vip(8) = kavs
            vip(9:8+self%ndimsi) = ev/rac2
        end if

999     continue
    end subroutine Integrate

! =====================================================================
!  PLASTICITY COMPUTATION AND TANGENT OPERATORS
! =====================================================================

    subroutine ComputePlasticity(self, eps, ka, ep, sig, deps_sig)

        implicit none

        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8), intent(in)             :: eps(:)
        real(kind=8), intent(inout)          :: ep(:), ka
        real(kind=8), intent(out)            :: sig(:), deps_sig(:, :)
! --------------------------------------------------------------------------------------------------
! eps       deformation a la fin du pas de temps
! ka        deformation plastique cumulee (in=debut pas de temps, out=fin)
! ep        deformation plastique (in=debut, out=fin)
! sig       ontrainte en fin de pas de temps (resi) ou au debut du pas de temps (not resi)
! deps_sig  derivee dsig / deps
! --------------------------------------------------------------------------------------------------
        integer(kind=8), parameter:: ELAS = 0
        integer(kind=8), parameter:: PLAS = 1
! --------------------------------------------------------------------------------------------------
        real(kind=8)    :: krn(self%ndimsi), pdev(self%ndimsi, self%ndimsi)
        real(kind=8)    :: sige(self%ndimsi), ne(self%ndimsi), neq
        integer(kind=8)         :: state
        real(kind=8)    :: n(self%ndimsi), dep(self%ndimsi)
        integer(kind=8)         :: matr
        real(kind=8)    :: coef
! --------------------------------------------------------------------------------------------------

!   INITIALISATION

        ! Kronecker et projecteur deviatorique
        krn = kron(self%ndimsi)/sqrt(3.d0)
        pdev = identity(self%ndimsi)-proten(krn, krn)

        ! Contrainte elastique
        sige = self%mat%troiskp*sph_norm(eps)*krn+self%mat%deuxmup*deviator(eps-ep)
        ne = deviator(sige)-self%mat%prag*ep
        neq = sqrt(1.5d0*dot_product(ne, ne))
! --------------------------------------------------------------------------------------------------

!  CALCUL DES CONTRAINTES

        ! Elasticite
        if (neq .le. self%mat%sc) then
            state = ELAS
            sig = sige

            ! Plasticite
        else
            state = PLAS
            n = self%mat%sc*ne/neq
            dep = (ne-n)/(self%mat%deuxmup+self%mat%prag)
            sig = sige-self%mat%deuxmup*dep
            ep = ep+dep
            ka = ka+sqrt(2.d0/3.d0*dot_product(dep, dep))
        end if

! --------------------------------------------------------------------------------------------------

!  MATRICES TANGENTES

        if (.not. self%rigi) goto 999

        ! Selection de la matrice (elastique ou plastique)
        if (self%elas) then
            matr = ELAS
        else if (self%pred) then
            matr = merge(ELAS, PLAS, neq .le. (1-1.d-6)*self%mat%sc)
        else
            matr = state
        end if

        ! Matrice elastique
        if (matr .eq. ELAS) then
            deps_sig = self%mat%troiskp*proten(krn, krn)+self%mat%deuxmup*pdev
        end if

        ! Matrice plastique
        if (matr .eq. PLAS) then
            coef = self%mat%deuxmup/(self%mat%deuxmup+self%mat%prag)
            deps_sig = self%mat%troiskp*proten(krn, krn) &
                       +coef*(self%mat%sc/neq*self%mat%deuxmup+self%mat%prag)*pdev &
                       -1.5d0*self%mat%deuxmup*self%mat%sc/neq*coef*proten(ne/neq, ne/neq)
        end if

999     continue
    end subroutine ComputePlasticity

! =====================================================================
!  VISCOELASTICITY COMPUTATION AND TANGENT OPERATOR
! =====================================================================

    subroutine ComputeViscoElasticity(self, eps, ka, ev, sig, deps_sig)

        implicit none

        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8), intent(in)             :: eps(:)
        real(kind=8), intent(inout)          :: ev(:), ka
        real(kind=8), intent(out)            :: sig(:), deps_sig(:, :)
! --------------------------------------------------------------------------------------------------
! eps       deformation a la fin du pas de temps
! ka        deformation visqueuse cumulee (in=debut pas de temps, out=fin)
! ev        deformation viscoelastique (in=debut, out=fin)
! sig       ontrainte en fin de pas de temps (resi) ou au debut du pas de temps (not resi)
! deps_sig  derivee dsig / deps
! --------------------------------------------------------------------------------------------------
        integer(kind=8), parameter:: SING = -1
        integer(kind=8), parameter:: REGU = 0
        integer(kind=8), parameter:: LINE = 1
        integer(kind=8), parameter:: ELAS = 2
        integer(kind=8), parameter:: NONLIN = 3
! --------------------------------------------------------------------------------------------------
        real(kind=8)      :: krn(self%ndimsi), pdev(self%ndimsi, self%ndimsi)
        real(kind=8)      :: v_h, v_d, c_h, c_d
        real(kind=8)      :: sigi(self%ndimsi), sigi_h, sigi_d(self%ndimsi), nori
        integer(kind=8)           :: state
        real(kind=8)      :: sign_h, sign_d(self%ndimsi)
        real(kind=8)      :: a1, a2, b0, b1, b2
        real(kind=8)      :: dev(self%ndimsi)
        real(kind=8)      :: xmin, xgmin, xmax, x, xg, equ, xm, xp, xgm, xgp, equm, equp
        real(kind=8)      :: equmin, equmax, dx_xg, dxg_equ, dx_equ
        integer(kind=8)           :: ite
        type(newton_state):: mem
        integer(kind=8)           :: matr
        real(kind=8)      :: c_v_sign(self%ndimsi), h(self%ndimsi)
        real(kind=8)      :: ca1(self%ndimsi), ca2(self%ndimsi)
        real(kind=8)      :: cb0(self%ndimsi), cb1(self%ndimsi), cb2(self%ndimsi)
        real(kind=8)      :: coef0, coef1(self%ndimsi), cx(self%ndimsi), num, den
! --------------------------------------------------------------------------------------------------

!   Initialisation
        krn = kron(self%ndimsi)/sqrt(3.d0)
        pdev = identity(self%ndimsi)-proten(krn, krn)

        v_h = 1-2*self%mat%nud
        v_d = 1+self%mat%nud
        c_h = self%mat%troiskv
        c_d = self%mat%deuxmuv

!   Contrainte instantanee
        sigi_h = c_h*sph_norm(eps-ev)
        sigi_d = c_d*deviator(eps-ev)
        sigi = sigi_h*krn+sigi_d
        nori = sqrt(v_h*sigi_h**2+v_d*dot_product(sigi_d, sigi_d))

!   Cas a traiter parmi lineaire, non lineaire singulier et non lineaire regulier
        if (self%visc_lin) then
            state = LINE
        else if (sqrt(dot_product(sigi, sigi)) .le. self%small) then
            state = SING
        else
            state = REGU
        end if

        ! Cas lineaire
        if (state .eq. LINE) then
            b0 = self%mat%vsc
            b1 = c_h*v_h*b0
            b2 = c_d*v_d*b0
            sig = sigi_h/(1+b1)*krn+sigi_d/(1+b2)
            dev = b0*(v_h*sph_norm(sig)*krn+v_d*deviator(sig))
            ev = ev+dev
        end if

        ! Cas non lineaire singulier
        if (state .eq. SING) then
            sig = 0
            dev = 0
        end if

!   Cas non lineaire regulier
        if (state .eq. REGU) then

            ! Coefficients du probleme scalaire
            sign_h = sigi_h/nori
            sign_d = sigi_d/nori

            a1 = v_h*sign_h**2
            a2 = v_d*dot_product(sign_d, sign_d)
            b0 = self%mat%vsc*nori**(self%mat%ga-1)
            b1 = c_h*v_h*b0
            b2 = c_d*v_d*b0

            !   Resolution de l'equation scalaire
            xmin = sqrt(a1/(1+b1)**2+a2/(1+b2)**2)
            xmax = 1.d0

            ! On teste si la borne x=1 est solution
            if (abs(xmax-xmin) .le. self%cvuser) then
                x = 0.5d0*(xmin+xmax)
                goto 100
            end if
            xm = 1.d0-self%cvuser
            xgm = xm**(self%mat%ga-1)
            equm = a1/(1+b1*xgm)**2+a2/(1+b2*xgm)**2-xm**2
            if (equm .ge. 0) then
                x = (xm+1.d0)*0.5d0
                goto 100
            end if

            ! Conditions de robustesse numérique
            xgmin = xmin**(self%mat%ga-1)
            equmax = a1/(1+b1)**2+a2/(1+b2)**2-1.d0
            equmin = a1/(1+b1*xgmin)**2+a2/(1+b2*xgmin)**2-xmin**2
            ASSERT(xmax-xmin > self%cvuser)
            ASSERT(equmax .lt. 0.d0)
            ASSERT(equmin .gt. 0.d0)

            x = merge(xmin, xmax, self%mat%ga .lt. 2)
            do ite = 1, self%itemax
                xg = x**(self%mat%ga-1)
                equ = a1/(1+b1*xg)**2+a2/(1+b2*xg)**2-x**2

                ! Test de convergence
                if (equ .ge. 0) then
                    xm = x
                    equm = equ
                    xp = min(x+self%cvuser, xmax)
                    xgp = xp**(self%mat%ga-1)
                    equp = a1/(1+b1*xgp)**2+a2/(1+b2*xgp)**2-xp**2
                else
                    xp = x
                    equp = equ
                    xm = max(x-self%cvuser, xmin)
                    xgm = xm**(self%mat%ga-1)
                    equm = a1/(1+b1*xgm)**2+a2/(1+b2*xgm)**2-xm**2
                end if

                ! Convergence
                if (equm*equp .le. 0) then
                    ! Estimation robuste par une corde
                    num = (xm*equp-xp*equm)
                    den = equp-equm
                    if (abs(den) .le. abs(num)*self%small) then
                        x = (xm+xp)*0.5d0
                    else
                        x = num/den
                    end if
                    goto 100
                end if

                ! Calcul de la derivee
                dx_xg = (self%mat%ga-1)*x**(self%mat%ga-2)
                dxg_equ = -2*a1*b1/(1+b1*xg)**3-2*a2*b2/(1+b2*xg)**3
                dx_equ = dx_xg*dxg_equ-2*x

                ! Update
                x = utnewt(x, -equ, -dx_equ, ite, mem, xmin=xmin, xmax=xmax)
            end do
            self%exception = 1
            goto 999
100         continue

            xg = x**(self%mat%ga-1)
            sig = sigi_h/(1+b1*xg)*krn+sigi_d/(1+b2*xg)
            dev = b0*xg*(v_h*sph_norm(sig)*krn+v_d*deviator(sig))
            ev = ev+dev
        end if

        ! Deformation visqueuse cumulee
        ka = ka+sqrt(2.d0/3.d0*dot_product(dev, dev))

!   Calcul de la matrice tangente
        if (.not. self%rigi) goto 999

        ! Selection de la matrice (elastique ou plastique)
        if (self%elas .or. state .eq. SING) then
            matr = ELAS
        else
            matr = merge(LINE, NONLIN, state .eq. LINE)
        end if

        ! Cas lineaire
        if (matr .eq. LINE) then
            deps_sig = c_h/(1+b1)*proten(krn, krn)+c_d/(1+b2)*pdev
        end if

        ! Cas non lineaire singulier ou elastique
        if (matr .eq. ELAS) then
            deps_sig = c_h*proten(krn, krn)+c_d*pdev
        end if

        ! Cas non lineaire regulier
        if (matr .eq. NONLIN) then
            c_v_sign = c_h*v_h*sign_h*krn+c_d*v_d*sign_d

            h = c_h*v_h/(1+b1*xg)**2*sign_h*krn+c_d*v_d/(1+b2*xg)**2*sign_d

            ca1 = 2*v_h*(c_h*sign_h*krn-sign_h**2*c_v_sign)
            ca2 = 2*v_d*(c_d*sign_d-dot_product(sign_d, sign_d)*c_v_sign)
            cb0 = (self%mat%ga-1)*self%mat%vsc*nori**(self%mat%ga-1)*c_v_sign
            cb1 = c_h*v_h*cb0
            cb2 = c_d*v_d*cb0

            coef0 = (2*a1*(1+b1*self%mat%ga*xg)/(1+b1*xg)**3 &
                     +2*a2*(1+b2*self%mat%ga*xg)/(1+b2*xg)**3)
            coef1 = ((1+b1*xg)*ca1-2*a1*xg*cb1)/(1+b1*xg)**3 &
                    +((1+b2*xg)*ca2-2*a2*xg*cb2)/(1+b2*xg)**3
            cx = coef1/coef0

            deps_sig = c_h/(1+b1*xg)*proten(krn, krn)+c_d/(1+b2*xg)*pdev &
                       -xg*proten(h, cb0) &
                       -(self%mat%ga-1)*xg*b0*proten(h, cx)
        end if

999     continue
    end subroutine ComputeViscoElasticity

end module kichenin_nl_module
