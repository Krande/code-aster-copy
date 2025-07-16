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

function lcdp_compute(resi, rigi, elas, itemax, prec, m, eps, ep, ka, state, &
                      s, deps_s, vip) result(iret)

    use lcdp_module, only: dp_material

    implicit none
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/proten.h"
#include "asterfort/zerop2.h"
#include "asterfort/zerofr.h"
#include "asterfort/utmess.h"

    aster_logical, intent(in)     :: resi, rigi, elas
    type(dp_material), intent(in) :: m
    integer(kind=8), intent(in)           :: itemax
    real(kind=8), intent(in)      :: prec
    real(kind=8), intent(in)      :: eps(:)
    real(kind=8), intent(inout)   :: ep(:), ka
    integer(kind=8), intent(inout)        :: state
    real(kind=8), intent(out)     :: s(:), deps_s(:, :)
    integer(kind=8)                      :: iret
    real(kind=8) :: vip(9)
! ----------------------------------------------------------------------
!             LOI DRUCK_PRAG_N_A
! ----------------------------------------------------------------------
    real(kind=8), dimension(6), parameter:: kron = [1.d0, 1.d0, 1.d0, 0.d0, 0.d0, 0.d0]
! ----------------------------------------------------------------------
    integer(kind=8)     :: i, nrac
    real(kind=8):: kr(size(eps)), id(size(eps), size(eps)), prodev(size(eps), size(eps))
    real(kind=8):: sel(size(eps)), selh, selq, seld(size(eps))
    real(kind=8):: kam, fm, kas, fs, bs, bk, rks, fu, lbd
    real(kind=8):: deph, depd(size(eps)), depq, q
    real(kind=8):: p0, p1, p2, sol(2), coef, x0
    real(kind=8):: dselq_ka, dselh_ka, dka_blbd, dselh_sh, dselq_sh, dselh_q, dselq_q
    real(kind=8):: deps_selh(size(eps)), deps_selq(size(eps)), deps_seld(size(eps), size(eps))
    real(kind=8):: deps_sh(size(eps)), deps_q(size(eps))
    real(kind=8):: dka_s(size(eps)), deps_ka(size(eps))

    integer(kind=8) :: iter
! ----------------------------------------------------------------------

!   Initialisation
    iret = 0
    kr = kron(1:size(eps))
    kam = ka
    id = 0
    forall (i=1:size(eps)) id(i, i) = 1
    prodev = id-proten(kr, kr)/3.0

!   Contrainte elastique
    sel = m%lambda*sum(eps(1:3)-ep(1:3))*kr+m%deuxmu*(eps-ep)
    selh = sum(sel(1:3))/3.d0
    seld = sel-selh*kr
    selq = sqrt(1.5d0*dot_product(seld, seld))

! ======================================================================
!               INTEGRATION DE LA LOI DE COMPORTEMENT
! ======================================================================

    if (.not. resi) goto 500

! ----------------------------------------------------------------------
!  REGIME ELASTIQUE
! ----------------------------------------------------------------------

!   Borne du critere elastique
    fm = f_fhat(kam)
    if (fm .le. 0) then
        state = 0
        s = sel
        goto 500
    end if

! ----------------------------------------------------------------------
!  REGIME PLASTIQUE SINGULIER
! ----------------------------------------------------------------------

    depd = seld/m%deuxmu
    depq = selq/m%deuxmu
    kas = kam+2.0/3.0*depq
    bs = f_b(kas)

!   Test d'appartenance au cone
    fs = f_fhat(kas)
    if (fs .ge. 0) then
        rks = f_ecro(kas)
        deph = (selh-rks/m%troisa)/m%troisk
        ka = kas
        ep = ep+deph*kr+depd
        s = rks/m%troisa*kr
        state = 2
        goto 500
    end if

! ----------------------------------------------------------------------
!  REGIME PLASTIQUE REGULIER
! ----------------------------------------------------------------------

!   Resolution f_fhat(ka)=0 avec kam < ka < kas et f_fhat decroissante dans l'intervalle

!   Cas de la solution saturee
    fu = f_fhat(m%kau)

!   On cherche le plus petit zero de f_fhat
!   f_fhat est quadratique ou lineaire si l'ecrouissage est lineaire ou parabolique pour R et b
    if ((m%type_dp .eq. 1) .or. (m%type_dp .eq. 2)) then
        if (kam .ge. m%kau .or. fu .ge. 0) then
            ka = m%kau+fu/m%troismu
            lbd = ka-kam
            ! Cas non sature
        else
            !Ecrouissage lineaire
            if (m%type_dp .eq. 1) then
                p2 = m%troisk*m%troisa*m%b0/m%kau
                !Ecrouissage parabolique
            else if (m%type_dp .eq. 2) then
                p2 = m%troisk*m%troisa*m%b0/m%kau+1/m%kau*(m%syultm-m%sy)/m%kau
            end if
            p0 = fm/p2
            p1 = dka_fhat(kam)/p2
            call zerop2(p1, p0, sol, nrac)
            lbd = sol(nrac)
            ka = kam+lbd
        end if
!   Dans les autres cas, utilisation d'une methode de Brent pour trouver le zero de f_fhat
    else
        if (kam .ne. 0.0d0) then
            x0 = kam+0.5*kam
            call zerofr(2, 'AUTO', f_fhat, kam, x0, prec, abs(itemax), ka, iret, iter)
        else
            x0 = 1.0d0
            call zerofr(2, 'AUTO', f_fhat, kam, x0, prec, abs(itemax), ka, iret, iter)
        end if
        lbd = ka-kam
    end if

    bk = f_b(ka)
    q = m%troismu*lbd/selq
    ep = ep+lbd*(bk*kr+1.5d0*seld/selq)
    s = sel-lbd*m%troisk*bk*kr-q*seld
    state = 1

!  Condition pour garantir l'unicite de la solution au probleme d'integration locale du comportement
    if (-3.0d0*m%troisk*(dka_b(ka)*m%a*lbd+f_b(ka)*m%a) .ge. (dka_ecro(ka)+m%troismu)) then
        call utmess('F', 'COMPOR6_15')
    end if

! ======================================================================
!                           MATRICES TANGENTES
! ======================================================================

500 continue
    if (.not. rigi) goto 999

    ! Initialisation
    lbd = merge(ka-kam, 0.d0, resi)

    ! Regime elastique ou regime singulier en phase de prediction (non derivable)
    if (elas .or. state .eq. 0 .or. (.not. resi .and. state .eq. 2)) then

        deps_s = m%lambda*proten(kr, kr)+m%deuxmu*id

        ! Regime plastique regulier
    else if (state .eq. 1) then

        ! Initialisation
        q = m%troismu*lbd/selq

        ! Variation de kappa (fonction implicite)
        coef = dka_fhat(ka)
        dselq_ka = -dselq_fhat(ka)/coef
        dselh_ka = -dselh_fhat(ka)/coef

        ! Variation des invariants de contrainte
        dka_blbd = dka_b(ka)*lbd+f_b(ka)
        dselh_sh = 1-m%troisk*dselh_ka*dka_blbd
        dselq_sh = -m%troisk*dselq_ka*dka_blbd
        dselh_q = m%troismu*dselh_ka/selq
        dselq_q = m%troismu*(dselq_ka/selq-lbd/selq**2)

        ! Variation des invariants elastiques
        deps_selh = m%troisk/3.0*kr
        deps_selq = m%troismu*seld/selq

        ! Variation de la contrainte
        deps_seld = m%deuxmu*prodev
        deps_sh = deps_selh*dselh_sh+deps_selq*dselq_sh
        deps_q = deps_selh*dselh_q+deps_selq*dselq_q
        deps_s = proten(kr, deps_sh)+(1-q)*deps_seld-proten(seld, deps_q)

        ! Regime singulier (hors phase de prediction sinon pas derivable)
    else if (state .eq. 2 .and. resi) then
        depd = seld/m%deuxmu
        dka_s = dka_ecro(ka)/m%troisa*kr
        deps_ka = 2.d0/3.d0/lbd*depd
        deps_s = proten(dka_s, deps_ka)

    else
        ASSERT(.false.)
    end if

!    Fin de la routine
999 continue

contains

! ----------------------------------------------------------------------------------------
!  Liste des fonctions intermediaires et leurs derivees
! ----------------------------------------------------------------------------------------

! Variables globales partagees
!     type(dp_material):: m
!     real(kind=8)      :: selh,selq, kam
! --> Derivees: dselh_*, dselq_*

    real(kind=8) function f_b(ka)
        real(kind=8)::ka
        !Ecrouissage exponentiel pour b
        if (m%type_dp .eq. 3) then
            f_b = (m%b0-m%bultm)*exp(-ka/m%kac)+m%bultm
            !Ecrouissage lineaire pour b
        else
            f_b = m%b0*max(0.d0, 1-ka/m%kau)
        end if
    end function f_b

    real(kind=8) function dka_b(ka)
        real(kind=8)::ka
        !Ecrouissage exponentiel pour b
        if (m%type_dp .eq. 3) then
            dka_b = -(m%b0-m%bultm)/m%kac*exp(-ka/m%kac)
            !Ecrouissage lineaire pour b
        else
            dka_b = m%b0*merge(-1/m%kau, 0.d0, ka .le. m%kau)
        end if
    end function dka_b

    real(kind=8) function f_ecro(ka)
        real(kind=8)::ka
        real(kind=8):: P2, P1
        !Ecrouissage lineaire pour R
        if (m%type_dp .eq. 1) then
            f_ecro = m%sy+m%h*min(ka, m%kau)
            !Ecrouissage parabolique pour R
        else if (m%type_dp .eq. 2) then
            if (ka .le. m%kau) then
                P2 = -1/m%kau*(m%syultm-m%sy)/m%kau
                P1 = 2*(m%syultm-m%sy)/m%kau
                f_ecro = P2*ka**2+P1*ka+m%sy
            else
                f_ecro = m%syultm
            end if
            !Ecrouissage exponentiel pour R
        else if (m%type_dp .eq. 3) then
            f_ecro = (m%sy-m%syultm)*exp(-ka/m%kac)+m%syultm
        else
            ASSERT(.false.)
        end if
    end function f_ecro

    real(kind=8) function dka_ecro(ka)
        real(kind=8)::ka
        real(kind=8):: P2, P1
        !Ecrouissage lineaire pour R
        if (m%type_dp .eq. 1) then
            dka_ecro = m%h*merge(1.d0, 0.d0, ka .le. m%kau)
            !Ecrouissage parabolique pour R
        else if (m%type_dp .eq. 2) then
            if (ka .le. m%kau) then
                P2 = -1/m%kau*(m%syultm-m%sy)/m%kau
                P1 = 2*(m%syultm-m%sy)/m%kau
                dka_ecro = 2*P2*ka+P1
            else
                dka_ecro = 0
            end if
            !Ecrouissage exponentiel pour R
        else if (m%type_dp .eq. 3) then
            dka_ecro = -(m%sy-m%syultm)/m%kac*exp(-ka/m%kac)
        else
            ASSERT(.false.)
        end if
    end function dka_ecro

    real(kind=8) function f_fhat(ka)
        real(kind=8)::ka
        real(kind=8)::lbd
        lbd = ka-kam
        f_fhat = selq+m%troisa*selh-f_ecro(ka)-m%troismu*lbd-m%troisa*m%troisk*f_b(ka)*lbd
    end function f_fhat

    real(kind=8) function dka_fhat(ka)
        real(kind=8)::ka
        real(kind=8)::lbd, dka_blbd
        lbd = ka-kam
        dka_blbd = dka_b(ka)*lbd+f_b(ka)
        dka_fhat = -dka_ecro(ka)-m%troismu-m%troisa*m%troisk*dka_blbd
    end function dka_fhat

    real(kind=8) function dselq_fhat(ka)
        real(kind=8)::ka
        dselq_fhat = 1
    end function dselq_fhat

    real(kind=8) function dselh_fhat(ka)
        real(kind=8)::ka
        dselh_fhat = m%troisa
    end function dselh_fhat

end function
