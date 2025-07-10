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
subroutine pipeab(mat, dtau, vim, sup, sud, mup, mud, nsol, sol, sgn)

    implicit none
#include "asterfort/assert.h"
#include "asterc/r8gaem.h"
#include "asterfort/rcvalb.h"
#include "asterfort/zerop2.h"

    integer(kind=8)                  :: mat
    real(kind=8), intent(in) :: dtau, vim(:), sup(:), sud(:), mup(:), mud(:)
    integer(kind=8), intent(out)     :: nsol, sgn(2)
    real(kind=8), intent(out):: sol(2)

!
! ----------------------------------------------------------------------
!     pilotage pred_elas pour la loi d'interface czm_lab_mix
! ----------------------------------------------------------------------
! in  mat    : materiau
! in  dtau   : 2nd membre de l'equation f(eta)=dtau
! in  vim    : variables internes en t-
! in  sig_f  : partie cinematique correspondant a fixe_cste
! in  sig_p  : partie cinematique correspondant a fixe_pilo
! out nsol   : nombre de solutions (-1: le point ne contribue pas au pilotage)
! out sol    : solutions de l'equation de pilotage
! out sgn    : signe de la pente de la fonction de pilotage en chaque solution
! ----------------------------------------------------------------------
    integer(kind=8) :: cod(6), cine, nrac
    character(len=16) :: nom(6)
    real(kind=8) :: sc, dc, alpha, beta, s0, d0, r, ka, sr, val(6)
    real(kind=8) :: rac(2), p0, p1, p2
    real(kind=8), dimension(size(sup)) :: sig_f, sig_p, s_f, s_p
    integer(kind=8), parameter:: unilater = 0, glis_1d = 1, glis_2d = 2
! --------------------------------------------------------------------------------------------------
    data nom/'SIGM_C', 'GLIS_C', 'ALPHA', 'BETA', 'PENA_LAGR', 'CINEMATIQUE'/
! --------------------------------------------------------------------------------------------------

! --------------------------------------------------------------------------------------------------
!  Initialisation
! --------------------------------------------------------------------------------------------------

    ! Parametres du materiau
    call rcvalb('FPG1', 1, 1, '+', mat, &
                ' ', 'CZM_LAB_MIX', 0, ' ', [0.d0], &
                6, nom, val, cod, 2)

    sc = val(1)
    dc = val(2)
    alpha = val(3)
    beta = val(4)
    d0 = dc*beta/alpha
    s0 = sc*(alpha+beta)**(alpha+beta)/(alpha**alpha*beta**beta)
    r = val(5)*sc/dc
    cine = nint(val(6))

    ! Endommagement cible et seuil correspondant
    ka = vim(1)+dc*dtau
    sr = s0*(ka/d0)**alpha/(ka/d0+1.d0)**(alpha+beta)+r*ka

    ! Terme de pilotage
    sig_f = mup+r*sup
    sig_p = mud+r*sud

! --------------------------------------------------------------------------------------------------
!  Calcul des solutions en cinematique de glissement (1D ou 2D)
! --------------------------------------------------------------------------------------------------

    if (cine .eq. glis_1d .or. cine .eq. glis_2d) then

        ! termes de pilotage
        s_f = 0
        s_p = 0
        if (cine .eq. glis_1d) then
            s_f(2) = sig_f(2)
            s_p(2) = sig_p(2)
        else if (cine .eq. glis_2d) then
            s_f(2:) = sig_f(2:)
            s_p(2:) = sig_p(2:)
        end if

        ! Calcul des racines
        p2 = dot_product(s_p, s_p)
        p1 = 2*dot_product(s_p, s_f)
        p0 = dot_product(s_f, s_f)-sr**2

        if (abs(p2) .lt. abs(p0)/r8gaem()) then
            ! Point non pilotable ou ne contribuant pas au pilotage
            nsol = merge(-1, 0, p0 .le. 0)
            goto 999
        end if
        call zerop2(p1/p2, p0/p2, sol, nsol)

        ! Organisation des solutions dans le cas simple du glissement
        nsol = merge(2, 0, nsol .eq. 2)
        if (nsol .eq. 2) then
            sgn(1) = 1
            sgn(2) = -1
        end if
    end if

! --------------------------------------------------------------------------------------------------
!  Calcul des solutions en cinematique de glissement (1D ou 2D)
! --------------------------------------------------------------------------------------------------

    if (cine .eq. unilater) then

        s_f = sig_f
        s_p = sig_p

        ! Solutions correspondant a l'activation du terme sur la normale
        p2 = dot_product(s_p, s_p)
        p1 = 2*dot_product(s_p, s_f)
        p0 = dot_product(s_f, s_f)-sr**2

        ! Point non pilotable ou ne contribuant pas au pilotage
        if (abs(p2) .lt. abs(p0)/r8gaem()) then
            s_f(1) = max(0.d0, s_f(1))
            nsol = merge(-1, 0, dot_product(s_f, s_f)-sr**2 .lt. 0.d0)
            goto 999
        end if
        call zerop2(p1/p2, p0/p2, rac, nrac)

        nsol = 0
        if (nrac .eq. 2) then
            if (sig_f(1)+rac(2)*sig_p(1) .ge. 0.d0) then
                nsol = nsol+1
                sol(nsol) = rac(2)
                sgn(nsol) = -1
            end if
            if (sig_f(1)+rac(1)*sig_p(1) .ge. 0.d0) then
                nsol = nsol+1
                sol(nsol) = rac(1)
                sgn(nsol) = 1
            end if
        end if

        ! Solutions correspondant a l'annulation du terme sur la normale
        s_f(1) = 0
        s_p(1) = 0

        p2 = dot_product(s_p, s_p)
        p1 = 2*dot_product(s_p, s_f)
        p0 = dot_product(s_f, s_f)-sr**2

        if (abs(p2) .ge. abs(p0)/r8gaem()) then
            call zerop2(p1/p2, p0/p2, rac, nrac)

            if (nrac .eq. 2) then
                if (sig_f(1)+rac(2)*sig_p(1) .lt. 0.d0) then
                    nsol = nsol+1
                    ASSERT(nsol .le. 2)
                    sol(nsol) = rac(2)
                    sgn(nsol) = -1
                end if
                if (sig_f(1)+rac(1)*sig_p(1) .lt. 0.d0) then
                    nsol = nsol+1
                    ASSERT(nsol .le. 2)
                    sol(nsol) = rac(1)
                    sgn(nsol) = 1
                end if
            end if
        end if

    end if

999 continue
end subroutine
