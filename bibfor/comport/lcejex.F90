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
subroutine lcejex(fami, kpg, ksp, ndim, mate, &
                  option, am, da, sigma, dsidep, &
                  vim, vip)
!
!
    implicit none
#include "asterf_types.h"
#include "asterfort/r8inir.h"
#include "asterfort/rcvalb.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
    integer(kind=8) :: mate, ndim, kpg, ksp
    real(kind=8) :: am(ndim), da(ndim), sigma(6), dsidep(6, 6)
    real(kind=8) :: vim(*), vip(*)
    character(len=16) :: option
    character(len=*) :: fami
!
!-----------------------------------------------------------------------
!  LOI DE COMPORTEMENT COHESIVE : CZM_EXP_REG
!  POUR LES ELEMENTS DE JOINT 2D ET 3D
!
! IN : AM SAUT INSTANT - : AM(1) = SAUT NORMAL, AM(2) = SAUT TANGENTIEL
! IN : DA    INCREMENT DE SAUT
! IN : MATE, OPTION, VIM
! OUT : SIGMA , DSIDEP , VIP
!-----------------------------------------------------------------------
!
    aster_logical :: resi, rigi, elas
    integer(kind=8) :: diss, i, j
    real(kind=8) :: sc, gc, lc, k0, val(4), rtan
    real(kind=8) :: a(ndim), na, ka, kap, r0, rc, beta, rk, ra, coef, coef2
    integer(kind=8) :: cod(4)
    character(len=16) :: nom(4)
    character(len=1) :: poum
    blas_int :: b_incx, b_incy, b_n
!
! OPTION CALCUL DU RESIDU OU CALCUL DE LA MATRICE TANGENTE
!
    resi = option(1:9) .eq. 'FULL_MECA' .or. option .eq. 'RAPH_MECA'
    rigi = option(1:9) .eq. 'FULL_MECA' .or. option(1:9) .eq. 'RIGI_MECA'
    elas = option .eq. 'FULL_MECA_ELAS' .or. option .eq. 'RIGI_MECA_ELAS'
!
!
! CALCUL DU SAUT EN T+
!
    b_n = to_blas_int(ndim)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, am, b_incx, a, b_incy)
    if (resi) then
        b_n = to_blas_int(ndim)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, 1.d0, da, b_incx, a, &
                   b_incy)
    end if
!
!
! RECUPERATION DES PARAMETRES PHYSIQUES
!
    nom(1) = 'GC'
    nom(2) = 'SIGM_C'
    nom(3) = 'PENA_ADHERENCE'
    nom(4) = 'PENA_CONTACT'
!
    if (option .eq. 'RIGI_MECA_TANG') then
        poum = '-'
    else
        poum = '+'
    end if
!
    call rcvalb(fami, kpg, ksp, poum, mate, &
                ' ', 'RUPT_FRAG', 0, ' ', [0.d0], &
                4, nom, val, cod, 2)
!
    gc = val(1)
    sc = val(2)
    lc = gc/sc
    k0 = lc*val(3)
    r0 = sc*exp(-k0/lc)/k0
    beta = val(4)
!
! -- INITIALISATION
!
    ka = max(k0, vim(1))
    rtan = 0.d0
    do i = 2, ndim
        rtan = rtan+a(i)**2
    end do
    na = sqrt(max(0.d0, a(1))**2+rtan)
    rk = sc*exp(-ka/lc)/ka
    rc = rk+beta*(r0-rk)
!
!
! INITIALISATION COMPLEMENTAIRE POUR RIGI_MECA_TANG (SECANTE PENALISEE)
!
    if (.not. resi) then
        if (elas) then
            diss = 0
        else
            diss = nint(vim(2))
        end if
        goto 5000
    end if
!
!
! -- CALCUL DE LA CONTRAINTE
!
    call r8inir(6, 0.d0, sigma, 1)
!
!    CONTRAINTE DE CONTACT PENALISE
    sigma(1) = rc*min(0.d0, a(1))
!
!    CONTRAINTE DE FISSURATION
    if (na .le. ka) then
        diss = 0
        sigma(1) = sigma(1)+rk*max(0.d0, a(1))
        do i = 2, ndim
            sigma(i) = sigma(i)+rk*a(i)
        end do
    else
        diss = 1
        ra = sc*exp(-na/lc)/na
        sigma(1) = sigma(1)+ra*max(0.d0, a(1))
        do i = 2, ndim
            sigma(i) = sigma(i)+ra*a(i)
        end do
    end if
!
!
! -- ACTUALISATION DES VARIABLES INTERNES
!   V1 :  SEUIL, PLUS GRANDE NORME DU SAUT
!   V2 :  INDICATEUR DE DISSIPATION (0 : NON, 1 : OUI)
!   V3 :  INDICATEUR D'ENDOMMAGEMENT  (0 : SAIN, 1: ENDOM)
!   V4 :  POURCENTAGE D'ENERGIE DISSIPEE
!   V5 :  VALEUR DE L'ENERGIE DISSIPEE (V4*GC)
!   V6 :  VALEUR DE L'ENERGIE RESIDUELLE COURANTE
!   V7 A V9 : VALEURS DU SAUT
!
    kap = max(ka, na)
    vip(1) = kap
    vip(2) = diss
    if (na .le. k0) then
        vip(3) = 0.d0
    else
        vip(3) = 1.d0
    end if
    vip(4) = 1.d0-exp(-kap/lc)-0.5d0*kap*sc*exp(-kap/lc)/gc
    vip(5) = gc*vip(4)
    vip(6) = 0.5d0*(na**2)*sc*exp(-kap/lc)/kap
    vip(7) = a(1)
    vip(8) = a(2)
    if (ndim .eq. 3) then
        vip(9) = a(3)
    else
        vip(9) = 0.d0
    end if
!
! -- MATRICE TANGENTE
!
5000 continue
    if (.not. rigi) goto 999
!
    call r8inir(36, 0.d0, dsidep, 1)
!
!    MATRICE TANGENTE DE CONTACT PENALISE
    if (a(1) .le. 0.d0) dsidep(1, 1) = dsidep(1, 1)+rc
!
!    MATRICE TANGENTE DE FISSURATION
    if ((diss .eq. 0) .or. elas) then
!
        if (a(1) .gt. 0.d0) dsidep(1, 1) = dsidep(1, 1)+rk
        do i = 2, ndim
            dsidep(i, i) = dsidep(i, i)+rk
        end do
!
    else
!
        coef = (sc/lc+sc/na)/na**2
        coef2 = exp(-na/lc)
!
        if (a(1) .le. 0.d0) then
!
            do i = 2, ndim
                dsidep(i, i) = dsidep(i, i)+sc*coef2/na-coef*coef2*a(i)*a(i)
            end do
!
            if (ndim .eq. 3) then
                dsidep(2, 3) = dsidep(2, 3)-coef*coef2*a(2)*a(3)
                dsidep(3, 2) = dsidep(3, 2)-coef*coef2*a(3)*a(2)
            end if
!
        else
!
            do i = 1, ndim
                dsidep(i, i) = dsidep(i, i)+sc*coef2/na-coef*coef2*a(i)*a(i)
            end do
!
            do j = 1, ndim-1
                do i = j+1, ndim
                    dsidep(j, i) = dsidep(j, i)-coef*coef2*a(j)*a(i)
                    dsidep(i, j) = dsidep(i, j)-coef*coef2*a(i)*a(j)
                end do
            end do
!
        end if
!
    end if
!
999 continue
end subroutine
