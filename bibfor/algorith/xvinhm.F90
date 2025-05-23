! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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
! person_in_charge: daniele.colombo at ifpen.fr
!
subroutine xvinhm(ds_thm, jmate, ndim, &
                  cohes, dpf, saut, sautm, nd, lamb, &
                  w11m, rho11m, alpha, job, pf, &
                  rho11, w11, ipgf, rela, dsidep, &
                  delta, r, am)
!
    use THM_type
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterc/r8prem.h"
#include "asterfort/lceitc.h"
#include "asterfort/lceiou.h"
!
! ======================================================================
!
! ROUTINE MODELE HM-XFEM (CAS DE LA FRACTURE)
!
! CALCUL DES VARIABLES INTERNES (MECANIQUE ET HYDRAULIQUE)
!
! ----------------------------------------------------------------------
!
    type(THM_DS), intent(inout) :: ds_thm
    integer :: jmate, ndim, i, ipgf
    real(kind=8) :: cliq, vim(2), vip(2), cohes(5), rho11, rho11m
    real(kind=8) :: dsidep(6, 6), delta(6), eps, vim2(9), vip2(9), rela
    real(kind=8) :: w11, w11m, varbio, dpf, psp, psm, saut(3), lamb(3)
    real(kind=8) :: sautm(3), alpha(5), nd(3), rho110, r, pf, am(3)
    character(len=8)  :: job
    character(len=16) :: option
!
    vim(:) = 0.d0
    vip(:) = 0.d0
    vim2(:) = 0.d0
    vip2(:) = 0.d0
    dsidep(:, :) = 0.d0
    delta(:) = 0.d0
!
! - Get material parameters
!
    rho110 = ds_thm%ds_material%liquid%rho
    cliq = ds_thm%ds_material%liquid%unsurk
!
!   INITIALISATION DE LA VARIABLE INTERNE
!
    vim2(4) = cohes(1)
    vim2(2) = cohes(2)
    vim(1) = cohes(4)
    vim(2) = cohes(5)
!
    rho11 = vim(1)+rho110
    rho11m = vim(1)+rho110
!
    w11 = vim(2)
    w11m = vim(2)
!
!   PREDICTION: COHES(3)=1 ; CORRECTION: COHES(3)=2
!
    if (nint(cohes(3)) .eq. 1) then
        option = 'RIGI_MECA_TANG'
    else if (nint(cohes(3)) .eq. 2) then
        option = 'FULL_MECA'
    else
        option = 'FULL_MECA'
    end if
    if (job .eq. 'MATRICE' .and. option .eq. 'FULL_MECA') then
        eps = 100.*r8prem()
        vim2(4) = min(1.d0, vim2(4)*(1+eps))
    end if
!
!   UTILISATION DE LA LOI COHESIVE MIXTE TALON-CURNIER
    if (nint(rela) .eq. 3) then
        call lceitc('RIGI', ipgf, 1, jmate, option, &
                    lamb, am, delta, dsidep, vim2, &
                    vip2, r, pfluide=pf)
    else if (nint(rela) .eq. 4) then
        call lceiou('RIGI', ipgf, 1, jmate, option, &
                    lamb, am, delta, dsidep, vim2, &
                    vip2, r, pfluide=pf)
    end if
!
    if (option .eq. 'FULL_MECA') then
!   CALCUL DE LA VARIABLE INTERNE : MASSE VOLUMIQUE
        varbio = dpf*cliq
        if (varbio .gt. 5.d0) then
            ASSERT(.false.)
        end if
!
        vip(1) = -rho110+(vim(1)+rho110)*exp(varbio)
        rho11 = vip(1)+rho110
        rho11m = vim(1)+rho110
!
!   CALCUL DE LA VARIABLE INTERNE : APPORTS MASSIQUES
!   (SEULEMENT UTILE POUR LE CAS DU SECOND-MEMBRE)
!
        psp = 0.d0
        psm = 0.d0
        do i = 1, ndim
            psp = psp-saut(i)*nd(i)
            psm = psm-sautm(i)*nd(i)
        end do
!
        vip(2) = rho11*psp-rho11m*psm
        w11 = vip(2)+w11m
    end if
!
    alpha(1) = vip2(4)
    alpha(2) = vip2(2)
    alpha(4) = vip(1)
    alpha(5) = w11
!
    if (job .eq. 'ACTU_VI') then
        alpha(3) = 1.d0
    else if (job .eq. 'MATRICE') then
        alpha(3) = 2.d0
    end if
end subroutine
