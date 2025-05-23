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
! person_in_charge: daniele.colombo at ifpen.fr
! aslint: disable=W1504
!
subroutine xhmsa6(ds_thm, ndim, ipgf, imate, lamb, &
                  wsaut, nd, tau1, tau2, cohes, &
                  job, rela, alpha, dsidep, sigma, &
                  p, am, raug, wsautm, dpf, &
                  rho110)
!
    use THM_type
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/lcecli.h"
#include "asterfort/lcejex.h"
!
    type(THM_DS), intent(inout) :: ds_thm
    integer :: ndim, ipgf, imate
    real(kind=8) :: wsaut(3), lamb(3), am(3), dsidep(6, 6)
    real(kind=8) :: tau1(3), tau2(3), nd(3), wsautm(3)
    real(kind=8) :: alpha(5), p(3, 3), rho11, rho11m
    real(kind=8) :: cohes(5), rela, raug, dpf, w11, w11m
    character(len=8) :: job
! ======================================================================
!
!
! ROUTINE CONTACT (METHODE XFEM HPP - CALCUL ELEM.)
!
! --- CALCUL DU SAUT DE DEPLACEMENT EQUIVALENT [[UEG]]
!
! ----------------------------------------------------------------------
!
! IN  NDIM   : DIMENSION DE L'ESPACE
! IN  IPGF   : NUMÉRO DU POINTS DE GAUSS
! IN  IMATE  : ADRESSE DE LA SD MATERIAU
! IN  SAUT   : SAUT DE DEPLACEMENT
! IN  ND     : NORMALE À LA FACETTE ORIENTÉE DE ESCL -> MAIT
!                 AU POINT DE GAUSS
! IN  TAU1   : TANGENTE A LA FACETTE AU POINT DE GAUSS
! IN  TAU2   : TANGENTE A LA FACETTE AU POINT DE GAUSS
! IN  COHES  : VARIABLE INTERNE COHESIVE
! IN  JOB    : 'SAUT_EQ', 'MATRICE' OU 'VECTEUR'
! IN  RELA   : LOI COHESIVE 1:CZM_EXP_REG 2:CZM_LIN_REG
! OUT ALPHA  : SAUT DE DEPLACEMENT EQUIVALENT
! OUT DSIDEP : MATRICE TANGENTE DE CONTACT PENALISE ET DE FISSURATION
! OUT SIGMA  : CONTRAINTE
! OUT PP     : ND X ND
! OUT DNOR   : SAUT DEPLACEMENT NORMAL DANS LA BASE FIXE
! OUT DTANG  : SAUT DEPLACEMENT TANGENTIEL DANS LA BASE FIXE
! OUT P      : MATRICE DE PROJECTION SUR LE PLAN TANGENT
! OUT AM     : SAUT INSTANT - BASE LOCALE : AM(1) = SAUT NORMAL
!                                           AM(2) = SAUT TANGENTIEL
!
!
!
!
    integer :: i
    real(kind=8) :: vim(9), vip(9)
    real(kind=8) :: dsid2d(6, 6), dam(3)
    real(kind=8) :: sigma(6), cliq, varbio
    real(kind=8) :: rho110
    character(len=16) :: option
!
! ----------------------------------------------------------------------
!
    am(:) = 0.d0
    dam(:) = 0.d0
    p(:, :) = 0.d0
    dsidep(:, :) = 0.d0
    dsid2d(:, :) = 0.d0
    sigma(:) = 0.d0
    vim(:) = 0.d0
    vip(:) = 0.d0
!
! - Get material parameters
!
    rho110 = ds_thm%ds_material%liquid%rho
    cliq = ds_thm%ds_material%liquid%unsurk
!
    rho11 = cohes(4)+rho110
    rho11m = cohes(4)+rho110
!
    w11 = cohes(5)
    w11m = cohes(5)
!
! CONSTRUCTION DE LA MATRICE DE PASSAGE
!
! avec la nouvelle formulation, le saut est directement dans
! la bonne base
    do i = 1, ndim
        am(i) = wsaut(i)
    end do
!
    do i = 1, ndim
        p(1, i) = nd(i)
    end do
    do i = 1, ndim
        p(2, i) = tau1(i)
    end do
    if (ndim .eq. 3) then
        do i = 1, ndim
            p(3, i) = tau2(i)
        end do
    end if
!
! --- CALCUL DU SAUT DE DEPLACEMENT AM EN BASE LOCALE
! attention on ne fait plus d inversion de convention
! donc prevu pour fonctionner avec w mais pas avec u
!
!
! --- CALCUL VECTEUR ET MATRICE TANGENTE EN BASE LOCALE
!
    vim(1) = cohes(1)
    if (nint(rela) .eq. 1) then
        vim(2) = cohes(2)
    else
        if (cohes(2) .le. 0.d0) then
            vim(2) = 0.d0
        else
            vim(2) = 1.d0
        end if
        vim(3) = abs(cohes(2))-1.d0
    end if
!
! PREDICTION: COHES(3)=1 ; CORRECTION: COHES(3)=2
!
    if (nint(cohes(3)) .eq. 1) then
        option = 'RIGI_MECA_TANG'
    else if (nint(cohes(3)) .eq. 2) then
        option = 'FULL_MECA'
    else
        option = 'FULL_MECA'
    end if
!
! VIM = VARIABLES INTERNES UTILISEES DANS LCEJEX
!.............VIM(1): SEUIL, PLUS GRANDE NORME DU SAUT
!
    call lcecli('RIGI', ipgf, 1, ndim, imate, &
                option, lamb, wsaut, sigma, dsidep, &
                vim, vip, raug)
!
    alpha(1) = vip(1)
    if (nint(rela) .eq. 1) then
        alpha(2) = vip(2)
    else
        if (vip(2) .eq. 0.d0) then
            alpha(2) = -vip(3)-1.d0
        else if (vip(2) .eq. 1.d0) then
            alpha(2) = vip(3)+1.d0
        else
            ASSERT(.false.)
        end if
    end if
! --- VARIABLES INTERNES HYDRAULIQUES
    if (option .eq. 'FULL_MECA') then
!   CALCUL DE LA VARIABLE INTERNE : MASSE VOLUMIQUE
        varbio = dpf*cliq
        if (varbio .gt. 5.d0) then
            ASSERT(.false.)
        end if
!
        rho11 = rho11m*exp(varbio)
!
!   CALCUL DE LA VARIABLE INTERNE : APPORTS MASSIQUES
!   (SEULEMENT UTILE POUR LE CAS DU SECOND-MEMBRE)
!
        w11 = w11m+rho11*wsaut(1)-rho11m*wsautm(1)
    end if
    alpha(4) = rho11-rho110
    alpha(5) = w11
! ici on a enleve la securite numerique
!
    if (job .eq. 'ACTU_VI') then
        alpha(3) = 1.d0
    else if (job .eq. 'MATRICE') then
        alpha(3) = 2.d0
    end if
!
end subroutine
