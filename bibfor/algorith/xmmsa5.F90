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

subroutine xmmsa5(ndim, ipgf, imate, saut, lamb, &
                  nd, tau1, tau2, cohes, job, &
                  rela, alpha, dsidep, delta, p, &
                  am, r)
    implicit none
#include "jeveux.h"
#include "asterc/r8prem.h"
#include "asterfort/lceiou.h"
#include "asterfort/lceitc.h"
#include "asterfort/prmave.h"
    integer :: ndim, ipgf, imate
    real(kind=8) :: saut(3), am(3), dsidep(6, 6)
    real(kind=8) :: tau1(3), tau2(3), nd(3)
    real(kind=8) :: alpha(3), p(3, 3)
    real(kind=8) :: cohes(3), rela, r
    character(len=8) :: job
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
    integer :: i, ier
!
    real(kind=8) :: vim(9), vip(9), lamb(3)
    real(kind=8) :: delta(6), eps
!
    character(len=16) :: option
!
! ----------------------------------------------------------------------
!
! --- INIIALISATIONS
!
    am(:) = 0.d0
    p(:, :) = 0.d0
    dsidep(:, :) = 0.d0
    delta(:) = 0.d0
    vim(:) = 0.d0
    vip(:) = 0.d0
!
! --- ON CONSTRUIT P MATRICE DE PASSAGE BASE FIXE --> BASE COVARIANTE
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
! --- CALCUL SAUT DE DEPLACEMENT EN BASE LOCALE {AM}=[P]{SAUT}
! --- ON UTILISE L UTILITAIRE PRMAVE : PRODUIT MATRICE-VECTEUR
!
    call prmave(0, p, 3, ndim, ndim, &
                saut, ndim, am, ndim, ier)
!
! --- INVERSION DE CONVENTIONS ENTRE X-FEM ET ROUTINE COMPORTEMENT
!
    do i = 1, ndim
        am(i) = -am(i)
    end do
!
! SI ON VEUT SIMPLEMENT LE SAUT LOCAL, ON S ARRETE ICI
!
    if (job .ne. 'SAUT_LOC') then
!
! --- CALCUL VECTEUR ET MATRICE TANGENTE EN BASE LOCALE
!
        vim(4) = cohes(1)
        vim(2) = cohes(2)
!
! --- PREDICTION: COHES(3)=1, CORRECTION: COHES(3)=2
!
        if (cohes(3) .eq. 1.d0) then
            option = 'RIGI_MECA_TANG'
        else if (cohes(3) .eq. 2.d0) then
            option = 'FULL_MECA'
        else
            option = 'FULL_MECA'
        end if

        if (job .eq. 'MATRICE' .and. option .eq. 'FULL_MECA') then
            eps = 100.*r8prem()
            vim(4) = min(1.d0, vim(4)*(1+eps))
        end if
!
! VIM = VARIABLES INTERNES UTILISEES DANS LCEJEX
!.............VIM(1): SEUIL, PLUS GRANDE NORME DU SAUT
!
        if (rela .eq. 3.d0) then
            call lceitc('RIGI', ipgf, 1, imate, option, &
                        lamb, am, delta, dsidep, vim, &
                        vip, r)
        else if (rela .eq. 4.d0) then
            call lceiou('RIGI', ipgf, 1, imate, option, &
                        lamb, am, delta, dsidep, vim, &
                        vip, r)
!
        end if
!
! VARIABLES INTERNES ACTUALISEES
!
        alpha(1) = vip(4)
        alpha(2) = vip(2)
! SI ACTUALISATION: NOUVEAU PAS DONC PREDICTION EN PERSPECTIVE
! SINON, DESCENTE
        if (job .eq. 'ACTU_VI') then
            alpha(3) = 1.d0
        else if (job .eq. 'MATRICE') then
            alpha(3) = 2.d0
        end if
!
    end if
end subroutine
