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
!
subroutine xmmsa2(ndim, ipgf, imate, saut, nd, &
                  tau1, tau2, cohes, job, rela, &
                  alpha, dsidep, sigma, pp, dnor, &
                  dtang, p, am)
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/lcejex.h"
#include "asterfort/lcejli.h"
#include "asterfort/xmafr1.h"
    integer :: ndim, ipgf, imate
    real(kind=8) :: saut(3), am(3), pp(3, 3), dsidep(6, 6)
    real(kind=8) :: tau1(3), tau2(3), nd(3)
    real(kind=8) :: alpha(3), p(3, 3)
    real(kind=8) :: dtang(3), dnor(3), cohes(3), rela
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
!
!
    integer :: i, k
!
    real(kind=8) :: vim(9), vip(9)
    real(kind=8) :: am2d(2), dam2d(2), dsid2d(6, 6), dam(3)
    real(kind=8) :: sigma(6)
!
    character(len=16) :: option
!
! ----------------------------------------------------------------------
!
    am(:) = 0.d0
    dam(:) = 0.d0
    pp(:, :) = 0.d0
    p(:, :) = 0.d0
    dsidep(:, :) = 0.d0
    dsid2d(:, :) = 0.d0
    sigma(:) = 0.d0
    vim(:) = 0.d0
    vip(:) = 0.d0
    dtang(:) = 0.d0
    dnor(:) = 0.d0
!
    call xmafr1(3, nd, p)
!
! --- CALCUL DU SAUT DE DEPLACEMENT AM EN BASE LOCALE
!
    do i = 1, ndim
        dtang(i) = 0.d0
        dnor(i) = 0.d0
        do k = 1, ndim
            pp(i, k) = nd(i)*nd(k)
            dtang(i) = dtang(i)+p(i, k)*saut(k)
            dnor(i) = dnor(i)+pp(i, k)*saut(k)
        end do
!
! --- L'INTERPENETRATION CORRESPOND A SAUT<0 DANS LCEJEX
!
        am(1) = am(1)-dnor(i)*nd(i)
        am(2) = am(2)-dtang(i)*tau1(i)
        am(3) = am(3)-dtang(i)*tau2(i)
    end do
!
! --- CALCUL VECTEUR ET MATRICE TANGENTE EN BASE LOCALE
!
    if (job .ne. 'SAUT_LOC') then
        vim(1) = cohes(1)
        if (cohes(2) .le. 0.d0) then
            vim(2) = 0.d0
        else
            vim(2) = 1.d0
        end if
        vim(3) = abs(cohes(2))-1.d0
!
! PREDICTION: COHES(3)=1 ; CORRECTION: COHES(3)=2
!
        if (cohes(3) .eq. 1.d0) then
            option = 'RIGI_MECA_TANG'
        else if (cohes(3) .eq. 2.d0) then
            option = 'FULL_MECA'
        else
            option = 'FULL_MECA'
        end if
!
! VIM = VARIABLES INTERNES UTILISEES DANS LCEJEX
!.............VIM(1): SEUIL, PLUS GRANDE NORME DU SAUT
!
        if (ndim .eq. 2) then
            am2d(1) = am(1)
            am2d(2) = am(2)
            dam2d(1) = dam(1)
            dam2d(2) = dam(2)
            if (rela .eq. 1.d0) then
                call lcejex('RIGI', ipgf, 1, 2, imate, &
                            option, am2d, dam2d, sigma, dsid2d, &
                            vim, vip)
            else if (rela .eq. 2.d0) then
                call lcejli('RIGI', ipgf, 1, 2, imate, &
                            option, am2d, dam2d, sigma, dsid2d, &
                            vim, vip)
            end if
            dsidep(1, 1) = dsid2d(1, 1)
            dsidep(1, 2) = dsid2d(1, 2)
            dsidep(2, 1) = dsid2d(2, 1)
            dsidep(2, 2) = dsid2d(2, 2)
        else if (ndim .eq. 3) then
            if (rela .eq. 1.d0) then
                call lcejex('RIGI', ipgf, 1, ndim, imate, &
                            option, am, dam, sigma, dsidep, &
                            vim, vip)
            else if (rela .eq. 2.d0) then
                call lcejli('RIGI', ipgf, 1, ndim, imate, &
                            option, am, dam, sigma, dsidep, &
                            vim, vip)
            end if
        end if
!
        alpha(1) = vip(1)
        if (vip(2) .eq. 0.d0) then
            alpha(2) = -vip(3)-1.d0
        else if (vip(2) .eq. 1.d0) then
            alpha(2) = vip(3)+1.d0
        else
            ASSERT(.false.)
        end if
! ici on a enleve la securite numerique
!
        if (job .eq. 'ACTU_VI') then
            alpha(3) = 1.d0
        else if (job .eq. 'MATRICE') then
            alpha(3) = 2.d0
        end if
!
    end if
end subroutine
