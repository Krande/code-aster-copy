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
! person_in_charge: mickael.abbas at edf.fr
! aslint: disable=W1504
!
subroutine mmmvee(phase, l_pena_cont, l_pena_fric, l_large_slip, &
                  ndim, nne, &
                  norm, tau1, tau2, mprojt, &
                  wpg, ffe, jacobi, jeu, &
                  coefac, coefaf, lambda, coefff, &
                  dlagrc, dlagrf, djeu, &
                  rese, nrese, &
                  mprt11, mprt12, mprt21, mprt22, kappa, &
                  vectce, vectfe)
!
    implicit none
!
#include "asterf_types.h"
!
    character(len=4), intent(in) :: phase
    aster_logical, intent(in) :: l_pena_cont, l_pena_fric, l_large_slip
    integer(kind=8), intent(in) :: ndim, nne
    real(kind=8), intent(in) :: norm(3), tau1(3), tau2(3), mprojt(3, 3)
    real(kind=8), intent(in) :: wpg, ffe(9), jacobi, jeu
    real(kind=8), intent(in) :: coefac, coefaf, lambda, coefff
    real(kind=8), intent(in) :: dlagrc, dlagrf(2), djeu(3)
    real(kind=8), intent(in) :: rese(3), nrese
    real(kind=8), intent(in) :: mprt11(3, 3), mprt12(3, 3), mprt21(3, 3), mprt22(3, 3), kappa(2, 2)
    real(kind=8), intent(out) :: vectce(27), vectfe(27)
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Elementary computations
!
! Compute vector for DOF [slave]
!
! --------------------------------------------------------------------------------------------------
!
! In  phase            : phase to compute
!                        'SANS' - No contact
!                        'ADHE' - Stick
!                        'GLIS' - Slip
! In  l_pena_cont      : flag for penalized contact
! In  l_pena_fric      : flag for penalized friction
! In  l_large_slip     : flag for GRAND_GLISSEMENT
! In  ndim             : dimension of problem (2 or 3)
! In  nne              : number of slave nodes
! In  norm             : normal at current contact point
! In  tau1             : first tangent at current contact point
! In  tau2             : second tangent at current contact point
! In  mprojt           : matrix of tangent projection
! In  wpg              : weight for current Gauss point
! In  ffe              : shape function for slave nodes
! In  jacobi           : jacobian at integration point
! In  jeu              : normal gap
! In  coefac           : coefficient for updated Lagrangian method (contact)
! In  coefaf           : coefficient for updated Lagrangian method (friction)
! In  lambda           : contact pressure
! In  coefff           : friction coefficient (Coulomb)
! In  dlagrc           : increment of contact Lagrange from beginning of time step
! In  dlagrf           : increment of friction Lagrange from beginning of time step
! In  djeu             : increment of gap
! In  rese             : Lagrange (semi) multiplier for friction
! In  nrese            : norm of Lagrange (semi) multiplier for friction
! In  mprt11           : projection matrix first tangent/first tangent
!                        tau1*TRANSPOSE(tau1)(matrice 3*3)
! In  mprt12           : projection matrix first tangent/second tangent
!                        tau1*TRANSPOSE(tau2)(matrice 3*3)
! In  mprt21           : Projection matrix second tangent/first tangent
!                        tau2*TRANSPOSE(tau1)(matrice 3*3)
! In  mprt22           : Projection matrix second tangent/second tangent
!                        tau2*TRANSPOSE(tau2)(matrice 3*3)
! In  kappa            : scalar matrix for sliding kinematic
! Out vectce           : vector for DOF [slave] - For contact part
! Out vectfe           : vector for DOF [slave] - For friction part
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: inoe, idim, ii, i, j, k
    real(kind=8) :: dlagft(3), plagft(3), prese(3)
    real(kind=8) :: dvitet(3), pdvitt(3), g(3, 3)
!
! --------------------------------------------------------------------------------------------------
!
    plagft = 0.d0
    dlagft = 0.d0
    prese = 0.d0
    dvitet = 0.d0
    pdvitt = 0.d0
    vectce = 0.d0
    vectfe = 0.d0
!
! - PROJECTION DU LAGRANGE DE FROTTEMENT SUR LE PLAN TANGENT
!
    do i = 1, ndim
        dlagft(i) = dlagrf(1)*tau1(i)+dlagrf(2)*tau2(i)
    end do
!
! - PRODUIT LAGR. FROTTEMENT. PAR MATRICE P
!
    do i = 1, ndim
        do j = 1, ndim
            plagft(i) = mprojt(i, j)*dlagft(j)+plagft(i)
        end do
    end do
!
! - PRODUIT SEMI MULT. LAGR. FROTTEMENT. PAR MATRICE P
!
    if (phase .eq. 'GLIS') then
        do i = 1, ndim
            do j = 1, ndim
                if (l_large_slip) then
                    g(i, j) = kappa(1, 1)*mprt11(i, j)+ &
                              kappa(1, 2)*mprt12(i, j)+ &
                              kappa(2, 1)*mprt21(i, j)+ &
                              kappa(2, 2)*mprt22(i, j)
                    prese(i) = g(i, j)*rese(j)/nrese+prese(i)
                else
                    prese(i) = mprojt(i, j)*rese(j)/nrese+prese(i)
                end if
            end do
        end do
    end if
!
! - PROJECTION DU SAUT SUR LE PLAN TANGENT
!
    do i = 1, ndim
        do k = 1, ndim
            dvitet(i) = mprojt(i, k)*djeu(k)+dvitet(i)
        end do
    end do
!
! - PRODUIT SAUT PAR MATRICE P
!
    do i = 1, ndim
        do j = 1, ndim
            pdvitt(i) = mprojt(i, j)*dvitet(j)+pdvitt(i)
        end do
    end do
!
! - Compute terms
!
    if (phase .ne. 'SANS') then
        if (l_pena_cont) then
            do inoe = 1, nne
                do idim = 1, ndim
                    ii = ndim*(inoe-1)+idim
                    vectce(ii) = vectce(ii)+ &
                                 wpg*ffe(inoe)*jacobi*norm(idim)*(jeu*coefac)
                end do
            end do
        else
            do inoe = 1, nne
                do idim = 1, ndim
                    ii = ndim*(inoe-1)+idim
                    vectce(ii) = vectce(ii)- &
                                 wpg*ffe(inoe)*jacobi*norm(idim)*(dlagrc-jeu*coefac)
                end do
            end do
        end if
    end if
    if (phase .eq. 'GLIS') then
        do inoe = 1, nne
            do idim = 1, ndim
                ii = ndim*(inoe-1)+idim
                vectfe(ii) = vectfe(ii)- &
                             wpg*ffe(inoe)*jacobi*prese(idim)*(lambda-0.*jeu)*coefff
            end do
        end do
    end if
    if (phase .eq. 'ADHE') then
        if (l_pena_fric) then
            do inoe = 1, nne
                do idim = 1, ndim
                    ii = ndim*(inoe-1)+idim
                    vectfe(ii) = vectfe(ii)- &
                                 wpg*ffe(inoe)*jacobi*lambda*coefff*pdvitt(idim)*coefaf
                end do
            end do
        else
            do inoe = 1, nne
                do idim = 1, ndim
                    ii = ndim*(inoe-1)+idim
                    vectfe(ii) = vectfe(ii)- &
                                 wpg*ffe(inoe)*jacobi*(lambda-0.*jeu)*coefff* &
                                 (plagft(idim)+pdvitt(idim)*coefaf)
                end do
            end do
        end if
    end if
!
end subroutine
