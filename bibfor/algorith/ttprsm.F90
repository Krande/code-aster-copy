! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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
subroutine ttprsm(ndim, ddeple, ddeplm, dlagrf, coeffr,&
                  tau1, tau2, mprojt, inadh, rese,&
                  nrese, coeffp, lpenaf, dvitet)
!
! person_in_charge: mickael.abbas at edf.fr
!
    implicit none
#include "asterf_types.h"
#include "asterfort/assert.h"
    integer :: ndim
    real(kind=8) :: ddeple(3), ddeplm(3), dlagrf(2)
    real(kind=8) :: coeffr, coeffp
    real(kind=8) :: tau1(3), tau2(3), mprojt(3, 3)
    integer :: inadh
    real(kind=8) :: rese(3), nrese
    aster_logical :: lpenaf
!
! ----------------------------------------------------------------------
!
! ROUTINE CONTACT (METHODE CONTINUE - UTILITAIRE)
!
! ETAT D'ADHERENCE DU POINT DE CONTACT
!
! ----------------------------------------------------------------------
!
!
!  CALCUL DE P_B(0,1)(LAMDBA+ RHO C [[DELTA X]])
!
! IN  NDIM   : DIMENSION DE L'ESPACE
! IN  MPROJT : MATRICE DE PROJECTION TANGENTE
! IN  DLAGRF : INCREMENT DEPDEL DES LAGRANGIENS DE FROTTEMENT
! IN  DDEPLE : INCREMENT DEPDEL DU DEPL. DU POINT DE CONTACT
! IN  DDEPLM : INCREMENT DEPDEL DU DEPL. DU PROJETE DU POINT DE CONTACT
! IN  COEFFR : COEF_REGU_FROT
! IN  TAU1   : PREMIER VECTEUR TANGENT
! IN  TAU2   : SECOND VECTEUR TANGENT
! OUT INADH  : INDICE D'ADHERENCE
!               1 - ADHERENT
!               0 - GLISSANT
! OUT RESE   : SEMI-MULTIPLICATEUR GTK DE FROTTEMENT
!               GTK = LAMBDAF + COEFFR*VITESSE
! OUT NRESE  : NORME DU SEMI-MULTIPLICATEUR GTK DE FROTTEMENT
! OUT DVITET : VITESSE DE GLISSEMENT
!
! ----------------------------------------------------------------------
!
    integer :: i, k
    real(kind=8) :: dvite(3), dvitet(3)
!
! ----------------------------------------------------------------------
!
! --- INITIALISATIONS
!
    nrese = 0.d0
    do i = 1, 3
        rese(i) = 0.d0
        dvitet(i) = 0.d0
    end do
!
! --- CALCUL DU SAUT DE "VITESSE" [[DELTA X]]
!
    do i = 1, ndim
        dvite(i) = ddeple(i) - ddeplm(i)
    end do
!
! --- PROJECTION DU SAUT SUR LE PLAN TANGENT
!
    do i = 1, ndim
        do k = 1, ndim
            dvitet(i) = mprojt(i,k)*dvite(k)+dvitet(i)
        end do
    end do
!
! --- SEMI-MULTIPLICATEUR DE FROTTEMENT RESE
!
    if (lpenaf) then
        do i = 1, 3
            rese(i) = coeffp*dvitet(i)
        end do
    else
        if (ndim .eq. 2) then
            do i = 1, 2
                rese(i) = dlagrf(1)*tau1(i)+coeffr*dvitet(i)
            end do
        else if (ndim.eq.3) then
            do i = 1, 3
                rese(i) = dlagrf(1)*tau1(i)+ dlagrf(2)*tau2(i)+ coeffr*dvitet(i)
            end do
        else
            ASSERT(.false.)
        endif
    endif
!
! -- CALCUL DU COEF D'ADHERENCE
!
    do i = 1, 3
        nrese = rese(i)*rese(i) + nrese
    end do
    nrese = sqrt(nrese)
!
! --- ON TESTE SI    NRESE < 1 OU NON
! --- SI OUI ADHERENCE
!
    if (nrese .le. 1.d0) then
        inadh = 1
    else
        inadh = 0
    endif
!
end subroutine
