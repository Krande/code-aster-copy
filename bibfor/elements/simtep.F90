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
subroutine simtep(fami, nno, ndim, nbsig, npg, &
                  jvGaussWeight, jvBaseFunc, jvDBaseFunc, &
                  nodeCoor, nodeDisp, &
                  time, anglNaut, jvMaterCode, nharm, &
                  sigmEner)
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/sigmmc.h"
#include "asterfort/sigtmc.h"
#include "jeveux.h"
!
    character(len=*), intent(in) :: fami
    integer(kind=8), intent(in) :: nno, ndim, nbsig, npg
    integer(kind=8), intent(in) :: jvGaussWeight, jvBaseFunc, jvDBaseFunc
    real(kind=8), intent(in) :: nodeCoor(ndim*nno), nodeDisp(ndim*nno)
    real(kind=8), intent(in) :: time, anglNaut(3)
    integer(kind=8), intent(in) :: jvMaterCode
    real(kind=8), intent(in)  :: nharm
    real(kind=8), intent(out) :: sigmEner(nbsig*npg)
!
! --------------------------------------------------------------------------------------------------
!
! Compute "real" stress tensor at Gauss points
!
! CALCUL DES  CONTRAINTES 'VRAIES' POUR LE CALCUL DE L'ENERGIE POTENTIELLE
!                  (I.E.  1/2*SIGMA_MECA - SIGMA_THERMIQUES)
!                  AUX POINTS D'INTEGRATION POUR LES ELEMENTS
!                  ISOPARAMETRIQUES
!
! --------------------------------------------------------------------------------------------------
!
! In  fami             : Gauss family for integration point rule
! In  nno              : number of nodes of element
! In  ndim             : dimension of element (2 ou 3)
! In  nbsig            : number of components for stress tensors (4 or 6)
! In  npg              : number of Gauss points
! In  jvGaussWeight    : adresse to Gauss point weight
! In  jvBaseFunc       : adresse to shape functions
! In  jvDBaseFunc      : adresse to derivative of shape functions
! In  nodeCoor         : coordinates of nodes
! In  nodeDisp         : displacements at nodes
! In  time             : current time
! In  anglNaut         : nautical angles for definition of basis for non-isotropic elasticity
! In  jvMaterCode      : adress for material parameters
! In  nharm            : Fourier mode
! Out sigmEner         : "real" stress tensor at Gauss points
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: sigmTher(162), sigm(162)
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(nbsig .le. 6)
    ASSERT(npg .le. 27)
    sigmEner = 0.d0

! - Compute stress at Gauss points
    call sigmmc(fami, nno, ndim, nbsig, npg, &
                jvGaussWeight, jvBaseFunc, jvDBaseFunc, &
                nodeCoor, nodeDisp, &
                time, anglNaut, jvMaterCode, nharm, &
                sigm)

! - Compute stresses from external state variables
    call sigtmc(fami, nbsig, npg, ndim, &
                time, jvMaterCode, anglNaut, &
                VARC_STRAIN_TEMP, sigmTher)

! - Compute "real" stress
    sigmEner = 0.5d0*sigm-sigmTher

!
end subroutine
