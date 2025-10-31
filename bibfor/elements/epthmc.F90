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
subroutine epthmc(fami, nbEpsi, npg, ndim, &
                  time, anglNaut, jvMaterCode, &
                  indxVarcStrain, epsiVarc)
!
    use BehaviourStrain_module
    use BehaviourStrain_type
    implicit none
!
#include "asterc/r8vide.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/epstmc.h"
#include "asterfort/get_elas_id.h"
!
    character(len=*), intent(in) :: fami
    integer(kind=8), intent(in) ::  nbEpsi, npg, ndim
    real(kind=8), intent(in) :: time, anglNaut(3)
    integer(kind=8), intent(in) :: jvMaterCode, indxVarcStrain
    real(kind=8), intent(out) :: epsiVarc(nbEpsi*npg)
!
! --------------------------------------------------------------------------------------------------
!
! Compute anelastic strains from external state variables
!
! --------------------------------------------------------------------------------------------------
!
! In  fami             : Gauss family for integration point rule
! In  nno              : number of nodes
! In  ndim             : dimension of space
! In  nbEpsi           : number of strain tensor components
! In  npg              : number of Gauss points
! In  anglNaut         : nautical angles (for non-isotropic materials)
! In  time             : given time
! In  jvMaterCode      : coded material address
! In  indxVarcStrain   : index of external state variable
! Out epsiVarc         : anelastic strains from all external state variables
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: ksp = 1
    integer(kind=8) :: kpg, iEpsi, elasID
    character(len=16) :: elasKeyword
    real(kind=8) :: epsiVarcKpg(6)
    type(All_Varc_Strain) :: allVarcStrain
!
! --------------------------------------------------------------------------------------------------
!
    epsiVarc = 0.d0
    ASSERT(nbEpsi .le. 6)

! - Current time
    allVarcStrain%time = time
    if (time .eq. r8vide()) then
        allVarcStrain%hasTime = ASTER_FALSE
    else
        allVarcStrain%hasTime = ASTER_TRUE
    end if

! - Get type of elasticity (Isotropic/Orthotropic/Transverse isotropic)
    call get_elas_id(jvMaterCode, elasID, elasKeyword)

! - Loop on Gauss points
    do kpg = 1, npg
        epsiVarcKpg = 0.d0
        call epstmc(fami, "+", kpg, ksp, ndim, &
                    time, anglNaut, jvMaterCode, &
                    indxVarcStrain, allVarcStrain, &
                    epsiVarcKpg)
        do iEpsi = 1, nbEpsi
            epsiVarc(nbEpsi*(kpg-1)+iEpsi) = epsiVarcKpg(iEpsi)
        end do
    end do
!
end subroutine
