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

subroutine lc0054(fami, kpg, ksp, ndim, imate, &
                  compor, crit, instam, instap, epsm, &
                  deps, sigm, vim, option, angmas, &
                  sigp, vip, typmod, icomp, &
                  nvi, dsidep, codret)
! aslint: disable=W1504,W0104
    implicit none
#include "asterfort/nmtevp.h"
#include "asterfort/rccoma.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: imate, ndim, kpg, ksp, codret
    real(kind=8) :: crit(*), angmas(3)
    real(kind=8) :: instam, instap
    integer(kind=8) :: icomp, nvi
    real(kind=8) :: epsm(6), deps(6)
    real(kind=8) :: sigm(6), sigp(6)
    real(kind=8) :: vim(*), vip(*)
    real(kind=8) :: dsidep(6, 6)
    character(len=16) :: compor(*), option
    character(len=8) :: typmod(*)
    character(len=*) :: fami
!
!
!.......................................................................
!
!     BUT: LOI DE COMPORTEMENT A ECROUISSAGE ISOTROPE DEPENDANT DE LA
!          VITESSE DE DEFORMATION
!
!          RELATIONS : 'VMIS_JOHN_COOK'
!
!       IN      FAMI    FAMILLE DE POINT DE GAUSS (RIGI,MASS,...)
!       IN      KPG,KSP NUMERO DU (SOUS)POINT DE GAUSS
!       IN      NDIM    DIMENSION DE L ESPACE (3D=3,2D=2,1D=1)
!               TYPMOD  TYPE DE MODELISATION
!               IMATE    ADRESSE DU MATERIAU CODE
!               COMPOR    COMPORTEMENT DE L ELEMENT
!               INSTAM   INSTANT T
!               INSTAP   INSTANT T+DT
!               EPSM   DEFORMATION TOTALE A T
!               DEPS   INCREMENT DE DEFORMATION TOTALE
!               SIGM    CONTRAINTE A T
!               VIM    VARIABLES INTERNES A T    + INDICATEUR ETAT T
!    ATTENTION  VIM    VARIABLES INTERNES A T MODIFIEES SI REDECOUPAGE
!               OPTION     OPTION DE CALCUL A FAIRE
!               ANGMAS
!       OUT     SIGP    CONTRAINTE A T+DT
!               VIP    VARIABLES INTERNES A T+DT + INDICATEUR ETAT T+DT
!               DSIDEP    MATRICE DE COMPORTEMENT TANGENT A T+DT OU T
!.......................................................................
!               CODRET
    character(len=16) :: mcmate
    integer(kind=8) :: iret
    real(kind=8) :: r8bid
!
    call rccoma(imate, 'ELAS', 1, mcmate, iret)
!
    if (mcmate .eq. 'ELAS') then
        call nmtevp(fami, kpg, ksp, ndim, typmod, &
                    imate, compor, crit, instam, instap, &
                    deps, sigm, vim, option, sigp, &
                    vip, dsidep, r8bid, r8bid, codret)
    else
        call utmess('F', 'ALGORITH6_88')
    end if
!
end subroutine
