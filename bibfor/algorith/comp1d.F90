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
subroutine comp1d(BEHinteg, &
                  fami, kpg, ksp, option, sigx, &
                  epsx, depx, angmas, vim, vip, &
                  sigxp, etan, codret)
!
    use Behaviour_type
    use Behaviour_module
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/jevech.h"
#include "asterfort/nmcomp.h"
!
    type(Behaviour_Integ), intent(in) :: BEHinteg
    character(len=*) :: fami
    character(len=16) :: option
    integer(kind=8) :: codret, kpg, ksp
    real(kind=8) :: angmas(3)
    real(kind=8) :: vim(*), vip(*), sigx, sigxp, epsx, depx, etan
!
!     INTEGRATION DE LOIS DE COMPORTEMENT NON LINEAIRES
!     POUR DES ELEMENTS 1D PAR UNE METHODE INSPIREE DE CELLE DE DEBORST
!     EN CONTRAINTES PLANES.
!
!     PERMET D'UTILISER TOUS LES COMPORTEMENTS DEVELOPPES EN AXIS POUR
!     TRAITER DES PROBLEMES 1D (BARRES, PMF,...)
!
!     POUR POUVOIR UTILISER CETTE METHODE, IL FAUT FOURNIR SOUS LE
!     MOT-CLES COMP_INCR : ALGO_1D='DEBORST'
!
!     EN ENTREE ON DONNE LES VALEURS UNIAXIALES A L'INSTANT PRECEDENT :
!      - SIGX(T-),EPSX(T-),VIM(T-) ET L'INCREMENT DEPSX
!      - ET LES VARIABLES INTERNES A L'ITERATION PRECEDENTE,
!     RECOPIEES DANS VIP
!     EN SORTIE, ON OBTIENT LES CONTRAINTES SIGXP(T+), LE MODULES
!     TANGENT ETAN, LES VARAIBLES INTERNES VIP(T+) ET LE
!     CODE RETOUR CODRET.
!   EN SORTIE DE COMP1D, CE CODE RETOUR EST TRANSMIS A LA ROUTINE NMCONV
!   POUR AJOUTER DES ITERATIONS SI SIGYY OU SIGZZ NE SONT PAS NULLES.
!
! ----------------------------------------------------------------------
! IN  FAMI      : FAMILLE DU POINT DE GAUSS
!     KPG       : NUMERO DU POINT DE GAUSS
!     KSP       : NUMERO DU SOUS-POINT DE GAUSS
!     OPTION    : NOM DE L'OPTION A CALCULER
!     SIGM      : SIGMA XX A L'INSTANT MOINS
!     EPSX      : EPSI XX A L'INSTANT MOINS
!     DEPX      : DELTA-EPSI XX A L'INSTANT ACTUEL
!     VIM       : VARIABLES INTERNES A L'INSTANT MOINS
! VAR VIP       : VARIABLES INTERNES A L'INSTANT PLUS EN SORTIE
!                 VARIABLES INTERNES A L'ITERATION PRECEDENTE EN ENTREE
! OUT SIGXP     : CONTRAINTES A L'INSTANT PLUS
!     ETAN      : MODULE TANGENT DIRECTION X
!     CODRET    : CODE RETOUR NON NUL SI SIGYY OU SIGZZ NON NULS
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: imate, iinstm
    integer(kind=8) :: iinstp, icarcr
    integer(kind=8), parameter :: ndimLdc = 2
    real(kind=8) :: dsidep(6, 6)
    real(kind=8) :: sigm(6), sigp(6), eps(6), deps(6)
    character(len=8) :: typmod(2)
    character(len=16), pointer :: compor(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    eps = 0.d0
    sigm = 0.d0
    sigp = 0.d0
    deps = 0.d0
    dsidep = 0.d0
    eps(1) = epsx
    deps(1) = depx
    sigm(1) = sigx
    typmod(1) = 'COMP1D '
    typmod(2) = '        '
!
!
! ---    PARAMETRES EN ENTREE
    call jevech('PMATERC', 'L', imate)
    call jevech('PINSTMR', 'L', iinstm)
    call jevech('PINSTPR', 'L', iinstp)
    call jevech('PCOMPOR', 'L', vk16=compor)
    call jevech('PCARCRI', 'L', icarcr)

! - Integrator
    sigp = 0.d0
    call nmcomp(BEHinteg, &
                fami, kpg, ksp, ndimLdc, typmod, &
                zi(imate), compor, zr(icarcr), zr(iinstm), zr(iinstp), &
                6, eps, deps, 6, sigm, &
                vim, option, angmas, &
                sigp, vip, 36, dsidep, codret)
!
    sigxp = sigp(1)
    etan = dsidep(1, 1)
!
!
end subroutine
