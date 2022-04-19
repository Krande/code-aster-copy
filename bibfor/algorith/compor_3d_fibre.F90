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
subroutine compor_3d_fibre(fami,    kpg,    ksp,  option, sigx,   &
                           instam,  instap, crit, icdmat, materi, &
                           pcompor, epsx,   depx, angmas, vim,    &
                           vip,     sigxp,  etan, codret)
!
! --------------------------------------------------------------------------------------------------
!
!       Intégration de lois de comportement non linéaires pour des éléments 1D
!       par une méthode inspirée de celle de DE BORST en contraintes planes.
!
!       Permet d'utiliser les comportements développés en axis pour traiter des problèmes
!       1D (barres, pmf,...)
!
!       En entrée on donne les valeurs uni-axiales à l'instant précédent :
!       - sigx(t-),epsx(t-),vim(t-)
!       - l'incrément depsx
!       - les variables internes à l'itération précédente,
!
!       En sortie, on obtient :
!       - les contraintes sigxp(t+),
!       - le module tangent etan(t+)
!       - les variables internes vip(t+)
!       - le code retour codret.
!
!       Le code retour est transmis à la routine nmconv pour ajouter des itérations
!       si les contraintes sigyy ou sigzz ne sont pas nulles.
!
! --------------------------------------------------------------------------------------------------
!
use Behaviour_type
use Behaviour_module
!
implicit none
!
#include "asterfort/nmcomp.h"
!
    integer             :: codret, kpg, ksp, icdmat
    real(kind=8)        :: instam, instap
    real(kind=8)        :: sigx, sigxp, epsx, depx, etan
    real(kind=8)        :: angmas(3), vim(*), vip(*), crit(*)
    character(len=*)    :: fami
    character(len=8)    :: materi
    character(len=16)   :: option, pcompor(*)
!
! --------------------------------------------------------------------------------------------------
!
! IN    fami    : famille du point de gauss
!       kpg     : numéro du point de gauss
!       ksp     : numéro du sous-point de gauss
!       option  : nom de l'option à calculer
!       sigx    : sigma xx à l'instant moins
!       epsx    : epsi xx à l'instant moins
!       depx    : delta-epsi xx à l'instant actuel
!       instam  : instant moins
!       instap  : instant plus
!       crit    : critères de convergence locaux
!       icdmat  : code matériau
!       materi  : nom du matériau
!       vim     : variables internes à l'instant moins
!       angmas  : angle massif, loi orthotrope
!
! IN/OUT
!   IN  vip     : variables internes à l'itération précédente
!   OUT vip     : variables internes à l'instant plus
!
! OUT   sigxp   : contraintes à l'instant plus
!       etan    : module tangent à l'instant plus
!       codret  : code retour non nul si les contraintes sigyy ou sigzz sont non nulles
!
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: dsidep(6, 6)
    real(kind=8) :: sigm(6), sigp(6), eps(6), deps(6)
!
    character(len=8)        :: typmod(3)
!
    type(Behaviour_Integ)   :: BEHinteg
!
! --------------------------------------------------------------------------------------------------
!
    eps(:) = 0.0; sigm(:) = 0.0; sigp(:) = 0.0; deps(:) = 0.0; dsidep(:,:) = 0.0
    !
    eps(1) = epsx; deps(1) = depx; sigm(1) = sigx
    typmod(1) = 'COMP1D  '
    typmod(2) = '        '
    !
    ! Initialisation de la SD du comportement
    call behaviourInit(BEHinteg)
    !
    ! Appel à la loi de comportement
    call nmcomp(BEHinteg, &
                fami,   kpg,     ksp,    2,      typmod, &
                icdmat, pcompor, crit,   instam, instap, &
                6,      eps,     deps,   6,      sigm,   &
                vim,    option,  angmas, sigp,   vip,    &
                36,     dsidep,  codret, materi_= materi)
    !
    sigxp = sigp(1)
    etan  = dsidep(1,1)
!
end subroutine
