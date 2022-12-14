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
! aslint: disable=W1504
!
subroutine nmcomp(BEHinteg, &
                  fami,   kpg,    ksp,    ndim,       typmod,        &
                  imate,  compor, carcri, instam,     instap,        &
                  neps,   epsm,   deps,   nsig,       sigm,          &
                  vim,    option, angmas, sigp,       vip,           &
                  ndsde,  dsidep, codret, mult_comp_, l_epsi_varc_,  &
                  materi_)
!
use Behaviour_type
!
implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/lcvali.h"
#include "asterfort/nmcpl1.h"
#include "asterfort/nmcpl2.h"
#include "asterfort/nmcpl3.h"
#include "asterfort/redece.h"
#include "asterfort/lcidbg.h"
!
type(Behaviour_Integ) :: BEHinteg
!
integer :: kpg, ksp, ndim, imate, codret, neps, nsig, ndsde
!
character(len=*)    :: fami
character(len=8)    :: typmod(*)
character(len=16)   :: compor(*), option
!
real(kind=8) :: instam, instap
real(kind=8) :: epsm(*), deps(*), dsidep(*), carcri(*), sigm(*), vim(*), sigp(*), vip(*), angmas(*)
!
character(len=8),  optional, intent(in) :: materi_
character(len=16), optional, intent(in) :: mult_comp_
aster_logical,     optional, intent(in) :: l_epsi_varc_
!
! --------------------------------------------------------------------------------------------------
!
!     INTEGRATION DES LOIS DE COMPORTEMENT NON LINEAIRE POUR LES
!     ELEMENTS ISOPARAMETRIQUES EN PETITES OU GRANDES DEFORMATIONS
!
! --------------------------------------------------------------------------------------------------
!
! In  BEHinteg       : parameters for integration of behaviour
! IN  FAMI,KPG,KSP  : FAMILLE ET NUMERO DU (SOUS)POINT DE GAUSS
!     NDIM    : DIMENSION DE L'ESPACE
!               3 : 3D , 2 : D_PLAN ,AXIS OU  C_PLAN
!     TYPMOD(2): MODELISATION ex: 1:3D, 2:INCO
!     IMATE   : ADRESSE DU MATERIAU CODE
!     COMPOR  : COMPORTEMENT :  (1) = TYPE DE RELATION COMPORTEMENT
!                               (2) = NB VARIABLES INTERNES / PG
!                               (3) = HYPOTHESE SUR LES DEFORMATIONS
!                               (4) etc... (voir grandeur COMPOR)
!     CRIT    : CRITERES DE CONVERGENCE LOCAUX (voir grandeur CARCRI)
!     INSTAM  : INSTANT DU CALCUL PRECEDENT
!     INSTAP  : INSTANT DU CALCUL
!     NEPS    : NOMBRE DE CMP DE EPSM ET DEPS (SUIVANT MODELISATION)
!     EPSM    : DEFORMATIONS A L'INSTANT DU CALCUL PRECEDENT
!     DEPS    : INCREMENT DE DEFORMATION TOTALE :
!                DEPS(T) = DEPS(MECANIQUE(T)) + DEPS(DILATATION(T))
!     NSIG    : NOMBRE DE CMP DE SIGM ET SIGP (SUIVANT MODELISATION)
!     SIGM    : CONTRAINTES A L'INSTANT DU CALCUL PRECEDENT
!     VIM     : VARIABLES INTERNES A L'INSTANT DU CALCUL PRECEDENT
!     OPTION  : OPTION DEMANDEE : RIGI_MECA_TANG , FULL_MECA , RAPH_MECA
!     ANGMAS  : LES TROIS ANGLES DU MOT_CLEF MASSIF (AFFE_CARA_ELEM),
!               + UN REEL QUI VAUT 0 SI NAUTIQUIES OU 2 SI EULER
!               + LES 3 ANGLES D'EULER
!
! OUT SIGP    : CONTRAINTES A L'INSTANT ACTUEL
! VAR VIP     : VARIABLES INTERNES
!                IN  : ESTIMATION (ITERATION PRECEDENTE OU LAG. AUGM.)
!                OUT : EN T+
!     NDSDE   : DIMENSION DE DSIDEP
!     DSIDEP  : OPERATEUR TANGENT DSIG/DEPS OU DSIG/DF
!     CODRET  : CODE RETOUR LOI DE COMPORMENT :
!               CODRET=0 : TOUT VA BIEN
!               CODRET=1 : ECHEC DANS L'INTEGRATION DE LA LOI
!               CODRET=3 : SIZZ NON NUL (CONTRAINTES PLANES DEBORST)
!
! PRECISIONS :
! -----------
!  LES TENSEURS ET MATRICES SONT RANGES DANS L'ORDRE :
!         XX YY ZZ SQRT(2)*XY SQRT(2)*XZ SQRT(2)*YZ
!
! -SI DEFORMATION = SIMO_MIEHE
!   EPSM(3,3)    GRADIENT DE LA TRANSFORMATION EN T-
!   DEPS(3,3)    GRADIENT DE LA TRANSFORMATION DE T- A T+
!
!  OUTPUT SI RESI (RAPH_MECA, FULL_MECA_*)
!   VIP      VARIABLES INTERNES EN T+
!   SIGP(6)  CONTRAINTE DE KIRCHHOFF EN T+ RANGES DANS L'ORDRE
!         XX YY ZZ SQRT(2)*XY SQRT(2)*XZ SQRT(2)*YZ
!
!  OUTPUT SI RIGI (RIGI_MECA_*, FULL_MECA_*)
!   DSIDEP(6,3,3) MATRICE TANGENTE D(TAU)/D(FD) * (FD)T
!                 (AVEC LES RACINES DE 2)
!
! -SINON (DEFORMATION = PETIT OU PETIT_REAC OU GDEF_...)
!   EPSM(6), DEPS(6)  SONT LES DEFORMATIONS (LINEARISEES OU GREEN OU ..)
!
! --------------------------------------------------------------------------------------------------
!
!   Pour les utilitaires de calcul tensoriel
    integer :: ndt, ndi
    common /tdim/ ndt,ndi
!
    integer :: icp, numlc, cpl, nvv, ncpmax
!
    aster_logical :: cp, convcp, l_epsi_varc
    character(len=8)  :: materi
    character(len=16) :: optio2, mult_comp
!
! --------------------------------------------------------------------------------------------------
!
    codret = 0
    l_epsi_varc = ASTER_TRUE
    if (present(l_epsi_varc_)) then
        l_epsi_varc = l_epsi_varc_
    endif
!
!   Contraintes planes ?
    call nmcpl1(compor, typmod, option, vip, deps, optio2, cpl, nvv)
    cp=(cpl.ne.0)
!
!   Dimensionnement pour le calcul tensoriel
    ndt = 2*ndim
    ndi = ndim
!
    if (cp) then
        convcp = ASTER_FALSE
        ncpmax = nint(carcri(ITER_DEBORST_MAX))
    else
        convcp = ASTER_TRUE
        ncpmax = 1
    endif
!
!   Les param??tres optionnels
    mult_comp = ' '
    if (present(mult_comp_)) then
        mult_comp = mult_comp_
    endif
    materi = ' '
    if (present(materi_)) then
        materi = materi_
    endif
!
!   Num??ro de la loi de comportement : numlc
    read(compor(NUME),'(I16)') numlc
!
!   Boucle pour ??tablir les contraintes planes
    do icp = 1, ncpmax
        call redece(BEHinteg,&
                    fami,        kpg,    ksp,    ndim,   typmod,    &
                    l_epsi_varc, imate,  materi, compor, mult_comp, &
                    carcri,      instam, instap, neps,   epsm,      &
                    deps,        nsig,   sigm,   vim,    option,    &
                    angmas,      cp,     numlc,  sigp,   vip,       &
                    ndsde,       dsidep, codret)
        !
        ! V??rifier la convergence des contraintes planes et sortir de la boucle si n??cessaire
        if (cp) then
            ASSERT(ndsde.eq.36)
            call nmcpl3(compor, option, carcri, deps, dsidep, &
                        ndim,   sigp,   vip,    cpl,  icp,    &
                        convcp)
        endif
        !
        if (convcp) then
            exit
        endif
        !
    enddo
!
!   Contraintes planes m??thode DE BORST
    if (cp) then
        if (codret .eq. 0) then
            ASSERT(ndsde.eq.36)
            call nmcpl2(compor, typmod, option, optio2, cpl,  &
                        nvv,    carcri, deps,   dsidep, ndim, &
                        sigp,   vip,    codret)
        else
            option=optio2
        endif
    endif
!   Examen du domaine de validit??
    if (codret .eq. 0) then
        call lcvali(fami,   kpg,  ksp,  imate,  materi, &
                    compor, ndim, epsm, deps,   instam, &
                    instap, codret)
    else if (codret .eq. 1) then
        call lcidbg(fami,   kpg,    ksp,    typmod, compor, &
                    carcri, instam, instap, neps,   epsm,   &
                    deps,   nsig,   sigm,   vim,    option)
    endif
!
end subroutine
