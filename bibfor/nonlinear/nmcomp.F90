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
! aslint: disable=W1504,C1505
!
subroutine nmcomp(BEHinteg, &
                  fami, kpg, ksp, ndim, typmod, &
                  imate, compor, carcri, instam, instap, &
                  neps, epsm, deps, nsig, sigm, &
                  vim, option, angmas, sigp, vip, &
                  ndsde, dsidep, codret, mult_comp_, l_epsi_varc_, &
                  materi_)
!
    use Behaviour_type
    use Behaviour_module
    implicit none
!
#include "asterc/r8vide.h"
#include "asterc/r8gaem.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/lcvali.h"
#include "asterfort/redece.h"
!
    type(Behaviour_Integ) :: BEHinteg
    integer :: kpg, ksp, ndim, imate, codret, neps, nsig, ndsde
    character(len=*)    :: fami
    character(len=8)    :: typmod(*)
    character(len=16)   :: compor(*), option
    real(kind=8) :: instam, instap
    real(kind=8) :: epsm(neps), deps(neps)
    real(kind=8) :: dsidep(merge(nsig, 6, nsig*neps .eq. ndsde), &
                           merge(neps, 6, nsig*neps .eq. ndsde))
    real(kind=8) :: carcri(*), sigm(nsig), vim(*), sigp(nsig), vip(*), angmas(*)
    character(len=8), optional, intent(in) :: materi_
    character(len=16), optional, intent(in) :: mult_comp_
    aster_logical, optional, intent(in) :: l_epsi_varc_
! --------------------------------------------------------------------------------------------------
!     INTEGRATION DES LOIS DE COMPORTEMENT NON LINEAIRE
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
    aster_logical :: conv_cp, l_epsi_varc, lMatr, lVari, lSigm, lMatrPred, lPred, invert
    aster_logical :: l_defo_meca, l_czm, l_large
    integer :: icp, numlc, nvi_all, nvi, k, l, ndimsi
    real(kind=8):: prec
    real(kind=8):: epsm_meca(neps), deps_meca(neps)
    real(kind=8) :: dsidep_cp(merge(nsig, 6, nsig*neps .eq. ndsde), &
                              merge(neps, 6, nsig*neps .eq. ndsde))
    real(kind=8), allocatable:: vip_cp(:), ka3_max, k3a_max, c_max
    character(len=8)  :: materi
    character(len=8)  :: typmod_cp(2), typ_crit
    character(len=16) :: option_cp, mult_comp, defo_ldc, defo_comp
! --------------------------------------------------------------------------------------------------

    ! Controles
    ASSERT(neps*nsig .eq. ndsde .or. (ndsde .eq. 36 .and. neps .le. 6 .and. nsig .le. 6))

!   Les paramètres optionnels
    mult_comp = ' '
    materi = ' '
    l_epsi_varc = ASTER_TRUE
    if (present(mult_comp_)) mult_comp = mult_comp_
    if (present(materi_)) materi = materi_
    if (present(l_epsi_varc_)) l_epsi_varc = l_epsi_varc_

    ! Initialisation
    codret = 0
    read (compor(NUME), '(I16)') numlc
    read (compor(NVAR), '(I16)') nvi_all
    read (compor(DEFO_LDC), '(A16)') defo_ldc
    read (compor(DEFO), '(A16)') defo_comp
    lVari = L_VARI(option)
    lSigm = L_SIGM(option)
    lMatr = L_MATR(option)
    lMatrPred = L_MATR_PRED(option)
    lPred = L_PRED(option)
    l_defo_meca = defo_ldc .eq. 'MECANIQUE'
    l_czm = typmod(2) .eq. 'ELEMJOIN'
    l_large = defo_comp .eq. 'SIMO_MIEHE' .or. defo_comp .eq. 'GROT_GDEP'

    ! En phase de prediction / defo_meca, deps est tel que deps_meca = 0 (strucure additive defos)
    if (l_defo_meca .and. lPred) then

        call behaviourPrepESVAGauss(carcri, defo_ldc, imate, fami, kpg, ksp, neps, instap, BEHinteg)

        epsm_meca = epsm
        deps_meca = 0
        call behaviourPrepStrain(lPred, l_czm, l_large, l_defo_meca, imate, fami, kpg, ksp, &
                                 neps, BEHinteg%esva, epsm_meca, deps_meca)
        deps = -deps_meca
    end if

! --------------------------------------------------------------------------------------------------
!   Integration standard du comportement
! --------------------------------------------------------------------------------------------------

    if (compor(PLANESTRESS) (1:7) .ne. 'DEBORST') then

        call redece(BEHinteg, &
                    fami, kpg, ksp, ndim, typmod, &
                    l_epsi_varc, imate, materi, compor, mult_comp, &
                    carcri, instam, instap, neps, epsm, &
                    deps, nsig, sigm, nvi_all, vim, option, &
                    angmas, numlc, sigp, vip, &
                    ndsde, dsidep, codret)

        if (codret .eq. 1) goto 999

! --------------------------------------------------------------------------------------------------
!  Resolution des contraintes planes sizz=0 par une methode de Newton pour les lois non equipees
! --------------------------------------------------------------------------------------------------
    else

        ! Controles
        ASSERT(ndim .eq. 2)
        ASSERT(nsig .ge. 2*ndim)
        ASSERT(neps .ge. 2*ndim)
        ASSERT(compor(DEFO) .eq. 'PETIT')

        ! Modification des parametres
        nvi = nvi_all-1
        typmod_cp(1) = 'AXIS'
        typmod_cp(2) = typmod(2)

        ! Definition du critere de convergence
        prec = carcri(RESI_DEBORST_MAX)
        ASSERT(prec .ne. r8vide())
        if (prec .ge. 0.d0) then
            typ_crit = 'ABSOLU'
        else
            typ_crit = 'RELATIF'
            prec = -prec
        end if

        ! Estimation de la deformation epzz en utilisant le mecanisme vip_k contient vip_(k-1)
        epsm(3) = vim(nvi+1)
        if (.not. lVari) then
            deps(3) = 0.d0
        else
            deps(3) = vip(nvi+1)-vim(nvi+1)
        end if

        ! S'il faut calculer les contraintes, determination de epzz par methode de Newton
        if (lSigm) then

            ! Creation de l'espace des variables internes si necessaire
            allocate (vip_cp(nvi))
            if (lVari) then
                vip_cp(1:nvi) = vip(1:nvi)
            else
                vip_cp(1:nvi) = vim(1:nvi)
            end if

            do icp = 1, nint(carcri(ITER_DEBORST_MAX))

                ! Choix de l'option pour accéder à la matrice tangente pour methode de Newton
                if (icp .eq. 1 .and. lMatrPred) then
                    option_cp = 'RIGI_MECA_TANG'
                else
                    option_cp = 'FULL_MECA'
                end if

                ! Integration du comportement
                call redece(BEHinteg, &
                            fami, kpg, ksp, ndim, typmod_cp, &
                            l_epsi_varc, imate, materi, compor, mult_comp, &
                            carcri, instam, instap, neps, epsm, &
                            deps, nsig, sigm, nvi, vim, option_cp, &
                            angmas, numlc, sigp, vip_cp, &
                            ndsde, dsidep_cp, codret)

                if (codret .eq. 1) then
                    deallocate (vip_cp)
                    goto 999
                end if

                ! Test de convergence
                if (typ_crit .eq. 'ABSOLU') then
                    conv_cp = abs(sigp(3)) .le. prec
                else
                    conv_cp = abs(sigp(3)) .le. prec*maxval(abs(sigp(1:2*ndim)))
                end if
                if (conv_cp) exit

                ! Reactualisation de la deformation EPZZ en verifiant l'inversibilite
                invert = ASTER_TRUE
                if (abs(dsidep_cp(3, 3)) .lt. abs(sigp(3))) then
                    invert = abs(dsidep_cp(3, 3))/abs(sigp(3)) .gt. 1.d0/r8gaem()
                end if

                if (invert) then
                    deps(3) = deps(3)-sigp(3)/dsidep_cp(3, 3)
                else
                    ! Pivot nul
                    exit
                end if
            end do
            deallocate (vip_cp)
        end if

        ! Integration du comportement avec le bon epzz et l'option reelle
        call redece(BEHinteg, &
                    fami, kpg, ksp, ndim, typmod_cp, &
                    l_epsi_varc, imate, materi, compor, mult_comp, &
                    carcri, instam, instap, neps, epsm, &
                    deps, nsig, sigm, nvi, vim, option, &
                    angmas, numlc, sigp, vip, &
                    ndsde, dsidep, codret)

        if (codret .eq. 1) goto 999

        ! Test de convergence pour le code retour
        if (lSigm) then
            if (typ_crit .eq. 'ABSOLU') then
                codret = merge(0, 3, abs(sigp(3)) .le. prec)
            else
                codret = merge(0, 3, abs(sigp(3)) .le. prec*maxval(abs(sigp(1:2*ndim))))
            end if
        end if

        ! Correction de la matrice tangente pour tenir compte des contraintes planes
        if (lMatr) then

            ! pivot nul -> on ne corrige pas la matrice
            ka3_max = max(maxval(abs(dsidep(1:2, 3))), maxval(abs(dsidep(4:, 3))))
            k3a_max = max(maxval(abs(dsidep(3, 1:2))), maxval(abs(dsidep(3, 4:))))
            c_max = ka3_max*k3a_max
            invert = ASTER_TRUE
            if (abs(dsidep(3, 3)) .lt. c_max) then
                invert = abs(dsidep(3, 3))/c_max .gt. 1.d0/r8gaem()
            end if

            if (invert) then
                do k = 1, nsig
                    if (k .eq. 3) goto 136
                    do l = 1, neps
                        if (l .eq. 3) goto 137
                        dsidep(k, l) = dsidep(k, l)-dsidep(k, 3)*dsidep(3, l)/dsidep(3, 3)
137                     continue
                    end do
136                 continue
                end do
                dsidep(:, 3) = 0
                dsidep(3, :) = 0
            end if
        end if

        ! Actualisation de la deformation epzz dans les variables internes
        if (lVari) then
            vip(nvi+1) = epsm(3)+deps(3)
        end if

    end if

! - Prediction: contribution of the thermal stress to the Taylor expansion if needed
    if (l_defo_meca .and. lPred) then
        if (.not. l_czm) then
            ndimsi = 2*ndim
            ! A remettre suite à la fiche issue32329
            !ASSERT(.not. l_large)
            ASSERT(typmod(2) .eq. ' ' .or. typmod(2) .eq. 'GRADVARI' .or. typmod(2) .eq. 'HHO')
            ASSERT(nsig .ge. ndimsi)
            ASSERT(size(dsidep, 1) .ge. ndimsi)
            ASSERT(size(dsidep, 2) .ge. ndimsi)
            ASSERT(lSigm .and. lMatr)

            call behaviourPredictionStress(BEHinteg%esva, dsidep, sigp)
        end if
    end if

!   Examen du domaine de validité
    call lcvali(fami, kpg, ksp, imate, materi, &
                compor, ndim, epsm, deps, instam, &
                instap, codret)

999 continue
end subroutine
