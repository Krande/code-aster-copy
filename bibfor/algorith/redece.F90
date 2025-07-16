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
! aslint: disable=W1306,W1504,C1505
!
subroutine redece(BEHinteg, &
                  fami, kpg, ksp, ndim, typmod, &
                  l_epsi_varc, imate, materi, compor, mult_comp, &
                  carcri, instam, instap, neps, epsm, &
                  deps, nsig, sigm, nvi, vim, option, &
                  angmas, numlc, sigp, vip, &
                  ndsde, dsidep, codret)

    use calcul_module, only: ca_iredec_, ca_td1_, ca_tf1_, ca_timed1_, ca_timef1_
    use Behaviour_type
    use Behaviour_module

    implicit none

#include "asterf_types.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/assert.h"
#include "asterfort/lc0000.h"
#include "asterfort/utmess.h"

    type(Behaviour_Integ) :: BEHinteg
    character(len=*) :: fami
    integer(kind=8) :: imate, ndim, kpg, ksp, numlc
    integer(kind=8) :: neps, nsig, ndsde
    aster_logical, intent(in) :: l_epsi_varc
    integer(kind=8), intent(in):: nvi
    real(kind=8) :: carcri(*), angmas(*)
    real(kind=8) :: instam, instap
    real(kind=8) :: epsm(neps), deps(neps)
    real(kind=8) :: sigm(nsig), sigp(nsig)
    real(kind=8) :: vim(nvi), vip(nvi)
    real(kind=8) :: dsidep(merge(nsig, 6, nsig*neps .eq. ndsde), &
                           merge(neps, 6, nsig*neps .eq. ndsde))
    character(len=8)  :: typmod(*)
    character(len=16) :: compor(*), option
    character(len=8), intent(in) :: materi
    character(len=16), intent(in) :: mult_comp
    integer(kind=8), intent(out):: codret

! --------------------------------------------------------------------------------------------------
! INTEGRATION DU COMPORTEMENT : PRISE EN CHARGE D'UN EVENTUEL REDECOUPAGE LOCAL DU PAS DE TEMPS
! --------------------------------------------------------------------------------------------------
! in  fami,kpg,ksp  : famille et numero du (sous)point de gauss
!     ndim    : dimension de l'espace
!               3 : 3d , 2 : d_plan ,axis ou  c_plan
!     typmod  : modelisation ex: 1:3d, 2:inco
!     imate   : adresse du materiau code
!     compor  : comportement :  (1) = type de relation comportement
!                               (2) = nb variables internes / pg
!                               (3) = hypothese sur les deformations
!                               (4) etc... (voir grandeur compor)
!     mult_comp : multi-comportement (pour polycristal)
!     carcri  : criteres de convergence locaux (voir grandeur carcri)
!     instam  : instant du calcul precedent
!     instap  : instant du calcul
!     neps    : nombre de cmp de epsdt et depst (suivant modelisation)
!     epsm   : deformations totales a l'instant du calcul precedent
!     deps   : increment de deformation totale :
!                depst(t) = depst(mecanique(t)) + depst(dilatation(t))
!     nsig    : nombre de cmp de sigd et sigf (suivant modelisation)
!     sigm    : contraintes au debut du pas de temps
!     vim    : variables internes au debut du pas de temps
!     option  : option demandee : rigi_meca_tang , full_meca , raph_meca
!     angmas  : les trois angles du mot_clef massif (affe_cara_elem),
!               + un reel qui vaut 0 si nautiquies ou 2 si euler
!               + les 3 angles d'euler
!     numlc   : numero de loi de comportement issue du catalogue de lc
!
! out sigp    : contraintes a la fin du pas de temps
! var vip     : variables internes
!                in  : estimation (iteration precedente ou lag. augm.)
!                out : a la fin du pas de temps t+
!     ndsde   : dimension de dsidep
!     dsidep    : operateur tangent dsig/deps ou dsig/df
!     codret  : code retour loi de comporment :
!               codret=0 : tout va bien
!               codret=1 : echec dans l'integration de la loi
!               codret=2 : ok loi de comportement mais condition non respectee -> pas de cvg globale
! --------------------------------------------------------------------------------------------------
    aster_logical:: lMatrPred, lMatr, lSigm, lVari
    integer(kind=8), parameter:: SANS = 0, FORCE = 1, AUTO = 2, NBR_DECOUP_MAX = 5
    integer(kind=8):: decoup
    integer(kind=8) :: niv_dec, niv_ini, npas, iip, codret_sub, pas
    real(kind=8) :: epsm_sub(neps), deps_sub(neps), sigm_sub(nsig), vim_sub(nvi)
    real(kind=8) :: deltat, tm, tp
    real(kind=8) :: dsidep_sub(merge(nsig, 6, nsig*neps .eq. ndsde), &
                               merge(neps, 6, nsig*neps .eq. ndsde))
    character(len=16)::defo_comp
! --------------------------------------------------------------------------------------------------

    iip = nint(carcri(ITER_INTE_PAS))
    lSigm = L_SIGM(option)
    lMatr = L_MATR(option)
    lVari = L_VARI(option)
    lMatrPred = L_MATR_PRED(option)
    read (compor(DEFO), '(A16)') defo_comp

    ! Pour les variables de commande
    ca_iredec_ = 1
    ca_timed1_ = instam
    ca_timef1_ = instap
    ca_td1_ = instam
    ca_tf1_ = instap

    ! Strategie de redecoupage local
    if (abs(iip) .le. 1 .or. lMatrPred) then
        decoup = SANS
    else if (iip .le. -2) then
        decoup = AUTO
    else if (iip .ge. 2) then
        decoup = FORCE
    end if

! --------------------------------------------------------------------------------------------------
!  Integration du comportement sans redecoupage
! --------------------------------------------------------------------------------------------------

    if (decoup .eq. SANS) then
        niv_dec = 2
        codret = 0
        call lc0000(BEHinteg, &
                    fami, kpg, ksp, ndim, typmod, &
                    l_epsi_varc, imate, materi, compor, mult_comp, &
                    carcri, instam, instap, neps, epsm, &
                    deps, nsig, sigm, vim, option, &
                    angmas, numlc, sigp, vip, &
                    ndsde, dsidep, niv_dec, nvi, codret)
        goto 999
    end if

! --------------------------------------------------------------------------------------------------
!  Integration du comportement avec redecoupage
! --------------------------------------------------------------------------------------------------

    ! Tests de compatibilite
    if (typmod(2) .eq. 'GRADVARI') call utmess('F', 'COMPOR2_10', sk=typmod(2))
    if (numlc .ge. 8000 .and. numlc .lt. 9000) call utmess('F', 'COMPOR2_10', sk='KIT_DDI')
    if (defo_comp .eq. 'SIMO_MIEHE') call utmess('F', 'COMPOR2_10', sk='SIMO_MIEHE')
    ASSERT(lSigm)
    ASSERT(lVari)

    ! niveau de decoupe initial (si FORCE, on commence avec npas, si AUTO avec 1 seul pas)
    niv_ini = merge(0, 1, decoup .eq. AUTO)

    ! Boucle sur les niveaux de decoupage
    do niv_dec = niv_ini, NBR_DECOUP_MAX
        codret = 0
        npas = max(1, merge(1, abs(iip)*(2**(niv_dec-1)), niv_dec .eq. 0))
        deltat = (instap-instam)/npas
        deps_sub = deps/npas

        ! Initialisation des parametres d'entree
        vim_sub = vim
        sigm_sub = sigm

        if (lMatr) dsidep = 0

        ! Boucle sur les sous-pas
        do pas = 1, npas

            epsm_sub = epsm+(pas-1)*deps_sub
            tm = instam+(pas-1)*deltat
            tp = instam+pas*deltat
            ca_td1_ = tm
            ca_tf1_ = tp
            codret_sub = 0
            vip = vim_sub

            call lc0000(BEHinteg, &
                        fami, kpg, ksp, ndim, typmod, &
                        l_epsi_varc, imate, materi, compor, mult_comp, &
                        carcri, tm, tp, neps, epsm_sub, &
                        deps_sub, nsig, sigm_sub, vim_sub, option, &
                        angmas, numlc, sigp, vip, &
                        ndsde, dsidep_sub, niv_dec, nvi, codret_sub)

            select case (codret_sub)
            case (0)
                continue
            case (1)
                codret = 1
                exit
            case (2)
                codret = 2
                if (pas .ne. npas) exit
            case default
                ASSERT(ASTER_FALSE)
            end select

            ! Matrice tangente finale = moyenne des matrices tangentes a chaque pas (faute de mieux)
            if (lMatr) dsidep = dsidep+dsidep_sub/npas

            ! Actualisation des parametres d'entree
            sigm_sub = sigp
            vim_sub = vip

        end do
        ! si codret=0, on a converge ; si codret=1 ou 2, on redecoupe
        if (codret .eq. 0) goto 999

    end do

999 continue
    ASSERT(codret .ge. 0 .and. codret .le. 2)
end subroutine
