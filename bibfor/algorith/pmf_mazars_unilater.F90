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

subroutine pmf_mazars_unilater(for_pmf, nf, nbvalc, &
                               compor, crit, defam, defap, varim, &
                               varimp, contm, defm, ddefp, modf, &
                               sigf, varip, codret)
!
! --------------------------------------------------------------------------------------------------
!
!               COMPORTEMENT DES ÉLÉMENTS DE POUTRE MULTI-FIBRES
!
! --------------------------------------------------------------------------------------------------
!
!   IN
!       compor  : information sur le comportement du groupe de fibres
!       crit    : critères de convergence locaux
!       nf      : nombre de fibres du groupe
!       nbvalc  : nombre de variable internes
!       defam   : déformations anélastiques a l'instant précédent
!       defap   : déformations anélastiques a l'instant du calcul
!       varim   : variables internes moins
!       varimp  : variables internes itération précédente (pour DE BORST)
!       contm   : contraintes moins par fibre
!       defm    : déformation à l'instant du calcul précédent
!       ddefp   : incrément de déformation
!
!   OUT
!       modf    : module tangent des fibres
!       sigf    : contrainte a l'instant actuel des fibres
!       varip   : variables internes a l'instant actuel
!       codret  : code retour (0 c'est ok)
!
! --------------------------------------------------------------------------------------------------
!
    use pmfcom_type
!
    implicit none
!
#include "MultiFiber_type.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/mazarsmu.h"
#include "asterfort/rcexistvarc.h"
#include "asterfort/rcvala.h"
#include "asterfort/rcvalb.h"
#include "asterfort/rcvarc.h"
#include "asterfort/utmess.h"
#include "asterfort/verift.h"
!
! --------------------------------------------------------------------------------------------------
!
    type(pmfcom_user), intent(in) :: for_pmf
    integer(kind=8)      :: nf, nbvalc, codret
    real(kind=8) :: contm(nf), defm(nf), ddefp(nf), modf(nf), sigf(nf)
    real(kind=8) :: varimp(nbvalc*nf), varip(nbvalc*nf), varim(nbvalc*nf)
    real(kind=8) :: crit(*), defap(*), defam(*)
!
    character(len=24) :: compor(*)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8)      :: ksp, fib, ivari, iret, nbvari_grfibre
    real(kind=8) :: ep, depsth, tempplus, tempmax
    real(kind=8) :: epsmeca(6), depsmeca(6), dsidep(6, 6), sigf_mu(6)
    real(kind=8) :: nu, bendo, kdess, valsech, valsechref, valhydr
    character(len=4)  :: fami
    character(len=8)  :: materi
    character(len=16) :: rela_comp, algo, nomres(2)
    character(len=30) :: valkm(3)
    aster_logical     :: istemp, ishydr, issech, resi, rigi
! --------------------------------------------------------------------------------------------------
    integer(kind=8), parameter :: nbval = 10
    integer(kind=8)             :: icodre(nbval)
    real(kind=8)        :: valres(nbval)
!   Index des coefficients de la loi de mazars
!       iepsd0 = 1, ik = 2, iac = 3, ibc = 4, iat = 5, ibt = 6
    integer(kind=8), parameter :: iepsd0 = 1
    integer(kind=8), parameter :: isigmlim = 7, iepsilim = 8, iepsc0 = 9, iepst0 = 10
    integer(kind=8), parameter :: iyoung = 11, inu = 12
    character(len=8), parameter :: mazars(nbval) = ['EPSD0   ', 'K       ', 'AC      ', &
                                                    'BC      ', 'AT      ', 'BT      ', &
                                                    'SIGM_LIM', 'EPSI_LIM', 'EPSC0   ', &
                                                    'EPST0   ']
    real(kind=8) :: valmazars(nbval+2)
!
!   Variables internes : ici on a seulement besoin de itemp
!       itemp  : température maximale atteinte par le matériau
    integer(kind=8), parameter :: itemp = 7
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8)             :: kpg
    integer(kind=8)             :: debsp
    integer(kind=8)             :: icdmat
    real(kind=8)        :: instam
    real(kind=8)        :: instap
    real(kind=8)        :: epsm
    character(len=16)   :: option
!
! --------------------------------------------------------------------------------------------------
    kpg = for_pmf%kpg
    icdmat = for_pmf%icdmat
    option = for_pmf%option
    debsp = for_pmf%debsp
    instam = for_pmf%instam
    instap = for_pmf%instap
    epsm = for_pmf%epsm
! --------------------------------------------------------------------------------------------------
    codret = 0
    fami = 'RIGI'
    materi = compor(MULTI_FIBER_MATER) (1:8)
    rela_comp = compor(MULTI_FIBER_RELA) (1:16)
    algo = compor(MULTI_FIBER_ALGO) (1:16)
    read (compor(MULTI_FIBER_NBVARI), '(I24)') nbvari_grfibre
!   Vérification de la loi
    ASSERT(rela_comp .eq. 'MAZARS_UNIL')
!   Vérification du nombre de fibre
    ASSERT(nbvari_grfibre .le. nbvalc)
!   Attention :
!       nbvari_grfibre : nombre de Vint du comportement
!       nbvalc         : nombre de Vint sur les pts de Gauss de la PMF
!                        nbvalc = max( Vint des comportements sur la PMF )
!       ==> Le décalage c'est donc avec "nbvalc"
!
!   Initialisations
    sigf(1:nf) = 0.d0
    ep = 0.0; nu = 0.0; depsth = 0.0
!   Par défaut c'est nul
    bendo = 0.0; kdess = 0.0
!   Valeur des champs : par défaut on considère qu'ils sont nuls
    valhydr = 0.0; valsech = 0.0; valsechref = 0.0; tempmax = 0.0
    !
    ! Température ou pas ?
    istemp = rcexistvarc('TEMP')
    ! Y a-t-il de HYDR ou SECH.
    ishydr = rcexistvarc('HYDR')
    issech = rcexistvarc('SECH')
    !
    if (ishydr .or. issech) then
        nomres(1) = 'B_ENDOGE'; nomres(2) = 'K_DESSIC'; valres(:) = 0.0
        call rcvala(icdmat, ' ', 'ELAS', 0, ' ', [0.0d0], 2, nomres, valres, icodre, 0)
        bendo = valres(1); kdess = valres(2)
        if ((icodre(1) .eq. 0) .and. (.not. ishydr)) then
            valkm(1) = 'MAZARS_UNIL'
            valkm(2) = 'ELAS/B_ENDOGE'
            valkm(3) = 'HYDR'
            call utmess('F', 'COMPOR1_74', nk=3, valk=valkm)
        end if
        if ((icodre(2) .eq. 0) .and. (.not. issech)) then
            valkm(1) = 'MAZARS_UNIL'
            valkm(2) = 'ELAS/K_DESSIC'
            valkm(3) = 'SECH'
            call utmess('F', 'COMPOR1_74', nk=3, valk=valkm)
        end if
    end if
    !
    ! On récupère les paramètres matériau si pas de température ==> ils sont constants
    if (.not. istemp) then
        ksp = debsp
        ! Élasticité
        nomres(1) = 'E'; nomres(2) = 'NU'; valres(:) = 0.0
        call rcvalb(fami, kpg, ksp, '+', icdmat, materi, 'ELAS', &
                    0, ' ', [0.0d0], 2, nomres, valres, icodre, 1)
        ep = valres(1); nu = valres(2)
        ! Mazars
        valres(:) = 0.0
        call rcvalb(fami, kpg, ksp, '+', icdmat, materi, 'MAZARS', &
                    0, ' ', [0.0d0], nbval, mazars, valres, icodre, 0)
        valmazars(1:nbval) = valres(1:nbval)
        if (icodre(isigmlim)+icodre(iepsilim) .ne. 0) then
            valkm(1) = 'MAZARS_UNIL'
            valkm(2) = mazars(isigmlim)
            valkm(3) = mazars(iepsilim)
            call utmess('F', 'COMPOR1_76', nk=3, valk=valkm)
        end if
        ! C'est soit iepsd0 soit iepst0 et iepsc0
        ASSERT(icodre(iepst0) .ne. icodre(iepsd0))
        ASSERT(icodre(iepst0) .eq. icodre(iepsc0))
        if (icodre(iepsd0) .eq. 0) then
            valmazars(iepst0) = valres(iepsd0)
            valmazars(iepsc0) = valres(iepsd0)/(nu*sqrt(2.d0))
        else
            valmazars(iepsd0) = valres(iepst0)
        end if
        ! Ajout de YOUNG, NU dans "valmazars"
        valmazars(iyoung) = ep
        valmazars(inu) = nu
    end if
    ! On mémorise varip(température) que si "resi"
    rigi = (option(1:4) .eq. 'RIGI' .or. option(1:4) .eq. 'FULL')
    resi = (option(1:4) .eq. 'RAPH' .or. option(1:4) .eq. 'FULL')
    ! Boucle sur chaque fibre
    do fib = 1, nf
        ivari = nbvalc*(fib-1)+1
        ksp = debsp-1+fib
        ! On récupère les paramètres matériau si variable de commande.
        ! Elles sont ELGA et peuvent donc être différentes d'un sous point à l'autre.
        if (istemp) then
            call verift(fami, kpg, ksp, '+', icdmat, materi_=materi, &
                        epsth_=depsth, temp_curr_=tempplus)
            ! Température maximale
            tempmax = max(tempplus, varim(ivari-1+itemp))
            ! Mémorise ou pas la température maximale atteinte
            if (resi) then
                varip(ivari-1+itemp) = tempmax
            end if
            ! Élasticité fct de la température max
            nomres(1) = 'E'; nomres(2) = 'NU'; valres(:) = 0.0
            call rcvalb(fami, kpg, ksp, '+', icdmat, materi, 'ELAS', &
                        1, 'TEMP', [tempmax], 2, nomres, valres, icodre, 1)
            ep = valres(1); nu = valres(2)
            ! Mazars fct de la température max
            valres(:) = 0.0
            call rcvalb(fami, kpg, ksp, '+', icdmat, materi, 'MAZARS', &
                        1, 'TEMP', [tempmax], nbval, mazars, valres, icodre, 0)
            valmazars(1:nbval) = valres(1:nbval)
            ! C'est soit iepsd0 soit iepst0 et iepsc0
            if (icodre(iepsd0) .eq. 0) then
                valmazars(iepst0) = valres(iepsd0)
                valmazars(iepsc0) = valres(iepsd0)/(nu*sqrt(2.d0))
            else
                valmazars(iepsd0) = valres(iepst0)
            end if
            ! Ajout de YOUNG, NU dans "valmazars"
            valmazars(iyoung) = ep
            valmazars(inu) = nu
        end if
        if (ishydr) then
            call rcvarc('F', 'HYDR', '+', fami, kpg, ksp, valhydr, iret)
        end if
        if (issech) then
            call rcvarc('F', 'SECH', '+', fami, kpg, ksp, valsech, iret)
            call rcvarc('F', 'SECH', 'REF', fami, kpg, ksp, valsechref, iret)
        end if
        !
        ! Initialisations
        epsmeca(:) = 0.0; depsmeca(:) = 0.0; dsidep(:, :) = 0.0
        !
        epsmeca(1) = defm(fib)-depsth-kdess*(valsech-valsechref)+bendo*valhydr
        depsmeca(1) = ddefp(fib)
        call mazarsmu(option, epsmeca, depsmeca, 1, valmazars, &
                      varim(ivari), varip(ivari), sigf_mu, dsidep)
        modf(fib) = dsidep(1, 1)
        sigf(fib) = sigf_mu(1)
    end do
!
end subroutine
