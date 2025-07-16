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

subroutine pmfcom(kpg, debsp, option, compor, crit, &
                  nf, instam, instap, icdmat, nbvalc, &
                  defam, defap, varim, varimp, contm, &
                  defm, ddefp, epsm, modf, sigf, &
                  varip, codret)
!
! aslint: disable=W1504
! --------------------------------------------------------------------------------------------------
!
!               COMPORTEMENT DES ÉLÉMENTS DE POUTRE MULTI-FIBRES
!
! --------------------------------------------------------------------------------------------------
!
!   IN
!       kpg     : numéro de point de gauss
!       debsp   : numéro de sous-point de la première fibre du groupe
!       option  : option de calcul
!       compor  : information sur le comportement du groupe de fibres
!       crit    : critères de convergence locaux
!       nf      : nombre de fibres du groupe
!       instam  : instant du calcul précédent
!       instap  : instant du calcul
!       icdmat  : code matériau
!       nbvalc  : nombre de variable internes
!       defam   : déformations anélastiques a l'instant précédent
!       defap   : déformations anélastiques a l'instant du calcul
!       varim   : variables internes moins
!       varimp  : variables internes itération précédente (pour DE BORST)
!       contm   : contraintes moins par fibre
!       defm    : déformation  a l'instant du calcul precedent
!       ddefp   : incrément de déformation
!       epsm    : déformation a l'instant précédent sur l'élément de structure
!
!   OUT
!       modf    : module tangent des fibres
!       sigf    : contrainte a l'instant actuel des fibres
!       varip   : variables internes a l'instant actuel
!       codret :
!
! --------------------------------------------------------------------------------------------------
!
    use pmfcom_type
    implicit none
!
#include "MultiFiber_type.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/pmf_mazars_unilater.h"
#include "asterfort/pmf_pinto_menegotto.h"
#include "asterfort/pmf_vmis.h"
#include "asterfort/nm1dco.h"
#include "asterfort/nm1vil.h"
#include "asterfort/paeldt.h"
#include "asterfort/rcexistvarc.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utmess.h"
!
#include "jeveux.h"
!
!
    integer(kind=8)      :: nf, icdmat, nbvalc, kpg, debsp, codret
    real(kind=8) :: contm(nf), defm(nf), ddefp(nf), modf(nf), sigf(nf)
    real(kind=8) :: varimp(nbvalc*nf), varip(nbvalc*nf), varim(nbvalc*nf)
    real(kind=8) :: instam, instap, epsm
    real(kind=8) :: crit(*), defap(*), defam(*)
!
    character(len=16) :: option
    character(len=24) :: compor(*)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbval = 1
    integer(kind=8)    :: icodre(nbval)
    real(kind=8)       :: valres(nbval)

    integer(kind=8)    ::  codrep, ksp, fib, ivari, nbvari_grfibre
    real(kind=8) :: ep, em, depsth, tref, tempm, tempp
    real(kind=8) :: angmas(3), depsm, nu
!
    character(len=4)  :: fami
    character(len=8)  :: materi
    character(len=16) :: rela_comp, algo
    character(len=30) :: valkm(3)
!
    aster_logical     :: istemp
!
    type(pmfcom_user) :: for_pmf
!
! --------------------------------------------------------------------------------------------------
    codret = 0
    codrep = 0
    fami = 'RIGI'
    materi = compor(MULTI_FIBER_MATER) (1:8)
    rela_comp = compor(MULTI_FIBER_RELA) (1:16)
    algo = compor(MULTI_FIBER_ALGO) (1:16)
    read (compor(MULTI_FIBER_NBVARI), '(I24)') nbvari_grfibre
!   Vérification du nombre de fibre
    ASSERT(nbvari_grfibre .le. nbvalc)
!   Attention :
!       nbvari_grfibre : nombre de Vint du comportement
!       nbvalc         : nombre de Vint sur les pts de Gauss de la PMF
!                        nbvalc = max( Vint des comportements sur la PMF )
!       ==> Le décalage c'est donc avec "nbvalc" et pas le nb de Vint du comportement
!
!   Initialisation
    sigf(1:nf) = 0.d0
    ! Angle du MOT_CLEF massif (AFFE_CARA_ELEM) initialise à 0.0 (on ne s'en sert pas)
    angmas(:) = 0.0
!
! --------------------------------------------------------------------------------------------------
    if (rela_comp .eq. 'ELAS') then
        ! Température ou pas ?
        istemp = rcexistvarc('TEMP')
        if (.not. istemp) then
            call rcvalb(fami, 1, 1, '+', icdmat, materi, 'ELAS', &
                        0, '', [0.d0], 1, ['E'], valres, icodre, 1)
            ep = valres(1)
            em = ep
            ! Boucle sur chaque fibre
            do fib = 1, nf
                ksp = debsp-1+fib
                modf(fib) = ep
                sigf(fib) = ep*(contm(fib)/em+ddefp(fib))
            end do
        else
            ! Boucle sur chaque fibre avec température
            do fib = 1, nf
                ksp = debsp-1+fib
                call paeldt(kpg, ksp, fami, 'T', icdmat, materi, em, ep, nu, depsth)
                modf(fib) = ep
                sigf(fib) = ep*(contm(fib)/em+ddefp(fib)-depsth)
            end do
        end if
!
! --------------------------------------------------------------------------------------------------
    else if (rela_comp .eq. 'MAZARS_UNIL') then
        for_pmf%kpg = kpg
        for_pmf%icdmat = icdmat
        for_pmf%option = option
        for_pmf%debsp = debsp
        for_pmf%instam = instam
        for_pmf%instap = instap
        for_pmf%epsm = epsm
        call pmf_mazars_unilater(for_pmf, nf, nbvalc, &
                                 compor, crit, defam, defap, varim, &
                                 varimp, contm, defm, ddefp, modf, &
                                 sigf, varip, codret)
!
! --------------------------------------------------------------------------------------------------
    else if (rela_comp .eq. 'PINTO_MENEGOTTO') then
        for_pmf%kpg = kpg
        for_pmf%icdmat = icdmat
        for_pmf%option = option
        for_pmf%debsp = debsp
        for_pmf%instam = instam
        for_pmf%instap = instap
        for_pmf%epsm = epsm
        call pmf_pinto_menegotto(for_pmf, nf, nbvalc, &
                                 compor, crit, defam, defap, varim, &
                                 varimp, contm, defm, ddefp, modf, &
                                 sigf, varip, codret)
!
! --------------------------------------------------------------------------------------------------
    else if ((rela_comp .eq. 'VMIS_CINE_GC') .or. &
             (rela_comp .eq. 'VMIS_CINE_LINE') .or. &
             (rela_comp .eq. 'VMIS_ISOT_LINE') .or. &
             (rela_comp .eq. 'VMIS_ISOT_TRAC')) then
        for_pmf%kpg = kpg
        for_pmf%icdmat = icdmat
        for_pmf%option = option
        for_pmf%debsp = debsp
        for_pmf%instam = instam
        for_pmf%instap = instap
        for_pmf%epsm = epsm
        call pmf_vmis(for_pmf, nf, nbvalc, &
                      compor, crit, defam, defap, varim, &
                      varimp, contm, defm, ddefp, modf, &
                      sigf, varip, codret)
!
! --------------------------------------------------------------------------------------------------
    else if (rela_comp .eq. 'CORR_ACIER') then
        ! Température ou pas ?
        istemp = rcexistvarc('TEMP')
        if (.not. istemp) then
            call rcvalb(fami, 1, 1, '+', icdmat, materi, 'ELAS', &
                        0, '', [0.d0], 1, ['E'], valres, icodre, 1)
            ep = valres(1)
            ! Boucle sur chaque fibre
            do fib = 1, nf
                ivari = nbvalc*(fib-1)+1
                ksp = debsp-1+fib
                depsm = ddefp(fib)
                call nm1dco('RIGI', kpg, fib, option, icdmat, &
                            materi, ep, contm(fib), defm(fib), depsm, &
                            varim(ivari), sigf(fib), varip(ivari), modf(fib), crit, &
                            codret)
                if (codret .ne. 0) goto 999
            end do
        else
            ! Boucle sur chaque fibre
            do fib = 1, nf
                ivari = nbvalc*(fib-1)+1
                ksp = debsp-1+fib
                call paeldt(kpg, ksp, fami, '+', icdmat, materi, em, ep, nu, depsth)
                depsm = ddefp(fib)-depsth
                call nm1dco('RIGI', kpg, fib, option, icdmat, &
                            materi, ep, contm(fib), defm(fib), depsm, &
                            varim(ivari), sigf(fib), varip(ivari), modf(fib), crit, &
                            codret)
                if (codret .ne. 0) goto 999
            end do
        end if
!
! --------------------------------------------------------------------------------------------------
    else if ((rela_comp .eq. 'GRAN_IRRA_LOG') .or. &
             (rela_comp .eq. 'VISC_IRRA_LOG')) then
        ! C'est le seul algo disponible dans les catalogues de ces 2 comportements
        if (algo(1:10) .ne. 'ANALYTIQUE') then
            valkm(1) = rela_comp
            valkm(2) = 'DEFI_COMPOR/MULTIFIBRE'
            valkm(3) = algo(1:10)
            call utmess('F', 'COMPOR5_81', nk=3, valk=valkm)
        end if
        ! Température ou pas ?
        istemp = rcexistvarc('TEMP')
        if (.not. istemp) then
            call utmess('F', 'COMPOR5_40', sk=rela_comp)
        end if
        ! Boucle sur chaque fibre
        do fib = 1, nf
            ivari = nbvalc*(fib-1)+1
            ksp = debsp-1+fib
            call paeldt(kpg, ksp, fami, 'T', icdmat, materi, em, ep, nu, depsth, &
                        tmoins=tempm, tplus=tempp, trefer=tref)
            depsm = ddefp(fib)-depsth
            call nm1vil('RIGI', kpg, fib, icdmat, materi, &
                        crit, instam, instap, tempm, tempp, &
                        tref, depsm, contm(fib), varim(ivari), option, &
                        defam(1), defap(1), angmas, sigf(fib), varip(ivari), &
                        modf(fib), codret, rela_comp, nbvalc)
            if (codret .ne. 0) goto 999
        end do
!
! --------------------------------------------------------------------------------------------------
    else
        call utmess('F', 'ELEMENTS2_39', sk=rela_comp)
    end if
!
999 continue
end subroutine
