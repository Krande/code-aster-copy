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

subroutine pmfcom(kpg,   debsp,  option, compor, crit,   &
                  nf,    instam, instap, icdmat, nbvalc, &
                  defam, defap,  varim,  varimp, contm,  &
                  defm,  ddefp,  epsm,   modf,   sigf,   &
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
!       epsm    : déformation a l'instant précédent
!
!   OUT
!       modf    : module tangent des fibres
!       sigf    : contrainte a l'instant actuel des fibres
!       varip   : variables internes a l'instant actuel
!       codret :
!
! --------------------------------------------------------------------------------------------------
!
!
implicit none

#include "MultiFiber_type.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/comp1d.h"
#include "asterfort/compor_3d_fibre.h"
#include "asterfort/mazu1d.h"
#include "asterfort/nm1dci.h"
#include "asterfort/nm1dco.h"
#include "asterfort/nm1dis.h"
#include "asterfort/nm1dpm.h"
#include "asterfort/nm1vil.h"
#include "asterfort/paeldt.h"
#include "asterfort/rcexistvarc.h"
#include "asterfort/rcvala.h"
#include "asterfort/rcvalb.h"
#include "asterfort/rcvarc.h"
#include "asterfort/utmess.h"
#include "asterfort/verift.h"
#include "asterfort/vmci1d.h"
#include "blas/dcopy.h"
!
#include "jeveux.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/jevech.h"
integer :: icompo
!
    integer :: nf, icdmat, nbvalc, kpg, debsp, codret
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
    integer , parameter :: nbval=12
    integer :: icodre(nbval)
    real(kind=8) :: valres(nbval)

    integer :: nbvari, codrep, ksp, fib, ivari, iret, nbvari_grfibre
    real(kind=8) :: ep, em, depsth, tref, tempm, tempp, sigx, epsx, depsx, tempplus, tempmax
    real(kind=8) :: cstpm(13), angmas(3), depsm, nu, bendo, kdess, valsech, valsechref, valhydr
    character(len=4) :: fami
    character(len=8) :: nompim(12), mazars(8), materi
    character(len=16) :: rela_comp, algo, nomres(2)
    character(len=30) :: valkm(3)
    aster_logical :: istemp, ishydr, issech, resi
!
    data nompim /'SY', 'EPSI_ULT', 'SIGM_ULT', 'EPSP_HAR', 'R_PM',&
                 'EP_SUR_E', 'A1_PM', 'A2_PM', 'ELAN', 'A6_PM', 'C_PM', 'A_PM'/
    data mazars /'EPSD0', 'K', 'AC', 'BC', 'AT', 'BT', 'SIGM_LIM', 'EPSI_LIM'/
!
! --------------------------------------------------------------------------------------------------
    codret = 0
    codrep = 0
    fami   = 'RIGI'
    materi      = compor(MULTI_FIBER_MATER)(1:8)
    rela_comp   = compor(MULTI_FIBER_RELA)(1:16)
    algo        = compor(MULTI_FIBER_ALGO)(1:16)
    read(compor(MULTI_FIBER_NBVARI),'(I24)') nbvari_grfibre
!   Vérification du nombre de fibre
    ASSERT( nbvari_grfibre .le. nbvalc )
!   Initialisation
    sigf(:) = 0.d0
!
!   Température ou pas ?
    istemp = rcexistvarc('TEMP')
!
    if (.not.istemp) then
        nomres(1) = 'E'
        nomres(2) = 'NU'
        call rcvalb(fami, 1, 1, '+', icdmat, materi, 'ELAS', &
                    0, '', [0.d0], 2, nomres, valres, icodre, 1)
        ep = valres(1)
        nu = valres(2)
        em=ep
        depsth=0.d0
    endif
!   Angle du MOT_CLEF massif (AFFE_CARA_ELEM) initialise à 0.0 (on ne s'en sert pas)
    angmas(:) = 0.0
!
! --------------------------------------------------------------------------------------------------
    if (rela_comp .eq. 'ELAS') then
        nomres(1) = 'E'
!       Boucle sur chaque fibre
        do fib = 1, nf
            if (istemp) then
                ksp=debsp-1+fib
                call paeldt(kpg, ksp, fami, 'T', icdmat, materi, em, ep, nu, depsth)
            endif
            modf(fib) = ep
            sigf(fib) = ep*(contm(fib)/em + ddefp(fib) - depsth)
        enddo
!
! --------------------------------------------------------------------------------------------------
    else if (rela_comp.eq.'MAZARS_GC') then
!       Y a-t-il de HYDR ou SECH
!       Par défaut c'est nul
        bendo = 0.0
        kdess = 0.0
!       Valeur des champs : par défaut on considère qu'ils sont nuls
        valhydr    = 0.0
        valsech    = 0.0
        valsechref = 0.0
!
        ishydr = rcexistvarc('HYDR')
        issech = rcexistvarc('SECH')
        if ( ishydr .or. issech ) then
            nomres(1)='B_ENDOGE'
            nomres(2)='K_DESSIC'
            valres(1:2) = 0.0
            call rcvala(icdmat, ' ', 'ELAS', 0, ' ', [0.d0], 2, nomres, valres, icodre,0)
            bendo = valres(1)
            kdess = valres(2)
            if ((icodre(1).eq.0).and.(.not.ishydr)) then
                valkm(1)='MAZARS_GC'
                valkm(2)='ELAS/B_ENDOGE'
                valkm(3)='HYDR'
                call utmess('F', 'COMPOR1_74', nk=3, valk=valkm)
            endif
            if ((icodre(2).eq.0).and.(.not.issech)) then
                valkm(1)='MAZARS_GC'
                valkm(2)='ELAS/K_DESSIC'
                valkm(3)='SECH'
                call utmess('F', 'COMPOR1_74', nk=3, valk=valkm)
            endif
        endif
!
!       On récupère les paramètres matériau si pas de variable de commande ==> ils sont constants
        valres(:) = 0.0
        call rcvalb(fami, 1, 1, '+', icdmat, materi, 'MAZARS', &
                    0, ' ', [0.0d0], 8, mazars, valres, icodre, 1)
        if (icodre(7)+icodre(8) .ne. 0) then
            valkm(1)='MAZARS_GC'
            valkm(2)=mazars(7)
            valkm(3)=mazars(8)
            call utmess('F', 'COMPOR1_76', nk=3, valk=valkm)
        endif
!       On mémorise varip(température) que si "resi"
        resi = (option(1:4).eq.'RAPH' .or. option(1:4).eq.'FULL')
!       Boucle sur chaque fibre
        do fib = 1, nf
            ivari = nbvalc*(fib-1) + 1
            ksp=debsp-1+fib
            if (istemp) then
                call verift(fami, kpg, ksp, '+', icdmat, materi, epsth_=depsth, temp_curr_=tempplus)
                ! Température maximale
                if ( tempplus .gt. varim(ivari+7-1) ) then
                    tempmax = tempplus
                else
                    tempmax = varim(ivari+7-1)
                endif
!               Mémorise ou pas la température maximale atteinte
                if ( resi ) then
                    varip(ivari+7-1) = tempmax
                endif
                nomres(1) = 'E'
                nomres(2) = 'NU'
                call rcvalb(fami, kpg, ksp, '+', icdmat, materi, 'ELAS', &
                            1, 'TEMP', [tempmax], 2, nomres, valres, icodre, 1)
                ep = valres(1)
                nu = valres(2)
            endif
            if ( ishydr ) then
                call rcvarc('F', 'HYDR', '+',   fami, kpg, ksp, valhydr, iret)
            endif
            if ( issech ) then
                call rcvarc('F', 'SECH', '+',   fami, kpg, ksp, valsech, iret)
                call rcvarc('F', 'SECH', 'REF', fami, kpg, ksp, valsechref, iret)
            endif
            epsm = defm(fib) - depsth - kdess*(valsech-valsechref) + bendo*valhydr
!           On récupère les paramètres matériau s'il y a une variable de commande.
!           Elles sont ELGA et peuvent donc être différentes d'un sous point à l'autre.
            if ( istemp .or. ishydr .or. issech ) then
                if ( istemp ) then
                    call rcvalb(fami, kpg, ksp, '+', icdmat, materi, 'MAZARS', &
                                1, 'TEMP', [tempmax], 8, mazars, valres, icodre, 1)
                else
                    call rcvalb(fami, kpg, ksp, '+', icdmat, materi, 'MAZARS', &
                                0, ' ', [0.0d0], 8, mazars, valres, icodre, 1)
                endif
            endif
!           Ajout de NU dans VALRES
            valres(9) = nu
            call mazu1d(ep,         valres,    contm(fib), varim(ivari), epsm, &
                        ddefp(fib), modf(fib), sigf(fib),  varip(ivari), option)
        enddo
!
! --------------------------------------------------------------------------------------------------
    else if (rela_comp.eq.'VMIS_CINE_GC') then
!       Boucle sur chaque fibre
        do fib = 1, nf
            ivari = nbvalc*(fib-1) + 1
            ksp=debsp-1+fib
            if (istemp) then
                call paeldt(kpg, ksp, fami, 'T', icdmat, materi, em, ep, nu, depsth)
            endif
            depsm = ddefp(fib)-depsth
            call vmci1d('RIGI', kpg,        ksp,          icdmat,       em,     &
                        ep,     contm(fib), depsm,        varim(ivari), option, &
                        materi, sigf(fib),  varip(ivari), modf(fib))
        enddo
!
! --------------------------------------------------------------------------------------------------
    else if (rela_comp.eq.'PINTO_MENEGOTTO') then
!       on récupère les paramètres matériau
        valres(:) = 0.0
        call rcvalb(fami, 1, 1, '-', icdmat, materi, 'PINTO_MENEGOTTO', &
                    0, ' ', [0.0d0], 12, nompim, valres, icodre, 0)
        if (icodre(7) .ne. 0) valres(7) = -1.0d0
        cstpm(1) = ep
        do fib = 1, 12
            cstpm(fib+1) = valres(fib)
        enddo
!       Boucle sur chaque fibre
        do fib = 1, nf
            ivari = nbvalc*(fib-1) + 1
            if (istemp) then
                ksp=debsp-1+fib
                call paeldt(kpg, ksp, fami, 'T', icdmat, materi, em, ep, nu, depsth)
                cstpm(1) = ep
            endif
            depsm = ddefp(fib)-depsth
            call nm1dpm('RIGI', kpg,          fib,       icdmat,     option,       &
                        nbvalc, 13,           cstpm,     contm(fib), varim(ivari), &
                        depsm,  varip(ivari), sigf(fib), modf(fib))
        enddo
!
! --------------------------------------------------------------------------------------------------
    else if (rela_comp.eq.'VMIS_CINE_LINE') then
!       Boucle sur chaque fibre
        do fib = 1, nf
            ivari = nbvalc*(fib-1) + 1
            if (istemp) then
                ksp=debsp-1+fib
                call paeldt(kpg, ksp, fami, 'T', icdmat, materi, em, ep, nu, depsth)
            endif
            depsm = ddefp(fib)-depsth
            call nm1dci('RIGI', kpg,        fib,          icdmat,       em,     &
                        ep,     contm(fib), depsm,        varim(ivari), option, &
                        materi, sigf(fib),  varip(ivari), modf(fib))
        enddo
!
! --------------------------------------------------------------------------------------------------
    else if ( (rela_comp.eq.'VMIS_ISOT_LINE') .or. &
              (rela_comp.eq.'VMIS_ISOT_TRAC') ) then
!       Boucle sur chaque fibre
        do fib = 1, nf
            ivari = nbvalc*(fib-1) + 1
            if (istemp) then
                ksp=debsp-1+fib
                call paeldt(kpg, ksp, fami, 'T', icdmat, materi, em, ep, nu, depsth)
            endif
            depsm = ddefp(fib)-depsth
            call nm1dis('RIGI',    kpg,        fib,       icdmat,       em,     &
                        ep,        contm(fib), depsm,     varim(ivari), option, &
                        rela_comp, materi,     sigf(fib), varip(ivari), modf(fib))
        enddo
!
! --------------------------------------------------------------------------------------------------
    else if (rela_comp.eq.'CORR_ACIER') then
!       Boucle sur chaque fibre
        do fib = 1, nf
            ivari = nbvalc*(fib-1) + 1
            if (istemp) then
                ksp=debsp-1+fib
                call paeldt(kpg, ksp, fami, '+', icdmat, materi, em, ep, nu, depsth)
            endif
            depsm = ddefp(fib)-depsth
            call nm1dco('RIGI',       kpg,       fib,          option,    icdmat, &
                        materi,       ep,        contm(fib),   defm(fib), depsm,  &
                        varim(ivari), sigf(fib), varip(ivari), modf(fib), crit,   &
                        codret)
            if (codret .ne. 0) goto 999
        enddo
!
! --------------------------------------------------------------------------------------------------
    else if ( (rela_comp.eq.'GRAN_IRRA_LOG') .or. &
              (rela_comp.eq.'VISC_IRRA_LOG') ) then
        ! C'est le seul algo disponible dans les catalogues de ces 2 comportements
        if (algo(1:10).ne.'ANALYTIQUE') then
            valkm(1) = rela_comp
            valkm(2) = 'DEFI_COMPOR/MULTIFIBRE'
            valkm(3) = algo(1:10)
            call utmess('F', 'COMPOR5_81', nk=3, valk=valkm)
        endif
!
        if (.not. istemp) then
            call utmess('F', 'COMPOR5_40',sk=rela_comp)
        endif
!       Boucle sur chaque fibre
        do fib = 1, nf
            ivari = nbvalc*(fib-1) + 1
            if (istemp) then
                ksp=debsp-1+fib
                call paeldt(kpg, ksp, fami, 'T', icdmat, materi, em, ep, nu, depsth, &
                            tmoins=tempm, tplus=tempp, trefer=tref)
            endif
            depsm = ddefp(fib)-depsth
            call nm1vil('RIGI',    kpg,      fib,        icdmat,       materi,       &
                        crit,      instam,   instap,     tempm,        tempp,        &
                        tref,      depsm,    contm(fib), varim(ivari), option,       &
                        defam(1),  defap(1), angmas,     sigf(fib),    varip(ivari), &
                        modf(fib), codret,   rela_comp,  nbvalc)
            if (codret .ne. 0) goto 999
        enddo
!
! --------------------------------------------------------------------------------------------------
    else if (rela_comp.eq.'BETON_GRANGER') then
        ! Algo de type double 'DE BORST'
        if (algo(1:7).ne.'DEBORST') then
            valkm(1) = rela_comp
            valkm(2) = 'DEFI_COMPOR/MULTIFIBRE'
            valkm(3) = algo(1:7)
            call utmess('F', 'COMPOR5_81', nk=3, valk=valkm)
        endif
        ! La LDC doit retourner :
        !   - sigf(fib)   : la contrainte sur la fibre
        !   - modf(fib)   : le module tangent de la fibre
        !   - varip(fib)  : les variables internes de la LdC sur la fibre
        !   - codrep    : code retour
        if ((option(1:9).eq.'FULL_MECA') .or. (option(1:9) .eq.'RAPH_MECA')) then
            nbvari = nbvalc*nf
            call dcopy(nbvari, varimp, 1, varip, 1)
        endif
        ! A faire en avant la boucle sur les fibres : elles ont toutes le même comportement
        call jevech('PCOMPOR', 'L', icompo)
!       Boucle sur chaque fibre
        do fib = 1, nf
            ivari = nbvalc*(fib-1) + 1
            sigx  = contm(fib)
            epsx  = defm(fib)
            depsx = ddefp(fib)
            call compor_3d_fibre('RIGI',       kpg,       fib,       option, sigx,         &
                                 instam,       instap,    crit,      icdmat, materi,       &
                                 zk16(icompo), epsx,      depsx,     angmas, varim(ivari), &
                                 varip(ivari), sigf(fib), modf(fib), codrep)
            if (codrep .ne. 0) then
                codret=codrep
                ! code 3 : on continue et on le renvoie à la fin. Autre codes: sortie immédiate
                if (codrep .ne. 3) goto 999
            endif
        enddo
    else
        call utmess('F', 'ELEMENTS2_39', sk=rela_comp)
    endif
!
999 continue
end subroutine
