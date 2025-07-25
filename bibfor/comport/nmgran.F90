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

subroutine nmgran(fami, kpg, ksp, typmod, imate, &
                  compor, instam, instap, tpmxm, tpmxp, &
                  depst, sigm, vim, option, sigp, &
                  vip, dsidep, materi)
!
    implicit none
#include "asterf_types.h"
#include "asterc/r8t0.h"
#include "asterfort/granvi.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcvalb.h"
#include "asterfort/rcvarc.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: imate, kpg, ksp
    character(len=8) :: typmod(*), materi
    character(len=16) :: compor(*), option, phenbid
    character(len=*) :: fami
    real(kind=8) :: instam, instap
    real(kind=8) :: tpmxm, tpmxp
    real(kind=8) :: depst(6)
    real(kind=8) :: sigm(6), vim(55), sigp(6), vip(55), dsidep(6, 6)
! ----------------------------------------------------------------------
!     REALISE LA LOI VISCOELASTIQUE DE GRANGER (FLUAGE PROPRE)
!     ELEMENTS ISOPARAMETRIQUES EN PETITES DEFORMATIONS
!
! IN  NDIM    : DIMENSION DE L'ESPACE
! IN  TYPMOD  : TYPE DE MODELISATION
! IN  IMATE   : ADRESSE DU MATERIAU CODE
! IN  COMPOR  : COMPORTEMENT : RELCOM ET DEFORM
! IN  CRIT    : CRITERES DE CONVERGENCE LOCAUX
! IN  INSTAM  : INSTANT DU CALCUL PRECEDENT
! IN  INSTAP  : INSTANT DU CALCUL
! IN  TPMXM   : TEMPERATURE MAX ATTEINTE AU COURS DE L'HISTORIQUE DE
!               CHARGEMENT A T (POUR LE COUPLAGE FLUAGE/FISSURATION)
! IN  TPMXP   : TEMPERATURE MAX ATTEINTE AU COURS DE L'HISTORIQUE DE
!               CHARGEMENT A T+DT (POUR LE COUPLAGE FLUAGE/FISSURATION)
! IN  DEPST   : INCREMENT DE DEFORMATION
!               SI C_PLAN DEPST(3) EST EN FAIT INCONNU (ICI:0)
!                 =>  ATTENTION LA PLACE DE DEPST(3) EST ALORS UTILISEE.
! IN  SIGM    : CONTRAINTES A L'INSTANT DU CALCUL PRECEDENT
! IN  VIM     : VARIABLES INTERNES A L'INSTANT DU CALCUL PRECEDENT
! IN  OPTION  : OPTION DEMANDEE : RIGI_MECA_TANG , FULL_MECA , RAPH_MECA
! OUT SIGP    : CONTRAINTES A L'INSTANT ACTUEL
! OUT VIP     : VARIABLES INTERNES A L'INSTANT ACTUEL
! OUT DSIDEP  : MATRICE CARREE (INUTILISE POUR RAPH_MECA)
!
!               ATTENTION LES TENSEURS ET MATRICES SONT RANGES DANS
!               L'ORDRE :  XX,YY,ZZ,SQRT(2)*XY,SQRT(2)*XZ,SQRT(2)*YZ
!
! ----------------------------------------------------------------------
!
    real(kind=8) :: valres(16), tm, tp, tref
    real(kind=8) :: e, nu, troisk, deuxmu
    real(kind=8) :: delta, dteqt, agem, agep, dage, tceq
    real(kind=8) :: temp, tabs, tkm, tkp, tmprim, tpprim, tkref
    real(kind=8) :: kron(6), hydrm, hydrp, sechm, sechp, sref
    real(kind=8) :: em, num, troikm, deumum
    real(kind=8) :: depsmo, depsdv(6), sigmmo, sigmdv(6), sigpmo, sigpdv(6)
    real(kind=8) :: smdv(6), spdv(6), smmo, spmo, dsdv(6), dsmo, deps(6), deps3
    real(kind=8) :: sigldv(6), siglmo, sigmp(6), sigmpo, sigmpd(6)
    integer(kind=8) :: ndimsi, ibid, iret, ibid2
    integer(kind=8) :: i, k, l, n, iret1, iret2, iret3
    integer(kind=8) :: icodre(16), ndt, ndi
    character(len=16) :: nomres(16)
    character(len=8) :: nompar, mod
    real(kind=8) :: valpam, valpap
    real(kind=8) :: bendom, bendop, kdessm, kdessp
    real(kind=8) :: j(8), taux(8), hygrm, hygrp, vieil
    real(kind=8) :: amdv(6, 9), apdv(6, 9), ammo(9), apmo(9), ap(6, 9), am(6, 9)
    real(kind=8) :: ther, coefa(9), coefc(9), coeff(9), coefb, coefd
    real(kind=8) :: coefg, coefh, coefi, coefj, coefv, coefk(9), epsthp, epsthm
    aster_logical :: cplan
    data kron/1.d0, 1.d0, 1.d0, 0.d0, 0.d0, 0.d0/

    common/tdim/ndt, ndi
! DEB -----------------------------------------------------------------
!
    call rcvarc('F', 'TEMP', '-', fami, kpg, ksp, tm, iret2)
    call rcvarc('F', 'TEMP', '+', fami, kpg, ksp, tp, iret3)
    call rcvarc('F', 'TEMP', 'REF', fami, kpg, ksp, tref, iret1)
    if ((iret1+iret2+iret3) .ge. 1) then
        call utmess('F', 'COMPOR5_40', sk=compor(1))
    end if
    nompar = 'TEMP'
!
!     -- 1 INITIALISATIONS :
!     ----------------------
!
!
    cplan = typmod(1) .eq. 'C_PLAN'
!
!      NDIMSI = 2*NDIM
    mod = typmod(1)
    call granvi(mod, ndimsi, ibid, ibid2)
    call granvi(mod, ndt, ndi)
    dsidep(:, :) = 0.d0
    if (.not. (compor(1) (1:13) .eq. 'BETON_GRANGER')) then
        call utmess('F', 'ALGORITH4_50', sk=compor(1))
    end if
    delta = instap-instam
    temp = (tp+tm)/2
    tabs = r8t0()
    temp = temp+tabs
    tkm = tm+tabs
    tkp = tp+tabs
    tkref = tref+tabs
    l = 0
    do i = 1, 9
        do k = 1, ndimsi
            l = l+1
            am(k, i) = vim(l)
        end do
        ammo(i) = 0.d0
        do n = 1, 3
            ammo(i) = ammo(i)+am(n, i)
        end do
        ammo(i) = ammo(i)/3.d0
        do k = 1, ndimsi
            amdv(k, i) = am(k, i)-ammo(i)*kron(k)
        end do
    end do
    agem = vim(l+1)
!
!     -- 2 RECUPERATION DES CARACTERISTIQUES
!     ---------------------------------------
    nomres(1) = 'E'
    nomres(2) = 'NU'
    nomres(3) = 'ALPHA'
!
    nompar = 'TEMP'
    valpam = tpmxm
    valpap = tpmxp
!
    epsthm = 0.0; epsthp = 0.0
    call rcvalb(fami, kpg, ksp, '-', imate, materi, 'ELAS', &
                1, nompar, [valpam], 2, nomres, valres, icodre, 2)
!
    call rcvalb(fami, kpg, ksp, '-', imate, materi, 'ELAS', &
                1, nompar, [valpam], 1, nomres(3), valres(3), icodre(3), 2)
    if ((.not. isnan(tp)) .and. (.not. isnan(tm))) then
        if ((icodre(3) .ne. 0) .or. (isnan(tref))) then
            call utmess('F', 'COMPOR5_42')
        else
            epsthm = valres(3)*(tm-tref)
        end if
    else
        valres(3) = 0.d0
        epsthm = 0.d0
    end if
    em = valres(1)
    num = valres(2)
    deumum = em/(1.d0+num)
    troikm = em/(1.d0-2.d0*num)
!
!
    call rcvalb(fami, kpg, ksp, '+', imate, materi, 'ELAS', &
                1, nompar, [valpap], 2, nomres, valres, icodre, 2)
!
    call rcvalb(fami, kpg, ksp, '+', imate, materi, 'ELAS', &
                1, nompar, [valpap], 1, nomres(3), valres(3), icodre(3), 0)
!
    if (icodre(3) .ne. 0) then
        valres(3) = 0.d0
        epsthp = 0.d0
    else if (((isnan(tp)) .or. (isnan(tm)) .or. (isnan(tref))) .and. (icodre(3) .eq. 0)) then
        call utmess('F', 'COMPOR5_42')
    else
        epsthp = valres(3)*(tp-tref)
    end if
    e = valres(1)
    nu = valres(2)
    deuxmu = e/(1.d0+nu)
    troisk = e/(1.d0-2.d0*nu)
!
! ------- CARAC. RETRAIT ENDOGENE ET RETRAIT DE DESSICCATION
!
    nomres(1) = 'B_ENDOGE'
    nomres(2) = 'K_DESSIC'
    call rcvalb(fami, kpg, ksp, '-', imate, materi, 'ELAS', &
                1, nompar, [valpam], 1, nomres, valres, icodre, 0)
!
    if (icodre(1) .ne. 0) valres(1) = 0.d0
    bendom = valres(1)
!
    call rcvalb(fami, kpg, ksp, '+', imate, materi, 'ELAS', &
                1, nompar, [valpap], 1, nomres, valres, icodre, 0)
!
    if (icodre(1) .ne. 0) valres(1) = 0.d0
    bendop = valres(1)
!
    call rcvalb(fami, kpg, ksp, '-', imate, materi, 'ELAS', &
                1, nompar, [valpam], 1, nomres(2), valres(2), icodre, 0)
!
    if (icodre(2) .ne. 0) valres(2) = 0.d0
    kdessm = valres(2)
!
    call rcvalb(fami, kpg, ksp, '+', imate, materi, 'ELAS', &
                1, nompar, [valpap], 1, nomres(2), valres(2), icodre, 0)
!
    if (icodre(2) .ne. 0) valres(2) = 0.d0
    kdessp = valres(2)
!
!  ------- CARACTERISTIQUES FONCTION DE FLUAGE
!
    nomres(1) = 'J1'
    nomres(2) = 'J2'
    nomres(3) = 'J3'
    nomres(4) = 'J4'
    nomres(5) = 'J5'
    nomres(6) = 'J6'
    nomres(7) = 'J7'
    nomres(8) = 'J8'
    nomres(9) = 'TAUX_1'
    nomres(10) = 'TAUX_2'
    nomres(11) = 'TAUX_3'
    nomres(12) = 'TAUX_4'
    nomres(13) = 'TAUX_5'
    nomres(14) = 'TAUX_6'
    nomres(15) = 'TAUX_7'
    nomres(16) = 'TAUX_8'
    coefj = 0.d0
    do i = 1, 8
        call rcvalb(fami, kpg, ksp, '+', imate, materi, 'BETON_GRANGER', &
                    0, ' ', [0.d0], 1, nomres(i), valres(i), icodre(i), 0)
        call rcvalb(fami, kpg, ksp, '+', imate, materi, 'BETON_GRANGER', &
                    0, ' ', [0.d0], 1, nomres(i+8), valres(i+8), icodre(i+8), 0)
        if ((icodre(i) .ne. 0) .and. (icodre(i+8) .ne. 0)) then
            valres(i) = 0.d0
            valres(i+8) = 1.d0
        else if (((icodre(i) .eq. 0) .and. (icodre(i+8) .ne. 0)) .or. &
                 ((icodre(i) .ne. 0) .and. (icodre(i+8) .eq. 0))) then
            call utmess('F', 'ALGORITH8_2')
        end if
        j(i) = valres(i)
        taux(i) = valres(i+8)
        coefj = coefj+j(i)
    end do
!
!  ------- CARACTERISTIQUES EFFET DE LA TEMPERATURE
! Eliminee la dependance a la temperature le 19 mai 2016
    dteqt = delta
!
    tmprim = 1.d0
    tpprim = 1.d0
!
!  ------- CARACTERISTIQUES EFFET DU VIEILLISSEMENT
!
    if (compor(1) (1:15) .eq. 'BETON_GRANGER_V') then
!       FONCTION MULTIPLICATIVE - VIEILLISSEMENT K
!           AGE EQUIVALENT DU BETON : AGE
!               coefv=(-qsrv)*(1/temp-1/tkref)
!               coefv=exp(coefv)
        coefv = 1.d0
        dage = coefv*delta
        agep = agem+dage
        tceq = (agem+agep)/2
        nomres(1) = 'FONC_V'
        call rcvalb(fami, kpg, ksp, '+', imate, materi, 'V_BETON_GRANGER', &
                    1, 'INST', [tceq], 1, nomres, valres, icodre, 0)
!
        vieil = valres(1)
    else
        vieil = 1.d0
        dage = delta
        agep = agem+dage
    end if
!
!  ------- CARACTERISTIQUES HYGROMETRIE H
!
    call rccoma(imate, 'BETON_DESORP', 0, phenbid, icodre(1))
    if (icodre(1) .ne. 0) then
        call utmess('F', 'ALGORITH7_97')
    end if

    nomres(1) = 'FONC_DESORP'
    call rcvalb(fami, kpg, ksp, '-', imate, materi, 'BETON_DESORP', &
                1, nompar, [valpam], 1, nomres, valres, icodre, 0)
!
    if (icodre(1) .ne. 0) then
        call utmess('F', 'ALGORITH7_98')
    end if
    hygrm = valres(1)
    call rcvalb(fami, kpg, ksp, '+', imate, materi, 'BETON_DESORP', &
                1, nompar, [valpap], 1, nomres, valres, icodre, 0)
!
    if (icodre(1) .ne. 0) then
        call utmess('F', 'ALGORITH7_98')
    end if
    hygrp = valres(1)
!
!
!--3  CALCUL DE SIGMP, SIGEL
!   -------------------------
    sigmmo = 0.d0
    do k = 1, 3
        sigmmo = sigmmo+sigm(k)
    end do
    sigmmo = sigmmo/3.d0
    do k = 1, ndimsi
        sigmdv(k) = sigm(k)-sigmmo*kron(k)
        sigmp(k) = deuxmu/deumum*sigmdv(k)+troisk/troikm*sigmmo*kron( &
                   k)
    end do
    sigmpo = 0.d0
    do k = 1, 3
        sigmpo = sigmpo+sigmp(k)
    end do
    sigmpo = sigmpo/3.d0
    do k = 1, ndimsi
        sigmpd(k) = sigmp(k)-sigmpo*kron(k)
    end do
!
! -------  QUELQUES COEFFICIENTS-------------------------
!
    coefb = 0.d0
    do i = 1, 8
        coefb = coefb+j(i)*(1-(taux(i)/delta)*(1-exp(-dteqt/taux(i))))
    end do
    coefd = (hygrp*tpprim*vieil)*coefb
    do k = 1, ndimsi
        coefa(k) = 0.d0
        do i = 1, 8
            coefa(k) = coefa(k)+amdv(k, i)*(1-exp(-dteqt/taux(i)))
        end do
        coefc(k) = (sigmdv(k)*hygrm*tmprim)*vieil*coefb
        coeff(k) = coefa(k)-coefc(k)
    end do
    coefa(9) = 0.d0
    do i = 1, 8
        coefa(9) = coefa(9)+ammo(i)*(1-exp(-dteqt/taux(i)))
    end do
    coefc(9) = (sigmmo*hygrm*tmprim)*vieil*coefb
    coeff(9) = coefa(9)-coefc(9)
    do k = 1, 3
        coefk(k) = coeff(k)+coeff(9)
    end do
!
! ------- CALCUL DE DEPSMO ET DEPSDV
!
!
    call rcvarc(' ', 'HYDR', '+', fami, kpg, ksp, hydrp, iret)
    if (iret .ne. 0) hydrp = 0.d0
    call rcvarc(' ', 'HYDR', '-', fami, kpg, ksp, hydrm, iret)
    if (iret .ne. 0) hydrm = 0.d0
    call rcvarc(' ', 'SECH', '+', fami, kpg, ksp, sechp, iret)
    if (iret .ne. 0) sechp = 0.d0
    call rcvarc(' ', 'SECH', '-', fami, kpg, ksp, sechm, iret)
    if (iret .ne. 0) sechm = 0.d0
    call rcvarc(' ', 'SECH', 'REF', fami, kpg, ksp, sref, iret)
    if (iret .ne. 0) sref = 0.d0
    ther = epsthp-epsthm-bendop*hydrp+bendom*hydrm-kdessp*(sref-sechp)+kdessm*(sref-sechm)
    if (cplan) then
        deps3 = -nu/(1.d0-nu)*(depst(1)+depst(2))+(1.d0+nu)/(1.d0-nu)*ther
        deps3 = deps3+(3.d0*e/(deuxmu*2.d0+troisk))*coefk(3)
        deps3 = deps3-(3.d0/(deuxmu*2.d0+troisk))*sigmp(3)
    end if
    depsmo = 0.d0
    do k = 1, 3
        deps(k) = depst(k)-ther
        deps(k+3) = depst(k+3)
        depsmo = depsmo+deps(k)
    end do
    if (cplan) then
        deps(3) = deps3-ther
        depsmo = deps(1)+deps(2)+deps(3)
    end if
    depsmo = depsmo/3.d0
    do k = 1, ndimsi
        depsdv(k) = deps(k)-depsmo*kron(k)
    end do
!
    do k = 1, ndimsi
        sigldv(k) = sigmpd(k)+deuxmu*depsdv(k)
    end do
    siglmo = sigmpo+troisk*depsmo
!
!--4-CALCUL DE SIGP
!------------------------
!
    if ((option(1:9) .eq. 'RAPH_MECA') .or. (option(1:9) .eq. 'FULL_MECA')) then
        do k = 1, ndimsi
            sigpdv(k) = sigldv(k)-e*coeff(k)
            sigpdv(k) = sigpdv(k)/(1+e*coefd)
        end do
        sigpmo = siglmo-e*coeff(9)
        sigpmo = sigpmo/(1+e*coefd)
        do k = 1, ndimsi
            sigp(k) = sigpdv(k)+sigpmo*kron(k)
        end do
!
!-- 6 CALCUL DE VIP
!   -------------------
!--------CALCUL DE DS
        do k = 1, ndimsi
            smdv(k) = sigmdv(k)*hygrm*tmprim
            spdv(k) = sigpdv(k)*hygrp*tpprim
            dsdv(k) = spdv(k)-smdv(k)
        end do
        smmo = sigmmo*hygrm*tmprim
        spmo = sigpmo*hygrp*tpprim
        dsmo = spmo-smmo
!--------CALCUL DE APDV, APMO
        do i = 1, 8
            do k = 1, ndimsi
                apdv(k, i) = amdv(k, i)*exp(-dteqt/taux(i))
                apdv(k, i) = apdv(k, i)+dsdv(k)*vieil*j(i)*(taux(i)/delta)*(1-exp(-dteqt/taux(i)))
            end do
            apmo(i) = ammo(i)*exp(-dteqt/taux(i))
            apmo(i) = apmo(i)+dsmo*vieil*j(i)*(taux(i)/delta)*(1-exp(-dteqt/taux(i)))
        end do
!
        do k = 1, ndimsi
            apdv(k, 9) = amdv(k, 9)
            do i = 1, 8
                apdv(k, 9) = apdv(k, 9)+j(i)*vieil*dsdv(k)
            end do
        end do
        apmo(9) = ammo(9)
        do i = 1, 8
            apmo(9) = apmo(9)+j(i)*vieil*dsmo
        end do
!--------CALCUL DES AP ET DES VIP
        l = 0
        do i = 1, 9
            do k = 1, ndimsi
                ap(k, i) = apdv(k, i)+apmo(i)*kron(k)
                l = l+1
                vip(l) = ap(k, i)
            end do
        end do
        vip(l+1) = agep
    end if
!
!-- 7 CALCUL DE DSIDEP POUR LA MATRICE TANGENTE
!   -------------------------------------------
    if ((option(1:14) .eq. 'RIGI_MECA_TANG') .or. (option(1:9) .eq. 'FULL_MECA')) then
        if (option(1:9) .eq. 'FULL_MECA') then
            coefg = coefd
            coefh = (1+e*coefg)
            coefi = (1+e*coefg)
        else
!        (METRICE ELASTIQUE)
            coefh = 1.d0
            coefi = 1.d0
        end if
        do k = 1, 3
            do l = 1, 3
                dsidep(k, l) = dsidep(k, l)+(troisk/(3.d0*coefi))
                dsidep(k, l) = dsidep(k, l)-deuxmu/(3.d0*coefh)
            end do
        end do
        do k = 1, ndimsi
            dsidep(k, k) = dsidep(k, k)+deuxmu/coefh
        end do
!
!----------- CORRECTION POUR LES CONTRAINTES PLANES :
        if (cplan) then
            do k = 1, ndimsi
                if (k .ne. 3) then
                    do l = 1, ndimsi
                        if (l .ne. 3) then
                            dsidep(k, l) = dsidep(k, l)-1.d0/dsidep(3, 3)*dsidep(k, 3)*dsidep(3, l)
                        end if
                    end do
                end if
            end do
        end if
    end if
end subroutine
