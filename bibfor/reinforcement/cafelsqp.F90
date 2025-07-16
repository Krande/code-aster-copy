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

subroutine cafelsqp(cequi, effm, effn, ht, bw, &
                    enrobi, enrobs, wmaxi, wmaxs, &
                    ferrcomp, precs, ferrsyme, slsyme, uc, um, &
                    kt, facier, fbeton, eys, sigelsqp, phiinf, phisup, &
                    dnsinf, dnssup, sigmsi, sigmss, sigmci, sigmcs, &
                    alpha, pivot, etat, &
                    wfini, wfins, kvarf, ierr)
!_______________________________________________________________________________________________
!
!   CAFELSQP
!
!   CALCUL DES ACIERS EN FLEXION COMPOSEE A L'ELS QP
!   CRITERE = LIMITATION DES OUVERTURES DES FISSURES
!
!   I CEQUI        COEFFICIENT D'EQUIVALENCE ACIER/BETON
!   I EFFM         MOMENT DE FLEXION
!   I EFFN         EFFORT NORMAL
!   I HT           EPAISSEUR DE LA SECTION
!   I BW           LARGEUR DE LA SECTION
!   I ENROBI       ENROBAGE DES ARMATURES INFERIEURES
!   I ENROBS       ENROBAGE DES ARMATURES SUPERIEURES
!   I WMAXI        OUVERTURE MAX DE FISSURES ADMISSIBLE EN FIBRE INFERIEURE
!   I WMAXS        OUVERTURE MAX DE FISSURES ADMISSIBLE EN FIBRE SUPERIEURE
!   I FERRCOMP     PRISE EN COMPTE DU FERRAILLAGE DE COMPRESSION
!                         FERRCOMP = 0 (NON)
!                         FERRCOMP = 1 (OUI)
!   I PRECS        PRECISION ITERATION
!   I FERRSYME     FERRAILLAGE SYMETRIQUE?
!                         FERRSYME = 0 (NON)
!                         FERRSYME = 1 (OUI)
!   I SLSYME       SECTION SEUIL DE TOLERANCE POUR UN FERRAILLAGE SYMETRIQUE
!   I UC           UNITE DES CONTRAINTES :
!                         UC = 0 CONTRAINTES EN Pa
!                         UC = 1 CONTRAINTES EN MPa
!   I UM           UNITE DES DIMENSIONS :
!                         UM = 0 DIMENSIONS EN m
!                         UM = 1 DIMENSIONS EN mm
!   I KT           COEFFICIENT DE DUREE DE CHARGEMENT
!   I FACIER       LIMITE D'ELASTICITÉ DES ACIERS (CONTRAINTE)
!   I FBETON       RESISTANCE EN COMPRESSION DU BÉTON (CONTRAINTE)
!   I EYS          MODULE D'YOUNG DE L'ACIER
!   I SIGELSQP     CONTRAINTE ADMISSIBLE DANS LE BETON À L'ELS QP
!   I PHIINF       DIAMETRE ESTIMATIF DES ARMATURES EN FIBRE INFERIEURE
!   I PHISUP       DIAMETRE ESTIMATIF DES ARMATURES EN FIBRE SUPERIEURE

!   O DNSINF       DENSITE DE L'ACIER INFERIEUR
!   O DNSSUP       DENSITE DE L'ACIER SUPERIEUR
!   O SIGMSI       CONTRAINTE AU NIVEAU DE L'ACIER INFERIEUR
!   O SIGMSS       CONTRAINTE AU NIVEAU DE L'ACIER SUPERIEUR
!   O SIGMCI       CONTRAINTE AU NIVEAU DE LA FIBRE SUPERIEURE DE BETON
!   O SIGMCS       CONTRAINTE AU NIVEAU DE LA FIBRE INFERIEURE DE BETON
!   O ALPHA        COEFFICIENT DE PROFONDEUR DE L'AN
!   O PIVOT        PIVOT DE FONCTIONNEMENT DE LA SECTION
!   O ETAT         ETAT DE FONCTIONNEMENT DE LA SECTION
!   O WFINI        OUVERTURE DES FISSURES FINALE EN FIBRE INFERIEURE
!   O WFINS        OUVERTURE DES FISSURES FINALE EN FIBRE SUPERIEURE
!   O KVARF        TAUX DE CONTRAINTE DE TRACTION LIMITE DANS L'ACIER
!   O IERR         CODE RETOUR (0 = OK)
!_______________________________________________________________________________________________
!

    implicit none
#include "asterfort/cafels.h"

    real(kind=8) :: cequi
    real(kind=8) :: effm
    real(kind=8) :: effn
    real(kind=8) :: ht
    real(kind=8) :: bw
    real(kind=8) :: enrobi
    real(kind=8) :: enrobs
    real(kind=8) :: wmaxi
    real(kind=8) :: wmaxs
    integer(kind=8) :: ferrcomp
    integer(kind=8) :: precs
    integer(kind=8) :: ferrsyme
    real(kind=8) :: slsyme
    integer(kind=8) :: uc
    integer(kind=8) :: um
    real(kind=8) :: kt
    real(kind=8) :: facier
    real(kind=8) :: fbeton
    real(kind=8) :: eys
    real(kind=8) :: sigelsqp
    real(kind=8) :: phiinf
    real(kind=8) :: phisup
    real(kind=8) :: dnsinf
    real(kind=8) :: dnssup
    real(kind=8) :: sigmsi
    real(kind=8) :: sigmss
    real(kind=8) :: sigmci
    real(kind=8) :: sigmcs
    real(kind=8) :: alpha
    integer(kind=8) :: pivot
    integer(kind=8) :: etat
    real(kind=8) :: wfini
    real(kind=8) :: wfins
    real(kind=8) :: kvarf
    integer(kind=8) :: ierr

!   COEFFICIENTS LIE A L'UNITE CHOISIE
    real(kind=8) :: unite_pa, unite_m
!
!   BRAS DE LEVIER DE CALCUL
    real(kind=8) :: d, d0

!   RESISTANCE MOYENNE DU BETON EN TRACTION
    real(kind=8) :: fctm

!   VARIABLES DE CALCUL DES OUVERTURES DES FISSURES A L'ELS QP
    real(kind=8) :: f1, f2, f3, f4, f5, ssmax
    real(kind=8) :: scmaxi, scmaxs
    real(kind=8) :: ScTracMAX, ScTracMIN, xAN
    real(kind=8) :: hceff, rhoeff, SrmaxSUP, SrmaxINF, DESUP, DEINF
    logical :: COND_SUP, COND_INF
    real(kind=8) :: kVAR_A, kVAR_B, kVAR_NEW
    integer(kind=8) :: COUNT_B
    logical :: COND_A, COND_B, CONSTAT, COND_FISS, COND_NEW
    integer(kind=8) :: Niter

!   INITIALISATION DU CODE RETOUR
    ierr = 0
    etat = 0
    pivot = 0
    alpha = -1000
    dnssup = -1
    dnsinf = -1
    sigmss = -1
    sigmsi = -1
    sigmcs = -1
    sigmci = -1
    wfini = -1
    wfins = -1
    kvarf = -1

    if (uc .eq. 0) then
        unite_pa = 1.e-6
    elseif (uc .eq. 1) then
        unite_pa = 1.
    end if
    if (um .eq. 0) then
        unite_m = 1.e3
    elseif (um .eq. 1) then
        unite_m = 1.
    end if

    if (effm .le. 0) then
        d = ht-enrobi
        d0 = ht-enrobs
    else
        d = ht-enrobs
        d0 = ht-enrobi
    end if

!   CALCULS INTERMEDIAIRES

    f1 = 0.8
    f4 = 0.425
    f5 = 1.0
    if ((fbeton*unite_pa) .le. 50) then
        fctm = 0.30*((fbeton*unite_pa)**(2.0/3.0))
    else
        fctm = 2.12*LOG(1.0+((fbeton*unite_pa)+8.0)/10.0)
    end if

    fctm = fctm/unite_pa

    kVAR_A = 1.0
    kVAR_B = 0.5
    COND_A = .false.
    COND_B = .false.
    COND_FISS = .false.
    COUNT_B = 0
    CONSTAT = .false.
    scmaxi = sigelsqp
    scmaxs = sigelsqp

!   1 - CALCUL DES DENSITES DE FERRAILLAGE A L'ELS CARA (kvarA = 1.0)
!   ----------------------------------------------------------------------

    ssmax = kVAR_A*facier
    call cafels(cequi, effm, effn, ht, bw, &
                enrobi, enrobs, scmaxs, scmaxi, ssmax, &
                ferrcomp, precs, ferrsyme, slsyme, uc, um, &
                dnsinf, dnssup, sigmsi, sigmss, &
                sigmci, sigmcs, &
                alpha, pivot, etat, ierr)

    if (ierr .ne. 0) then
        goto 998
    end if

!   2 - VERIFICATION DE L'OUVERTURE DES FISSURES À L'ELS QP (kvarA = 1.0)
!   ----------------------------------------------------------------------

    if ((alpha .ge. 0) .AND. (alpha .le. 1)) then
        if (effm .gt. 0) then
            xAN = alpha*(ht-enrobi)
        elseif (effm .lt. 0) then
            xAN = alpha*(ht-enrobs)
        end if
    else
        xAN = -1
    end if

    ScTracMAX = min(sigmci, sigmcs)
    ScTracMIN = max(sigmci, sigmcs)
    if (ScTracMIN .gt. 0) then
        ScTracMIN = 0
    end if
    if (ScTracMAX .ge. 0) then
        f2 = 0.5
    else
        f2 = (ScTracMIN+ScTracMAX)/(2*ScTracMAX)
    end if

    if ((sigmss .ge. 0) .OR. (dnssup .eq. 0) .OR. (etat .eq. 1)) then
        wfins = 0
        COND_SUP = .true.
    else
        if ((enrobs*unite_m) .lt. 25) then
            f3 = 3.4
        else
            f3 = 3.4*((25.0/(enrobs*unite_m))**(2.0/3.0))
        end if
        if (xAN .ne. -1) then
            hceff = min(2.5*enrobs, 0.5*ht, (ht-xAN)/3.0)
        else
            hceff = min(2.5*enrobs, 0.5*ht)
        end if
        rhoeff = dnssup/hceff
        SrmaxSUP = f3*enrobs+f1*f2*f4*(phisup/rhoeff)
        if (xAN .ne. -1) then
            SrmaxSUP = min(SrmaxSUP, 1.3*(ht-xAN))
        else
            SrmaxSUP = min(SrmaxSUP, 1.3*ht)
        end if
        DESUP = (-sigmss-kt*(fctm/rhoeff)*(1.0+cequi*rhoeff))/eys
        DESUP = max(DESUP, -0.6*sigmss/eys)
        wfins = SrmaxSUP*DESUP
        if (wfins .le. wmaxs) then
            COND_SUP = .true.
        else
            COND_SUP = .false.
        end if
    end if

    if ((sigmsi .ge. 0) .OR. (dnsinf .eq. 0) .OR. (etat .eq. 1)) then
        wfini = 0
        COND_INF = .true.
    else
        if ((enrobi*unite_m) .lt. 25) then
            f3 = 3.4
        else
            f3 = 3.4*((25.0/(enrobi*unite_m))**(2.0/3.0))
        end if
        if (xAN .ne. -1) then
            hceff = min(2.5*enrobi, 0.5*ht, (ht-xAN)/3.0)
        else
            hceff = min(2.5*enrobi, 0.5*ht)
        end if
        rhoeff = dnsinf/hceff
        SrmaxINF = f3*enrobi+f1*f2*f4*(phiinf/rhoeff)
        if (xAN .ne. -1) then
            SrmaxINF = min(SrmaxINF, 1.3*(ht-xAN))
        else
            SrmaxINF = min(SrmaxINF, 1.3*ht)
        end if
        DEINF = (-sigmsi-kt*(fctm/rhoeff)*(1.0+cequi*rhoeff))/eys
        DEINF = max(DEINF, -0.6*sigmsi/eys)
        wfini = SrmaxINF*DEINF
        if (wfini .le. wmaxi) then
            COND_INF = .true.
        else
            COND_INF = .false.
        end if
    end if

    if ((COND_SUP .eqv. (.false.)) .or. (COND_INF .eqv. (.false.))) then
        COND_A = .false.
    else
        COND_A = .true.
        COND_FISS = .true.
    end if

    Niter = 20

    do while ((COND_B .eqv. (.false.)) .and. (COND_FISS .eqv. (.false.)) .and. (COUNT_B .le. Niter))

!   3 - CALCUL DES DENSITES DE FERRAILLAGE A L'ELS CARA (kvarB)
!   ----------------------------------------------------------------------

        ssmax = kVAR_B*facier
        call cafels(cequi, effm, effn, ht, bw, &
                    enrobi, enrobs, scmaxs, scmaxi, ssmax, &
                    ferrcomp, precs, ferrsyme, slsyme, uc, um, &
                    dnsinf, dnssup, sigmsi, sigmss, &
                    sigmci, sigmcs, &
                    alpha, pivot, etat, ierr)

        if (ierr .ne. 0) then
            goto 998
        end if

!   4 - VERIFICATION DE L'OUVERTURE DES FISSURES À L'ELS QP (kvarB)
!   ----------------------------------------------------------------------

        if ((alpha .ge. 0) .AND. (alpha .le. 1)) then
            if (effm .gt. 0) then
                xAN = alpha*(ht-enrobi)
            elseif (effm .lt. 0) then
                xAN = alpha*(ht-enrobs)
            end if
        else
            xAN = -1
        end if

        ScTracMAX = min(sigmci, sigmcs)
        ScTracMIN = max(sigmci, sigmcs)
        if (ScTracMIN .gt. 0) then
            ScTracMIN = 0
        end if
        if (ScTracMAX .ge. 0) then
            f2 = 0.5
        else
            f2 = (ScTracMIN+ScTracMAX)/(2*ScTracMAX)
        end if
        if ((sigmss .ge. 0) .OR. (dnssup .eq. 0) .OR. (etat .eq. 1)) then
            wfins = 0
            COND_SUP = .true.
        else
            if ((enrobs*unite_m) .lt. 25) then
                f3 = 3.4
            else
                f3 = 3.4*((25.0/(enrobs*unite_m))**(2.0/3.0))
            end if
            if (xAN .ne. -1) then
                hceff = min(2.5*enrobs, 0.5*ht, (ht-xAN)/3.0)
            else
                hceff = min(2.5*enrobs, 0.5*ht)
            end if
            rhoeff = dnssup/hceff
            SrmaxSUP = f3*enrobs+f1*f2*f4*(phisup/rhoeff)
            if (xAN .ne. -1) then
                SrmaxSUP = min(SrmaxSUP, 1.3*(ht-xAN))
            else
                SrmaxSUP = min(SrmaxSUP, 1.3*ht)
            end if
            DESUP = (-sigmss-kt*(fctm/rhoeff)*(1.0+cequi*rhoeff))/eys
            DESUP = max(DESUP, -0.6*sigmss/eys)
            wfins = SrmaxSUP*DESUP
            if (wfins .le. wmaxs) then
                COND_SUP = .true.
            else
                COND_SUP = .false.
            end if
        end if

        if ((sigmsi .ge. 0) .OR. (dnsinf .eq. 0) .OR. (etat .eq. 1)) then
            wfini = 0
            COND_INF = .true.
        else
            if ((enrobi*unite_m) .lt. 25) then
                f3 = 3.4
            else
                f3 = 3.4*((25.0/(enrobi*unite_m))**(2.0/3.0))
            end if
            if (xAN .ne. -1) then
                hceff = min(2.5*enrobi, 0.5*ht, (ht-xAN)/3.0)
            else
                hceff = min(2.5*enrobi, 0.5*ht)
            end if
            rhoeff = dnsinf/hceff
            SrmaxINF = f3*enrobi+f1*f2*f4*(phiinf/rhoeff)
            if (xAN .ne. -1) then
                SrmaxINF = min(SrmaxINF, 1.3*(ht-xAN))
            else
                SrmaxINF = min(SrmaxINF, 1.3*ht)
            end if
            DEINF = (-sigmsi-kt*(fctm/rhoeff)*(1.0+cequi*rhoeff))/eys
            DEINF = max(DEINF, -0.6*sigmsi/eys)
            wfini = SrmaxINF*DEINF
            if (wfini .le. wmaxi) then
                COND_INF = .true.
            else
                COND_INF = .false.
            end if
        end if

        if ((COND_SUP .eqv. (.false.)) .or. (COND_INF .eqv. (.false.))) then
            COND_B = .false.
        else
            COND_B = .true.
        end if

        if (COND_B .eqv. (.false.)) then
            kVAR_B = kVAR_B/2.
            COUNT_B = COUNT_B+1
        end if

!   ----------------------------------------------------------------------

    end do

!   5 - ITERATION ET DETERMINATION DE LA SOLUTION DIMENSIONNANTE
!   ----------------------------------------------------------------------

    if ((COND_FISS .eqv. (.false.)) .and. (COND_B .eqv. (.false.))) then
        CONSTAT = .true.
        ierr = 3
        goto 998
    elseif (COND_FISS .eqv. (.true.)) then
        CONSTAT = .true.
        kVAR_NEW = kVAR_A
        COND_NEW = COND_A
    end if

    do while (CONSTAT .eqv. (.false.))

        kVAR_NEW = 0.5*(kVAR_A+kVAR_B)

        ssmax = kVAR_NEW*facier
        call cafels(cequi, effm, effn, ht, bw, &
                    enrobi, enrobs, scmaxs, scmaxi, ssmax, &
                    ferrcomp, precs, ferrsyme, slsyme, uc, um, &
                    dnsinf, dnssup, sigmsi, sigmss, &
                    sigmci, sigmcs, &
                    alpha, pivot, etat, ierr)

        if (ierr .ne. 0) then
            goto 998
        end if

        if ((alpha .ge. 0) .AND. (alpha .le. 1)) then
        if (effm .gt. 0) then
            xAN = alpha*(ht-enrobi)
        elseif (effm .lt. 0) then
            xAN = alpha*(ht-enrobs)
        end if
        else
        xAN = -1
        end if
        ScTracMAX = min(sigmci, sigmcs)
        ScTracMIN = max(sigmci, sigmcs)
        if (ScTracMIN .gt. 0) then
            ScTracMIN = 0
        end if
        if (ScTracMAX .ge. 0) then
            f2 = 0.5
        else
            f2 = (ScTracMIN+ScTracMAX)/(2*ScTracMAX)
        end if
        if ((sigmss .ge. 0) .OR. (dnssup .eq. 0) .OR. (etat .eq. 1)) then
            wfins = 0
            COND_SUP = .true.
        else
            if ((enrobs*unite_m) .lt. 25) then
                f3 = 3.4
            else
                f3 = 3.4*((25.0/(enrobs*unite_m))**(2.0/3.0))
            end if
            if (xAN .ne. -1) then
                hceff = min(2.5*enrobs, 0.5*ht, (ht-xAN)/3.0)
            else
                hceff = min(2.5*enrobs, 0.5*ht)
            end if
            rhoeff = dnssup/hceff
            SrmaxSUP = f3*enrobs+f1*f2*f4*(phisup/rhoeff)
            if (xAN .ne. -1) then
                SrmaxSUP = min(SrmaxSUP, 1.3*(ht-xAN))
            else
                SrmaxSUP = min(SrmaxSUP, 1.3*ht)
            end if
            DESUP = (-sigmss-kt*(fctm/rhoeff)*(1.0+cequi*rhoeff))/eys
            DESUP = max(DESUP, -0.6*sigmss/eys)
            wfins = SrmaxSUP*DESUP
            if (wfins .le. wmaxs) then
                COND_SUP = .true.
            else
                COND_SUP = .false.
            end if
        end if
        if ((sigmsi .ge. 0) .OR. (dnsinf .eq. 0) .OR. (etat .eq. 1)) then
            wfini = 0
            COND_INF = .true.
        else
            if ((enrobi*unite_m) .lt. 25) then
                f3 = 3.4
            else
                f3 = 3.4*((25.0/(enrobi*unite_m))**(2.0/3.0))
            end if
            if (xAN .ne. -1) then
                hceff = min(2.5*enrobi, 0.5*ht, (ht-xAN)/3.0)
            else
                hceff = min(2.5*enrobi, 0.5*ht)
            end if
            rhoeff = dnsinf/hceff
            SrmaxINF = f3*enrobi+f1*f2*f4*(phiinf/rhoeff)
            if (xAN .ne. -1) then
                SrmaxINF = min(SrmaxINF, 1.3*(ht-xAN))
            else
                SrmaxINF = min(SrmaxINF, 1.3*ht)
            end if
            DEINF = (-sigmsi-kt*(fctm/rhoeff)*(1.0+cequi*rhoeff))/eys
            DEINF = max(DEINF, -0.6*sigmsi/eys)
            wfini = SrmaxINF*DEINF
            if (wfini .le. wmaxi) then
                COND_INF = .true.
            else
                COND_INF = .false.
            end if
        end if

        if ((COND_SUP .eqv. (.false.)) .or. (COND_INF .eqv. (.false.))) then
            COND_NEW = .false.
        else
            COND_NEW = .true.
        end if

        if (COND_NEW .eqv. COND_A) then
            kVAR_A = kVAR_NEW
        else
            kVAR_B = kVAR_NEW
        end if

        if ((abs(kVAR_A-kVAR_B)) .lt. (0.00001)) then
            CONSTAT = .true.
        end if

    end do

!   6 - CHOIX DE LA SOLUTION FINALE
!   ----------------------------------------------------------------------

998 continue

    if (ierr .eq. 0) then
        if (COND_FISS .eqv. (.false.)) then
            if ((COND_A .eqv. (.true.)) .and. (COND_B .eqv. (.false.))) then
                kVAR_NEW = kVAR_A
            elseif ((COND_A .eqv. (.false.)) .and. (COND_B .eqv. (.true.))) then
                kVAR_NEW = kVAR_B
            else
                ierr = 1
            end if
        end if
        if (ierr .eq. 0) then
            ssmax = kVAR_NEW*facier
            call cafels(cequi, effm, effn, ht, bw, &
                        enrobi, enrobs, scmaxs, scmaxi, ssmax, &
                        ferrcomp, precs, ferrsyme, slsyme, uc, um, &
                        dnsinf, dnssup, sigmsi, sigmss, &
                        sigmci, sigmcs, &
                        alpha, pivot, etat, ierr)
            kvarf = kVAR_NEW
        end if
    end if

end subroutine
