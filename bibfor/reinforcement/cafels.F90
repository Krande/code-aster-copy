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

subroutine cafels(cequi, effm, effn, ht, bw, &
                  enrobi, enrobs, scmaxi, scmaxs, ssmax, &
                  ferrcomp, precs, ferrsyme, slsyme, uc, um, &
                  dnsinf, dnssup, sigmsi, sigmss, &
                  sigmci, sigmcs, &
                  alpha, pivot, etat, ierr)

!_______________________________________________________________________________________________
!
!      CAFELS
!
!      CALCUL DES ACIERS EN FLEXION COMPOSEE A L'ELS CARACTERISTIQUE
!      CRITERE = LIMITATION DES CONTRAINTES ELASTIQUES
!
!      I CEQUI        COEFFICIENT D'EQUIVALENCE ACIER/BETON
!      I EFFM         MOMENT DE FLEXION
!      I EFFN         EFFORT NORMAL
!      I HT           HAUTEUR DE LA SECTION
!      I BW           LARGEUR DE LA SECTION
!      I ENROBI       ENROBAGE DES ARMATURES INFERIEURES
!      I ENROBS       ENROBAGE DES ARMATURES SUPERIEURES
!      I SCMAXI       CONTRAINTE DE COMPRESSION MAXI DU BETON EN FIBRE INF
!      I SCMAXS       CONTRAINTE DE COMPRESSION MAXI DU BETON EN FIBRE SUP
!      I SSMAX        CONTRAINTE MAXI DE L'ACIER DE FLEXION
!      I FERRCOMP     PRISE EN COMPTE DU FERRAILLAGE DE COMPRESSION
!                         FERRCOMP = 0 (NON)
!                         FERRCOMP = 1 (OUI)
!      I PRECS        PRECISION ITERATION
!      I FERRSYME     FERRAILLAGE SYMETRIQUE?
!                         FERRSYME = 0 (NON)
!                         FERRSYME = 1 (OUI)
!      I SLSYME       SECTION SEUIL DE TOLERANCE POUR UN FERRAILLAGE SYMETRIQUE
!      I UC           UNITE DES CONTRAINTES :
!                         UC = 0 CONTRAINTES EN Pa
!                         UC = 1 CONTRAINTES EN MPa
!      I UM           UNITE DES DIMENSIONS :
!                         UM = 0 DIMENSIONS EN m
!                         UM = 1 DIMENSIONS EN mm

!      O DNSINF       DENSITE DE L'ACIER INFERIEUR
!      O DNSSUP       DENSITE DE L'ACIER SUPERIEUR
!      O SIGMSI       CONTRAINTE AU NIVEAU DE L'ACIER INFERIEUR
!      O SIGMSS       CONTRAINTE AU NIVEAU DE L'ACIER SUPERIEUR
!      O SIGMCI       CONTRAINTE AU NIVEAU DE LA FIBRE SUPERIEURE DE BETON
!      O SIGMCS       CONTRAINTE AU NIVEAU DE LA FIBRE INFERIEURE DE BETON
!      O ALPHA        COEFFICIENT DE PROFONDEUR DE L'AN
!      O PIVOT        PIVOT DE FONCTIONNEMENT DE LA SECTION
!      O ETAT         ETAT DE FONCTIONNEMENT DE LA SECTION
!      O IERR         CODE RETOUR (0 = OK)
!_______________________________________________________________________________________________
!

    implicit none

#include "asterfort/cafelsiter.h"
!
!-----------------------------------------------------------------------
!!!!TERMES PRINCIPAUX D'ENTREE
!-----------------------------------------------------------------------
    real(kind=8) :: cequi
    real(kind=8) :: effm
    real(kind=8) :: effn
    real(kind=8) :: ht
    real(kind=8) :: bw
    real(kind=8) :: enrobi
    real(kind=8) :: enrobs
    real(kind=8) :: scmaxi
    real(kind=8) :: scmaxs
    real(kind=8) :: ssmax
    integer(kind=8) :: ferrcomp
    integer(kind=8) :: precs
    integer(kind=8) :: ferrsyme
    real(kind=8) :: slsyme
    integer(kind=8) :: uc
    integer(kind=8) :: um
    real(kind=8) :: dnsinf
    real(kind=8) :: dnssup
    real(kind=8) :: sigmsi
    real(kind=8) :: sigmss
    real(kind=8) :: sigmci
    real(kind=8) :: sigmcs
    real(kind=8) :: alpha
    integer(kind=8) :: pivot
    integer(kind=8) :: etat
    integer(kind=8) :: ierr

!-----------------------------------------------------------------------
!!!!VARIABLES DE CALCUL
!-----------------------------------------------------------------------
    real(kind=8) :: d, d0, a00, b00, c00, Del, var
    real(kind=8) :: Mcalc, Ncalc, scmax, scmaxc
    real(kind=8) :: unite_pa, mu, mu_12, alpha_12, mu_max, mu_lim, unite_m
    real(kind=8) :: AsTEND, AsCOMP, SsTEND, SsCOMP, ScTEND, ScCOMP
    real(kind=8) :: alpha_A, alpha_B, alpha_MED, RESIDU_A, RESIDU_B, RESIDU_MED, DIFF
    real(kind=8) :: X, Ncc, Mcc, Calc
    logical :: COND_ITER
    logical :: COND_NS

    !Significations des pointeurs :
    !PIVOT = 1 ==> "A"
    !PIVOT = 2 ==> "B"
    !PIVOT = 3 ==> "C"

    !ETAT = 1 ==> "BETON RESISTANT SEUL"
    !ETAT = 2 ==> "TRACTION PURE"
    !ETAT = 3 ==> "PARTIELLEMENT COMPRIMEE"
    !ETAT = 4 ==> "ENTIEREMENT TENDUE"
    !ETAT = 5 ==> "PARTIELLEMENT COMPRIMEE AVEC ACIER DE COMPRESSION"
    !ETAT = 6 ==> "ENTIEREMENT COMPRIMEE"

!-----------------------------------------------------------------------
!!!!LANCEMENT DU CALCUL
!-----------------------------------------------------------------------
!      INITIALISATION DU CODE RETOUR
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

    Mcalc = abs(effm)
    Ncalc = effn

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

    if (effm .ge. 0) then
        d = ht-enrobi
        d0 = ht-enrobs
        scmax = scmaxs
    else
        d = ht-enrobs
        d0 = ht-enrobi
        scmax = scmaxi
    end if
    scmaxc = min(scmaxs, scmaxi)

    alpha_12 = 1.0/(1.0+(ssmax/cequi)/scmax)
    mu_12 = 0.5*alpha_12*(1.0-alpha_12/3.0)
    mu_max = 0.5*1.0*(1.0-1.0/3.0)
    if (ferrcomp .eq. 1) then
        mu_lim = mu_12
    else
        mu_lim = mu_max
    end if

    mu = (Mcalc+Ncalc*(d-ht/2.0))/(d*d*bw*ScMAX)
    COND_ITER = .false.
    COND_NS = .false.

    if ((Mcalc .eq. 0) .AND. (Ncalc .ge. 0)) then

        ScCOMP = scmaxc
        ScTEND = scmaxc
        SsCOMP = scmaxc*cequi
        SsTEND = scmaxc*cequi
        Ncc = scmaxc*ht*bw
        Mcc = 0
        if (Ncc .ge. Ncalc) then
            etat = 1
            pivot = 2
            alpha = -1000
            ScCOMP = Ncalc/(ht*bw)
            ScTEND = Ncalc/(ht*bw)
            SsCOMP = 0
            SsTEND = 0
            AsTEND = 0
            AsCOMP = 0
        else
            AsTEND = 0.5*(Ncalc-Ncc)/SsTEND
            AsCOMP = 0.5*(Ncalc-Ncc)/SsCOMP
        end if

    elseif ((Mcalc .eq. 0) .AND. (Ncalc .lt. 0)) then

        ScCOMP = -ssmax/cequi
        ScTEND = -ssmax/cequi
        SsCOMP = -ssmax
        SsTEND = -ssmax
        alpha = -1000
        pivot = 1
        etat = 2
        AsTEND = 0.5*Ncalc/SsTEND
        AsCOMP = 0.5*Ncalc/SsCOMP

    elseif ((mu .gt. 0) .AND. (mu .le. mu_lim)) then

        etat = 3

        if (mu .le. mu_12) then

            pivot = 1
            !Calc pour alpha = 0
            alpha_A = 0
            alpha = alpha_A
            RESIDU_A = Mcalc+Ncalc*(d-0.5*ht)
            var = bw*(d**2)*(ssmax/cequi)*0.5*(alpha**2)*((1.-alpha/3.)/(1.-alpha))
            RESIDU_A = RESIDU_A-var
            !Calc pour alpha = alpha_12
            alpha_B = alpha_12
            alpha = alpha_B
            RESIDU_B = Mcalc+Ncalc*(d-0.5*ht)
            var = bw*(d**2)*(ssmax/cequi)*0.5*(alpha**2)*((1.-alpha/3.)/(1.-alpha))
            RESIDU_B = RESIDU_B-var
            !Iteration
            DIFF = alpha_12
            do while (DIFF .gt. 0.01)
                alpha_MED = 0.5*(alpha_A+alpha_B)
                alpha = alpha_MED
                RESIDU_MED = Mcalc+Ncalc*(d-0.5*ht)
                var = bw*(d**2)*(ssmax/cequi)*0.5*(alpha**2)*((1.-alpha/3.)/(1.-alpha))
                RESIDU_MED = RESIDU_MED-var
                if ((RESIDU_MED*RESIDU_B) .gt. 0) then
                    alpha_B = alpha_MED
                    RESIDU_B = RESIDU_MED
                else
                    alpha_A = alpha_MED
                    RESIDU_A = RESIDU_MED
                end if
                DIFF = abs(alpha_B-alpha_A)
            end do
            alpha = 0.5*(alpha_A+alpha_B)
            X = alpha*d
            SsTEND = -ssmax
            ScCOMP = (X/(d-X))*(ssmax/cequi)
            ScTEND = ScCOMP*(1-ht/X)
            SsCOMP = ScCOMP*(1-(ht-d0)/X)*cequi
            AsTEND = (0.5*ScCOMP*X*bw-Ncalc)/(-SsTEND)
            AsCOMP = 0

            if (AsTEND .lt. 0) then
                if (abs(AsTEND) .lt. ((1.e-5)*(1.e6/(unite_m**2)))) then
                    AsTEND = 0
                else
                    COND_ITER = .true.
                    goto 998
                end if
            else
                COND_NS = .true.
            end if

            if (ferrsyme .eq. 1) then
                Calc = abs(AsCOMP-AsTEND)
                if (Calc .gt. slsyme) then
                    COND_ITER = .true.
                    goto 998
                end if
            end if

        else

            pivot = 2
            a00 = (1./2.)*(-1./3.)
            b00 = 1./2.
            c00 = -(Mcalc+Ncalc*(d-0.5*ht))/(d*d*bw*scmax)
            Del = b00*b00-4.0*a00*c00
            if (Del .ge. 0) Then
                alpha_A = (-b00+(Del)**(0.5))/(2.0*a00)
                alpha_B = (-b00-(Del)**(0.5))/(2.0*a00)
                if ((alpha_A .gt. 0) .And. (alpha_A .le. 1)) Then
                    alpha = alpha_A
                elseif ((alpha_B .gt. 0) .And. (alpha_B .le. 1)) Then
                    alpha = alpha_B
                else
                    COND_ITER = .True.
                end if
            else
                COND_ITER = .True.
            end if
            if (COND_ITER .eqv. (.False.)) Then
                X = alpha*d
                ScCOMP = scmax
                SsTEND = -cequi*scmax*(d-X)/X
                ScTEND = ScCOMP*(1-ht/X)
                SsCOMP = ScCOMP*(1-(ht-d0)/X)*cequi
                AsTEND = (0.5*ScCOMP*X*bw-Ncalc)/(-SsTEND)
                AsCOMP = 0
            end if

            if (AsTEND .lt. 0) then
                if (abs(AsTEND) .lt. ((1.e-5)*(1.e6/(unite_m**2)))) then
                    AsTEND = 0
                else
                    COND_ITER = .true.
                    goto 998
                end if
            else
                COND_NS = .true.
            end if

            if (ferrsyme .eq. 1) then
                Calc = abs(AsCOMP-AsTEND)
                if (Calc .gt. slsyme) then
                    COND_ITER = .true.
                    goto 998
                end if
            end if

        end if

    else
        COND_ITER = .TRUE.

    end if

998 continue

    if (COND_ITER .eqv. (.TRUE.)) then
        call cafelsiter(cequi, effm, effn, ht, bw, &
                        enrobi, enrobs, scmaxi, scmaxs, ssmax, ferrcomp, &
                        precs, ferrsyme, slsyme, uc, um, COND_NS, &
                        AsTEND, AsCOMP, SsTEND, SsCOMP, &
                        ScTEND, ScCOMP, &
                        alpha, pivot, etat, ierr)
    end if

!------------------------------------------------------------------
!Distinction Ferr Sup et Inf
!------------------------------------------------------------------

    if (effm .ge. 0) then
        dnssup = AsCOMP
        dnsinf = AsTEND
        sigmss = SsCOMP
        sigmsi = SsTEND
        sigmcs = ScCOMP
        sigmci = ScTEND
    else
        dnssup = AsTEND
        dnsinf = AsCOMP
        sigmss = SsTEND
        sigmsi = SsCOMP
        sigmcs = ScTEND
        sigmci = ScCOMP
    end if

    if (ferrcomp .eq. 0) then
        if (((sigmss .gt. 0) .and. (dnssup .gt. 0)) &
              & .or. ((sigmsi .gt. 0) .and. (dnsinf .gt. 0))) then
            etat = 0
            pivot = 0
            dnssup = -1
            dnsinf = -1
            sigmss = -1
            sigmsi = -1
            sigmcs = -1
            sigmci = -1
            ierr = 1
        end if
    end if

end subroutine
