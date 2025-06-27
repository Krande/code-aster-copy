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

subroutine cafelu(typco, alphacc, effm, effn, ht, bw, &
                  enrobi, enrobs, facier, fbeton, gammas, gammac, &
                  clacier, eys, typdiag, ferrcomp, precs, ferrsyme, slsyme, uc, um, &
                  dnsinf, dnssup, sigmsi, sigmss, ecinf, ecsup, &
                  alpha, pivot, etat, ierr)
!______________________________________________________________________
!
!      CAFELU

!      CALCUL DES ACIERS EN FLEXION COMPOSEE A L'ELU
!      CRITERE = LIMITATION DES DEFORMATIONS
!
!      I TYPCO     CODIFICATION UTILISEE (1 = BAEL91, 2 = EC2)
!      I ALPHACC   COEFFICIENT DE SECURITE SUR LA RESISTANCE
!                  DE CALCUL DU BETON EN COMPRESSION
!      I EFFM      MOMENT DE FLEXION
!      I EFFN      EFFORT NORMAL
!      I HT        HAUTEUR DE LA SECTION
!      I BW        LARGEUR DE LA SECTION
!      I ENROBI    ENROBAGE DES ARMATURES INFERIEURES
!      I ENROBS    ENROBAGE DES ARMATURES SUPERIEURES
!      I FACIER    LIMITE D'ELASTICITE DES ACIERS (CONTRAINTE)
!      I FBETON    RESISTANCE EN COMPRESSION DU BETON (CONTRAINTE)
!      I GAMMAS    COEFFICIENT DE SECURITE SUR LA RESISTANCE
!                  DE CALCUL DES ACIERS
!      I GAMMAC    COEFFICIENT DE SECURITE SUR LA RESISTANCE
!                  DE CALCUL DU BETON
!      I CLACIER   CLASSE DE DUCTILITE DES ACIERS (UTILISE POUR EC2) :
!                     CLACIER = 0 ACIER PEU DUCTILE (CLASSE A)
!                     CLACIER = 1 ACIER MOYENNEMENT DUCTILE (CLASSE B)
!                     CLACIER = 3 ACIER FORTEMENT DUCTILE (CLASSE C)
!      I EYS       MODULE D'YOUNG DE L'ACIER
!      I TYPDIAG   TYPE DE DIAGRAMME UTILISÉ POUR L'ACIER
!                     TYPDIAG = 1 ("B1" ==> PALIER INCLINÉ)
!                     TYPDIAG = 2 ("B2" ==> PALIER HORIZONTAL)
!      I FERRCOMP  PRISE EN COMPTE DU FERRAILLAGE DE COMPRESSION
!                     FERRCOMP = 0 (NON)
!                     FERRCOMP = 1 (OUI)
!      I PRECS     PRECISION ITERATION
!      I FERRSYME  FERRAILLAGE SYMETRIQUE?
!                     FERRSYME = 0 (NON)
!                     FERRSYME = 1 (OUI)
!      I SLSYME    SECTION SEUIL DE TOLERANCE POUR UN FERRAILLAGE SYMETRIQUE
!      I UC        UNITE DES CONTRAINTES :
!                     UC = 0 CONTRAINTES EN Pa
!                     UC = 1 CONTRAINTES EN MPa
!      I UM        UNITE DES DIMENSIONS :
!                     UM = 0 DIMENSIONS EN m
!                     UM = 1 DIMENSIONS EN mm
!
!      O DNSINF    DENSITE DE L'ACIER INFERIEUR
!      O DNSSUP    DENSITE DE L'ACIER SUPERIEUR
!      O SIGMSI    CONTRAINTE AU NIVEAU DE L'ACIER INFERIEUR
!      O SIGMSS    CONTRAINTE AU NIVEAU DE L'ACIER SUPERIEUR
!      O ECINF     DEFORMATION DU BETON EN FIBRE INFÉRIEURE
!      O ECSUP     DEFORMATION DU BETON EN FIBRE SUPERIEURE
!      O ALPHA     COEFFICIENT DE PROFONDEUR DE L'AN
!      O PIVOT     PIVOT DE FONCTIONNEMENT DE LA SECTION
!      O ETAT      ETAT DE FONCTIONNEMENT DE LA SECTION
!      O IERR      CODE RETOUR (0 = OK)
!
!______________________________________________________________________
!
    implicit none

#include "asterfort/cafeluiter.h"

    integer(kind=8) :: typco
    real(kind=8) :: alphacc
    real(kind=8) :: effm
    real(kind=8) :: effn
    real(kind=8) :: ht
    real(kind=8) :: bw
    real(kind=8) :: enrobi
    real(kind=8) :: enrobs
    real(kind=8) :: facier
    real(kind=8) :: fbeton
    real(kind=8) :: gammas
    real(kind=8) :: gammac
    integer(kind=8) :: clacier
    real(kind=8) :: eys
    integer(kind=8) :: typdiag
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
    real(kind=8) :: ecinf
    real(kind=8) :: ecsup
    real(kind=8) :: alpha
    integer(kind=8) :: pivot
    integer(kind=8) :: etat
    integer(kind=8) :: ierr

!-----------------------------------------------------------------------
!!!!VARIABLES DE CALCUL
!-----------------------------------------------------------------------

    real(kind=8) :: enrob, d, d0, fctm, fctd
    real(kind=8) :: unite_pa, unite_m
    real(kind=8) :: fyd, fcd, nC, ktys, xC, yC, xCt, D00, m1, m2
    real(kind=8) :: Mu, MccMAX, NccMAX, NccMIN, MuC, NuC
    real(kind=8) :: MuBC, MuBC_sup, MuBC_inf, MuAB, MuR, MuLIM
    real(kind=8) :: eta, lambda, Esu, Euk, Ecu, Ec2, Ese
    real(kind=8) :: Xsup
    real(kind=8) :: piv_a, piv_b, piv_c, alphaAB, alphaR, alphaBC
    real(kind=8) :: V, COEF, VAR_COEF
    logical :: COND_ITER, COND_NS
    real(kind=8) :: Calc, DELTA
    real(kind=8) :: Ec, EcTEND, EsCOMP, EsTEND
    real(kind=8) :: SigmAsCOMP, SigmAsTEND, AsCOMP, AsTEND
    real(kind=8) :: a00, b00, c00, alpha_1, alpha_2

    !Significations des pointeurs :
    !PIVOT = 0 ==> "?"
    !PIVOT = 1 ==> "A"
    !PIVOT = 2 ==> "B"
    !PIVOT = 3 ==> "C"

    !ETAT = 0 ==> "EQUILIBRE IMPOSSIBLE"
    !ETAT = 1 ==> "BETON RESISTANT SEUL"
    !ETAT = 2 ==> "TRACTION PURE"
    !ETAT = 3 ==> "PARTIELLEMENT COMPRIMEE"
    !ETAT = 4 ==> "ENTIEREMENT TENDUE"
    !ETAT = 5 ==> "PARTIELLEMENT COMPRIMEE AVEC ACIER DE COMPRESSION"
    !ETAT = 6 ==> "ENTIEREMENT COMPRIMEE"
    !ETAT = 7 ==> "COMPRESSION PURE"

!   INITIALISATION DU CODE RETOUR
    ierr = 0
    etat = 1
    pivot = 0
    alpha = -1000
    dnssup = -1
    dnsinf = -1
    sigmss = -1
    sigmsi = -1
    ecsup = -1
    ecinf = -1

!   CONVENTION D'ORIENTATION

    if (effm .ge. 0.) then
        enrob = enrobi
        d = ht-enrob
        d0 = enrobs
    else
        enrob = enrobs
        d = ht-enrob
        d0 = enrobi
    end if

    if (typco .eq. 1) then
!       CALCUL DES PARAMETRES POUR CODIFICATION = 'BAEL91'

        eta = 1.d0
        lambda = 0.8
        piv_a = 10.0E-3
        piv_b = 3.5E-3
        piv_c = 2.0E-3
        nC = 2
        fyd = facier/gammas
        fcd = fbeton*alphacc/gammac

    else if (typco .eq. 2) then
!       CALCUL DES PARAMETRES POUR CODIFICATION = 'EC2'

        if (uc .eq. 0) then
            unite_pa = 1.e-6
        elseif (uc .eq. 1) then
            unite_pa = 1.
        end if
        eta = min(1.d0, 1.d0-(fbeton*unite_pa-50.d0)/200.d0)
        lambda = min(0.8, 0.8-(fbeton*unite_pa-50.d0)/400.d0)
        if (clacier .eq. 0) then
            piv_a = 0.9*2.5e-2
            ktys = 1.05
        else if (clacier .eq. 1) then
            piv_a = 0.9*5.e-2
            ktys = 1.08
        else
            piv_a = 0.9*7.5e-2
            ktys = 1.15
        end if
        piv_b = min(3.5E-3, 0.26*0.01+3.5*0.01*(((90.d0-fbeton*unite_pa)/100.d0)**4))
        piv_c = 2.0E-3
        if ((fbeton*unite_pa) .ge. (50.d0)) then
            piv_c = 0.2*0.01+0.0085*0.01*((fbeton*unite_pa-50.d0)**(0.53))
        end if
        nC = min(2.0, 1.4+23.4*(((90.d0-fbeton*unite_pa)/100.d0)**4))
        fyd = facier/gammas
        fcd = fbeton*alphacc/gammac
        if ((fbeton*unite_pa) .le. 50) then
            fctm = 0.7*0.3*(fbeton**(2./3.))
        else
            fctm = 0.7*2.12*log(1+0.1*(fbeton+8))
        end if
        fctd = fctm/gammac

    end if

    if (um .eq. 0) then
        unite_m = 1.e3
    elseif (um .eq. 1) then
        unite_m = 1.
    end if

!   Paramètres de calcul
    Xsup = piv_b/piv_c
    xC = (1-piv_c/piv_b)*ht
    yC = ht-xC
    xCt = xC/ht
    D00 = (ht-d0)/ht
    m1 = (((1-xCt)**(nC+1))/(2.d0*(nC+1)))*(1-(2.d0*(1-xCt))/(nC+2))
    m2 = -((1-xCt)**(nC+1))/(nC+1)
    MuC = abs(effm)/(bw*(ht**2)*fcd)
    NuC = effn/(bw*ht*fcd)
    NccMAX = (bw*ht*fcd)
    NccMIN = (1+m2*((Xsup)**(nC)))*bw*ht*fcd
    MccMAX = m1*((Xsup)**(nC))*bw*(ht**2)*fcd
    Mu = (abs(effm)+effn*(d-0.5*ht))/(bw*(d**2)*eta*fcd)
    Esu = piv_a
    Euk = Esu/0.9
    Ecu = piv_b
    Ec2 = piv_c
    Ese = fyd/eys
    alphaAB = 1./(1+Esu/Ecu)
    alphaR = 1./(1+Ese/Ecu)
    alphaBC = 1.
    MuAB = lambda*alphaAB*(1-0.5*lambda*alphaAB)
    MuR = lambda*alphaR*(1-0.5*lambda*alphaR)
    MuBC_inf = lambda*alphaBC*(1-0.5*lambda*alphaBC)
    if (ferrcomp .eq. 0) then
        MuLIM = MuBC_inf
    elseif (ferrcomp .eq. 1) then
        MuLIM = MuR
    end if
    V = Ecu/Ec2
    COEF = (1-(V/ht)*yC)
    if (abs(COEF) .gt. epsilon(COEF)) then
        VAR_COEF = Abs(COEF)/COEF
    else
        VAR_COEF = 1
    end if
    MuBC_sup = xC*(d-0.5*xC)+0.5*yC*yC+(d-ht)*yC &
               & +(ht/V)*(1/(nC+1))*yC*(VAR_COEF*(Abs(COEF))**(nC+1)) &
               & +(1./((nC+1)*(nC+2)))*((ht/V)**2)*(VAR_COEF*(Abs(COEF))**(nC+2)-1) &
               & +(d-ht)*(1./(nC+1))*(ht/V)*(VAR_COEF*(Abs(COEF))**(nC+1)-1)
    MuBC_sup = MuBC_sup/(d*d*Eta)
    MuBC = Max(MuBC_inf, MuBC_sup)
    COND_ITER = .false.
    COND_NS = .false.

!!!Traitement des CAS TRIVIAUX (PIVOTS A et B sans armatures de compression)
!!!--------------------------------------------------------------------------

    if ((abs(effm) .lt. epsilon(effm)) .and. (effn .ge. 0)) then

        etat = 7
        pivot = 3
        Ec = Ec2
        EcTEND = Ec2
        EsTEND = Ec2
        EsCOMP = Ec2
        alpha = -1
        if (typdiag .eq. 1) then
        if (Abs(EsTEND) .lt. Ese) then
            SigmAsTEND = eys*(Abs(EsTEND))
        else
            SigmAsTEND = fyd+((ktys*fyd-fyd)/(Euk-Ese))*(Abs(EsTEND)-Ese)
        end if
        if (Abs(EsCOMP) .lt. Ese) then
            SigmAsCOMP = eys*(Abs(EsCOMP))
        else
            SigmAsCOMP = fyd+((ktys*fyd-fyd)/(Euk-Ese))*(Abs(EsCOMP)-Ese)
        end if
        else
        if (Abs(EsTEND) .lt. Ese) then
            SigmAsTEND = eys*(Abs(EsTEND))
        else
            SigmAsTEND = fyd
        end if
        if (Abs(EsCOMP) .lt. Ese) then
            SigmAsCOMP = eys*(Abs(EsCOMP))
        else
            SigmAsCOMP = fyd
        end if
        end if

        if (EsTEND .lt. 0) then
            SigmAsTEND = -SigmAsTEND
        end if
        if (EsCOMP .lt. 0) then
            SigmAsCOMP = -SigmAsCOMP
        end if

        if (effn .le. NccMAX) then
            etat = 1
            pivot = 3
            alpha = -1000
            Ec = 0
            EcTEND = -1
            EsTEND = -1
            EsCOMP = -1
            SigmAsCOMP = -1
            SigmAsTEND = -1
            AsCOMP = 0
            AsTEND = 0
        else
            AsTEND = (0.5*(effn-NccMAX))/SigmAsTEND
            AsCOMP = (0.5*(effn-NccMAX))/SigmAsCOMP
        end if

    elseif ((abs(effm) .lt. epsilon(effm)) .and. (effn .lt. 0)) then

        etat = 2
        pivot = 1
        Ec = -Esu
        EcTEND = -Esu
        EsTEND = -Esu
        EsCOMP = -Esu

        if (typdiag .eq. 1) then
        if (Abs(EsTEND) .lt. Ese) then
            SigmAsTEND = eys*(Abs(EsTEND))
        else
            SigmAsTEND = fyd+((ktys*fyd-fyd)/(Euk-Ese))*(Abs(EsTEND)-Ese)
        end if
        if (Abs(EsCOMP) .lt. Ese) Then
            SigmAsCOMP = eys*(Abs(EsCOMP))
        else
            SigmAsCOMP = fyd+((ktys*fyd-fyd)/(Euk-Ese))*(Abs(EsCOMP)-Ese)
        end if
        else
        if (Abs(EsTEND) .lt. Ese) then
            SigmAsTEND = eys*(Abs(EsTEND))
        else
            SigmAsTEND = fyd
        end if
        if (Abs(EsCOMP) .lt. Ese) then
            SigmAsCOMP = eys*(Abs(EsCOMP))
        else
            SigmAsCOMP = fyd
        end if
        end if

        if (EsTEND .lt. 0) then
            SigmAsTEND = -SigmAsTEND
        end if
        if (EsCOMP .lt. 0) then
            SigmAsCOMP = -SigmAsCOMP
        end if

        AsTEND = -(0.5*effn)/Abs(SigmAsTEND)
        AsCOMP = -(0.5*effn)/Abs(SigmAsCOMP)

    elseif ((Mu .gt. 0) .and. (Mu .le. MuLIM)) then

        etat = 3
        a00 = lambda*lambda*0.5
        b00 = -lambda
        c00 = Mu
        DELTA = b00*b00-4*a00*c00

        if (DELTA .ge. 0) then
            alpha_1 = (-b00+(DELTA)**0.5)/(2*a00)
            alpha_2 = (-b00-(DELTA)**0.5)/(2*a00)
            if ((alpha_1 .ge. 0) .and. (alpha_1 .le. 1)) then
                alpha = alpha_1
            elseif ((alpha_2 .ge. 0) .and. (alpha_2 .le. 1)) then
                alpha = alpha_2
            else
                COND_ITER = .true.
                goto 998
            end if
        else
            COND_ITER = .true.
            goto 998
        end if

        if (alpha .lt. alphaAB) then
            pivot = 1
            EsTEND = -Esu
            Ec = Esu*alpha/(1-alpha)
        elseif (alpha .gt. alphaAB) Then
            pivot = 2
            Ec = Ecu
            EsTEND = -Ecu*(1-alpha)/alpha
        else
            pivot = 1
            Ec = Ecu
            EsTEND = -Esu
        end if

        EsCOMP = ((EsTEND-Ec)/d)*(d0)+Ec
        EcTEND = ((EsTEND-Ec)/d)*(ht)+Ec

        if (typdiag .eq. 1) then
        if (Abs(EsTEND) .lt. Ese) then
            SigmAsTEND = eys*(Abs(EsTEND))
        else
            SigmAsTEND = fyd+((ktys*fyd-fyd)/(Euk-Ese))*(Abs(EsTEND)-Ese)
        end if
        if (Abs(EsCOMP) .lt. Ese) Then
            SigmAsCOMP = eys*(Abs(EsCOMP))
        else
            SigmAsCOMP = fyd+((ktys*fyd-fyd)/(Euk-Ese))*(Abs(EsCOMP)-Ese)
        end if
        else
        if (Abs(EsTEND) .lt. Ese) then
            SigmAsTEND = eys*(Abs(EsTEND))
        else
            SigmAsTEND = fyd
        end if
        if (Abs(EsCOMP) .lt. Ese) then
            SigmAsCOMP = eys*(Abs(EsCOMP))
        else
            SigmAsCOMP = fyd
        end if
        end if

        if (EsTEND .lt. 0) then
            SigmAsTEND = -SigmAsTEND
        end if
        if (EsCOMP .lt. 0) then
            SigmAsCOMP = -SigmAsCOMP
        end if

        AsCOMP = 0
        AsTEND = (lambda*Eta*bw*fcd*alpha*d-effn)/Abs(SigmAsTEND)

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

        COND_ITER = .true.

    end if

998 Continue

    !Pour les cas où COND_ITER = TRUE
    If (COND_ITER .eqv. (.true.)) Then
        call cafeluiter(typco, alphacc, effm, effn, ht, bw, &
                        enrobi, enrobs, facier, fbeton, gammas, gammac, &
                        clacier, eys, typdiag, ferrcomp, precs, ferrsyme, slsyme, uc, um, &
                        COND_NS, AsTEND, AsCOMP, SigmAsTEND, SigmAsCOMP, EcTEND, Ec, &
                        alpha, pivot, etat, ierr)
    end if

!------------------------------------------------------------------
!Distinction Ferr Sup et Inf ET Prise en compte COMPRESSION
!------------------------------------------------------------------

    if (effm .gt. 0) then
        dnssup = AsCOMP
        dnsinf = AsTEND
        sigmss = SigmAsCOMP
        sigmsi = SigmAsTEND
        ecsup = Ec
        ecinf = EcTEND
    else
        dnssup = AsTEND
        dnsinf = AsCOMP
        sigmss = SigmAsTEND
        sigmsi = SigmAsCOMP
        ecsup = EcTEND
        ecinf = Ec
    end if

    if (ferrcomp .eq. 0) then
        if (((sigmss .gt. 0) .and. (dnssup .gt. 0)) &
              & .or. ((sigmsi .gt. 0) .and. (dnsinf .gt. 0))) then
            etat = 0
            pivot = 0
            dnssup = -1
            dnsinf = -1
            sigmsi = -1
            sigmss = -1
            ecinf = -1
            ecsup = -1
            ierr = 1
        end if
    end if

end subroutine
