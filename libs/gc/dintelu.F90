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

subroutine dintelu(typco, alphacc, ht, bw, enrobi, enrobs, facier, fbeton, &
                   gammas, gammac, clacier, eys, typdiag, uc, &
                   ntot, dnsinf, dnssup, nrd, mrd, ndemi) bind(C)
!______________________________________________________________________
!
!      DINTELU

!      CONSTRUCTION DU DIAGRAMME D'INTERACTION D'UNE SECTION
!      FERRAILLEE - VERIFICATION D'UN FERRAILLAGE EXISTANT
!      CRITERE = LIMITATION DES DEFORMATIONS
!
!      I TYPCO     CODIFICATION UTILISEE (1 = BAEL91, 2 = EC2)
!      I ALPHACC   COEFFICIENT DE SECURITE SUR LA RESISTANCE
!                  DE CALCUL DU BETON EN SUPRESSION
!      I HT        HAUTEUR DE LA SECTION
!      I BW        LARGEUR DE LA SECTION
!      I ENROBI    ENROBAGE DES ARMATURES INFERIEURES
!      I ENROBS    ENROBAGE DES ARMATURES SUPERIEURES
!      I FACIER    LIMITE D'ELASTICITE DES ACIERS (CONTRAINTE)
!      I FBETON    RESISTANCE EN SUPRESSION DU BETON (CONTRAINTE)
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
!      I UC        UNITE DES CONTRAINTES :
!                     UC = 0 CONTRAINTES EN Pa
!                     UC = 1 CONTRAINTES EN MPa
!      I DNSINF    DENSITE DE L'ACIER INFERIEUR
!      I DNSSUP    DENSITE DE L'ACIER SUPERIEUR
!
!      I NTOT      DIMENSIONS DES VECTEURS
!      O NRD       VECTEUR DES EFFORTS NORMAUX RESISTANTS (DIAG INTER)
!      O MRD       VECTEUR DES MOMENTS RESISTANTS (DIAG INTER)
!
!______________________________________________________________________
!
!
    use, intrinsic :: iso_c_binding
    implicit none
!
    integer(c_long), intent(in) :: typco
    real(c_double), intent(in) :: alphacc
    real(c_double), intent(in) :: ht
    real(c_double), intent(in) :: bw
    real(c_double), intent(in) :: enrobi
    real(c_double), intent(in) :: enrobs
    real(c_double), intent(in) :: facier
    real(c_double), intent(in) :: fbeton
    real(c_double), intent(in) :: gammas
    real(c_double), intent(in) :: gammac
    integer(c_long), intent(in) :: clacier
    real(c_double), intent(in) :: eys
    integer(c_long), intent(in) :: typdiag
    integer(c_long), intent(in) :: uc
    integer(c_long), intent(inout) :: ntot
    real(c_double), intent(in), optional :: dnsinf
    real(c_double), intent(in), optional :: dnssup
    real(c_double), intent(out), optional :: nrd(1:ntot)
    real(c_double), intent(out), optional :: mrd(1:ntot)
    integer(c_long), intent(out), optional :: ndemi

!-----------------------------------------------------------------------
!!!!VARIABLES DE CALCUL
!-----------------------------------------------------------------------

    real(kind=8) :: d, d0, dneg, d0neg
    real(kind=8) :: unite_pa, DE, X, Calc
    real(kind=8) :: fyd, fcd, nC, ktys, xC, yC, xCt, m1, m2
    real(kind=8) :: Esu, Euk, Ecu, Ec2, Ese, Xsup
    real(kind=8) :: piv_a, piv_b, piv_c, alpha, alphaAB, alphaR, alphaBC
    real(kind=8) :: COEF1, COEF2, VAR_COEF1, VAR_COEF2
    real(kind=8) :: DELTA, x1, y1, Beta, yE
    real(kind=8) :: EcINF, EcSUP, EsSUP, EsINF
    real(kind=8) :: SigmAsSUP, SigmAsINF, Ncc, Mcc
    integer :: N_ET, N_PC, N_EC, N_PCN, k

    real(kind=8), allocatable :: N_P1(:), M_P1(:)
    real(kind=8), allocatable :: N_P2(:), M_P2(:)
    real(kind=8), allocatable :: N_P3(:), M_P3(:)
    real(kind=8), allocatable :: N_P4(:), M_P4(:)
    real(kind=8), allocatable :: N_P5(:), M_P5(:)
    real(kind=8), allocatable :: N_P6(:), M_P6(:)

    nrd = 0.0
    mrd = 0.0

    k = 1

    if (typco .eq. 1) then
!       CALCUL DES PARAMETRES POUR CODIFICATION = 'BAEL91'

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

    end if

    Esu = piv_a
    Euk = Esu/0.9
    Ecu = piv_b
    Ec2 = piv_c
    Ese = fyd/eys
    alphaAB = 1./(1+Esu/Ecu)
    alphaR = 1./(1+Ese/Ecu)
    alphaBC = 1.
    d = ht-enrobi
    d0 = enrobs
    dneg = ht-enrobs
    d0neg = enrobi

!   Paramètres de calcul
    Xsup = piv_b/piv_c
    xC = (1-piv_c/piv_b)*ht
    yC = ht-xC
    xCt = xC/ht
    m1 = (((1-xCt)**(nC+1))/(2.d0*(nC+1)))*(1-(2.d0*(1-xCt))/(nC+2))
    m2 = -((1-xCt)**(nC+1))/(nC+1)

    N_ET = floor(Esu*1000)+1
    N_PC = ceiling((ht/d)*100)+1
    N_EC = ceiling(Xsup*100)+1
    N_PCN = ceiling((ht/dneg)*100)+1
    if (ntot < 0) then
        ntot = N_ET+N_ET+N_PC+N_EC+N_EC+N_PCN
        if (present(ndemi)) then
            ndemi = N_ET+N_PC+N_EC
        end if
        return
    end if
    if (.not. present(dnsinf) .or. .not. present(dnssup)) then
        write (6, *) "SyntaxError: dnsinf and dnssup are required"
        return
    end if
    if (.not. present(nrd) .or. .not. present(mrd)) then
        write (6, *) "SyntaxError: nrd and mrd are required"
        return
    end if

!-----------------------------------------------------------------------
!Traitement des différents cas (Pivots A / B / C)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!ET = Entièrement Tendue
!PC = Partiellement Comprimée
!EC = Entièrement Comprimée

!Traitement en PIVOT A - Entièrement Tendu (ET) + Moment Positif
!---------------------------------------------------------------

    allocate (N_P1(N_ET))
    allocate (M_P1(N_ET))

    do k = 1, N_ET

        if (k .eq. 1) then
            EcSUP = -Esu
        else
            EcSUP = -(1.e-3)*(N_ET-k)
        end if

        EsINF = -Esu
        EsSUP = ((EsINF-EcSUP)/d)*(d0)+EcSUP
        EcINF = ((EsINF-EcSUP)/d)*(ht)+EcSUP

        if (typdiag .eq. 1) then
        if (Abs(EsINF) .lt. Ese) then
            SigmAsINF = eys*(Abs(EsINF))
        else
            SigmAsINF = fyd+((ktys*fyd-fyd)/(Euk-Ese))*(Abs(EsINF)-Ese)
        end if
        if (Abs(EsSUP) .lt. Ese) Then
            SigmAsSUP = eys*(Abs(EsSUP))
        else
            SigmAsSUP = fyd+((ktys*fyd-fyd)/(Euk-Ese))*(Abs(EsSUP)-Ese)
        end if
        else
        if (Abs(EsINF) .lt. Ese) then
            SigmAsINF = eys*(Abs(EsINF))
        else
            SigmAsINF = fyd
        end if
        if (Abs(EsSUP) .lt. Ese) then
            SigmAsSUP = eys*(Abs(EsSUP))
        else
            SigmAsSUP = fyd
        end if
        end if

        if (EsINF .lt. 0) then
            SigmAsINF = -SigmAsINF
        end if
        if (EsSUP .lt. 0) then
            SigmAsSUP = -SigmAsSUP
        end if

        N_P1(k) = dnsinf*SigmAsINF+dnssup*SigmAsSUP
        M_P1(k) = -dnsinf*SigmAsINF*(d-0.5*ht)+dnssup*SigmAsSUP*(0.5*ht-d0)

    end do

!
!Traitement en PIVOT A et B - Partiellement Comprimée (PC) + Moment Positif
!--------------------------------------------------------------------------

    allocate (N_P2(N_PC))
    allocate (M_P2(N_PC))

    do k = 1, N_PC

        if (k .lt. N_PC) then
            alpha = (k-1)*0.01
        else
            alpha = ht/d
        end if

        if (alpha .lt. alphaAB) then
            EsINF = -Esu
            EcSUP = Esu*alpha/(1-alpha)
        else
            EcSUP = Ecu
            EsINF = -Ecu*(1-alpha)/alpha
        end if

        EsSUP = ((EsINF-EcSUP)/d)*(d0)+EcSUP
        EcINF = ((EsINF-EcSUP)/d)*(ht)+EcSUP

        if (typdiag .eq. 1) then
        if (Abs(EsINF) .lt. Ese) then
            SigmAsINF = eys*(Abs(EsINF))
        else
            SigmAsINF = fyd+((ktys*fyd-fyd)/(Euk-Ese))*(Abs(EsINF)-Ese)
        end if
        if (Abs(EsSUP) .lt. Ese) Then
            SigmAsSUP = eys*(Abs(EsSUP))
        else
            SigmAsSUP = fyd+((ktys*fyd-fyd)/(Euk-Ese))*(Abs(EsSUP)-Ese)
        end if
        else
        if (Abs(EsINF) .lt. Ese) then
            SigmAsINF = eys*(Abs(EsINF))
        else
            SigmAsINF = fyd
        end if
        if (Abs(EsSUP) .lt. Ese) then
            SigmAsSUP = eys*(Abs(EsSUP))
        else
            SigmAsSUP = fyd
        end if
        end if

        if (EsINF .lt. 0) then
            SigmAsINF = -SigmAsINF
        end if
        if (EsSUP .lt. 0) then
            SigmAsSUP = -SigmAsSUP
        end if

        x1 = (EsINF-EcSUP)/Ec2
        y1 = EcSUP/Ec2
        DELTA = d/ht

        if (EcSUP .le. Ec2) then
            Beta = 0
        else
            yE = ((Ec2-EcSUP)/(EsINF-EcSUP))*d
            Beta = yE/d
        end if

        COEF1 = (1-y1-alpha*x1)
        COEF2 = (1-y1-Beta*x1)
        if (abs(COEF1) .gt. epsilon(COEF1)) then
            VAR_COEF1 = Abs(COEF1)/COEF1
        else
            VAR_COEF1 = 1
        end if
        if (abs(COEF2) .gt. epsilon(COEF2)) then
            VAR_COEF2 = Abs(COEF2)/COEF2
        else
            VAR_COEF2 = 1
        end if

        Ncc = (fcd*bw*d)*(alpha+(1/((nC+1)*x1))*(VAR_COEF1*((Abs(COEF1))**(nC+1)) &
               & -VAR_COEF2*((Abs(COEF2))**(nC+1))))
        Mcc = (bw*d*d*fcd)*(0.5*alpha*(1/DELTA-alpha) &
               & +(1/(2*DELTA))*(1/((nC+1)*x1))*(VAR_COEF1*((Abs(COEF1))**(nC+1)) &
               & -VAR_COEF2*((Abs(COEF2))**(nC+1))) &
               & -(1/((nC+1)*x1))*(alpha*VAR_COEF1*((Abs(COEF1))**(nC+1)) &
               & -Beta*VAR_COEF2*((Abs(COEF2))**(nC+1))) &
               & -(1/((nC+1)*(nC+2)*x1*x1))*(VAR_COEF1*((Abs(COEF1))**(nC+2)) &
               & -VAR_COEF2*((Abs(COEF2))**(nC+2))))

        N_P2(k) = dnsinf*SigmAsINF+dnssup*SigmAsSUP+Ncc
        M_P2(k) = -dnsinf*SigmAsINF*(d-0.5*ht)+dnssup*SigmAsSUP*(0.5*ht-d0)+Mcc

    end do

!Traitement en PIVOT C - Entièrement Comprimée (EC) + Moment Positif
!-------------------------------------------------------------------

    allocate (N_P3(N_EC))
    allocate (M_P3(N_EC))

    do k = 1, N_EC

        X = (N_EC-k)/100.0
        if (k .eq. 1) then
            X = Xsup
        end if

        DE = X*Ec2
        EcINF = Ec2-DE*(1-xCt)
        EcSUP = DE+EcINF
        EsINF = EcINF+(DE/ht)*(ht-d)
        EsSUP = EcINF+(DE/ht)*(ht-d0)

        Ncc = bw*ht*fcd*(1+m2*(X**(nC)))
        Mcc = bw*ht*ht*fcd*m1*(X**(nC))
        Calc = EcSUP-EcINF
        if (abs(Calc) .gt. epsilon(Calc)) then
            alpha = (1/(1-EcINF/EcSUP))*(ht/d)
        else
            alpha = -1000
        end if

        if (typdiag .eq. 1) then
        if (Abs(EsINF) .lt. Ese) then
            SigmAsINF = eys*(Abs(EsINF))
        else
            SigmAsINF = fyd+((ktys*fyd-fyd)/(Euk-Ese))*(Abs(EsINF)-Ese)
        end if
        if (Abs(EsSUP) .lt. Ese) Then
            SigmAsSUP = eys*(Abs(EsSUP))
        else
            SigmAsSUP = fyd+((ktys*fyd-fyd)/(Euk-Ese))*(Abs(EsSUP)-Ese)
        end if
        else
        if (Abs(EsINF) .lt. Ese) then
            SigmAsINF = eys*(Abs(EsINF))
        else
            SigmAsINF = fyd
        end if
        if (Abs(EsSUP) .lt. Ese) then
            SigmAsSUP = eys*(Abs(EsSUP))
        else
            SigmAsSUP = fyd
        end if
        end if

        if (EsINF .lt. 0) then
            SigmAsINF = -SigmAsINF
        end if
        if (EsSUP .lt. 0) then
            SigmAsSUP = -SigmAsSUP
        end if

        N_P3(k) = dnsinf*SigmAsINF+dnssup*SigmAsSUP+Ncc
        M_P3(k) = -dnsinf*SigmAsINF*(d-0.5*ht)+dnssup*SigmAsSUP*(0.5*ht-d0)+Mcc

    end do

!Traitement en PIVOT C - Entièrement Comprimée (EC) + Moment Negatif
!-------------------------------------------------------------------
    allocate (N_P4(N_EC))
    allocate (M_P4(N_EC))

    do k = 1, N_EC
        X = (1-k)/100.0
        if (k .eq. N_EC) then
            X = -Xsup
        end if

        DE = X*Ec2
        EcINF = Ec2-DE*(1-xCt)
        EcSUP = DE+EcINF
        EsINF = EcINF+(DE/ht)*(ht-d)
        EsSUP = EcINF+(DE/ht)*(ht-d0)

        Ncc = bw*ht*fcd*(1+m2*(abs(X)**(nC)))
        Mcc = bw*ht*ht*fcd*m1*(abs(X)**(nC))
        Calc = EcSUP-EcINF
        if (abs(Calc) .gt. epsilon(Calc)) then
            alpha = (1/(1-EcINF/EcSUP))*(ht/d)
        else
            alpha = -1000
        end if

        if (typdiag .eq. 1) then
        if (Abs(EsINF) .lt. Ese) then
            SigmAsINF = eys*(Abs(EsINF))
        else
            SigmAsINF = fyd+((ktys*fyd-fyd)/(Euk-Ese))*(Abs(EsINF)-Ese)
        end if
        if (Abs(EsSUP) .lt. Ese) Then
            SigmAsSUP = eys*(Abs(EsSUP))
        else
            SigmAsSUP = fyd+((ktys*fyd-fyd)/(Euk-Ese))*(Abs(EsSUP)-Ese)
        end if
        else
        if (Abs(EsINF) .lt. Ese) then
            SigmAsINF = eys*(Abs(EsINF))
        else
            SigmAsINF = fyd
        end if
        if (Abs(EsSUP) .lt. Ese) then
            SigmAsSUP = eys*(Abs(EsSUP))
        else
            SigmAsSUP = fyd
        end if
        end if

        if (EsINF .lt. 0) then
            SigmAsINF = -SigmAsINF
        end if
        if (EsSUP .lt. 0) then
            SigmAsSUP = -SigmAsSUP
        end if

        N_P4(k) = dnsinf*SigmAsINF+dnssup*SigmAsSUP+Ncc
        M_P4(k) = -dnsinf*SigmAsINF*(d-0.5*ht)+dnssup*SigmAsSUP*(0.5*ht-d0)-Mcc

    end do

!Traitement en PIVOT A et B - Partiellement Comprimée (PC) + Moment Negatif
!--------------------------------------------------------------------------

    allocate (N_P5(N_PCN))
    allocate (M_P5(N_PCN))

    do k = 1, N_PCN

        if (k .gt. 1) then
            alpha = (N_PCN-k)*0.01
        else
            alpha = ht/dneg
        end if

        if (alpha .lt. alphaAB) then
            EsSUP = -Esu
            EcINF = Esu*alpha/(1-alpha)
        else
            EcINF = Ecu
            EsSUP = -Ecu*(1-alpha)/alpha
        end if

        EsINF = ((EsSUP-EcINF)/dneg)*(d0neg)+EcINF
        EcSUP = ((EsSUP-EcINF)/dneg)*(ht)+EcINF

        if (typdiag .eq. 1) then
        if (Abs(EsINF) .lt. Ese) then
            SigmAsINF = eys*(Abs(EsINF))
        else
            SigmAsINF = fyd+((ktys*fyd-fyd)/(Euk-Ese))*(Abs(EsINF)-Ese)
        end if
        if (Abs(EsSUP) .lt. Ese) Then
            SigmAsSUP = eys*(Abs(EsSUP))
        else
            SigmAsSUP = fyd+((ktys*fyd-fyd)/(Euk-Ese))*(Abs(EsSUP)-Ese)
        end if
        else
        if (Abs(EsINF) .lt. Ese) then
            SigmAsINF = eys*(Abs(EsINF))
        else
            SigmAsINF = fyd
        end if
        if (Abs(EsSUP) .lt. Ese) then
            SigmAsSUP = eys*(Abs(EsSUP))
        else
            SigmAsSUP = fyd
        end if
        end if

        if (EsINF .lt. 0) then
            SigmAsINF = -SigmAsINF
        end if
        if (EsSUP .lt. 0) then
            SigmAsSUP = -SigmAsSUP
        end if

        x1 = (EsSUP-EcINF)/Ec2
        y1 = EcINF/Ec2
        DELTA = dneg/ht

        if (EcINF .le. Ec2) then
            Beta = 0
        else
            yE = ((Ec2-EcINF)/(EsSUP-EcINF))*dneg
            Beta = yE/dneg
        end if

        COEF1 = (1-y1-alpha*x1)
        COEF2 = (1-y1-Beta*x1)
        if (abs(COEF1) .gt. epsilon(COEF1)) then
            VAR_COEF1 = Abs(COEF1)/COEF1
        else
            VAR_COEF1 = 1
        end if
        if (abs(COEF2) .gt. epsilon(COEF2)) then
            VAR_COEF2 = Abs(COEF2)/COEF2
        else
            VAR_COEF2 = 1
        end if

        Ncc = (fcd*bw*dneg)*(alpha+(1/((nC+1)*x1))*(VAR_COEF1*((Abs(COEF1))**(nC+1)) &
               & -VAR_COEF2*((Abs(COEF2))**(nC+1))))
        Mcc = (bw*dneg*dneg*fcd)*(0.5*alpha*(1/DELTA-alpha) &
               & +(1/(2*DELTA))*(1/((nC+1)*x1))*(VAR_COEF1*((Abs(COEF1))**(nC+1)) &
               & -VAR_COEF2*((Abs(COEF2))**(nC+1))) &
               & -(1/((nC+1)*x1))*(alpha*VAR_COEF1*((Abs(COEF1))**(nC+1)) &
               & -Beta*VAR_COEF2*((Abs(COEF2))**(nC+1))) &
               & -(1/((nC+1)*(nC+2)*x1*x1))*(VAR_COEF1*((Abs(COEF1))**(nC+2)) &
               & -VAR_COEF2*((Abs(COEF2))**(nC+2))))

        N_P5(k) = dnsinf*SigmAsINF+dnssup*SigmAsSUP+Ncc
        M_P5(k) = -dnsinf*SigmAsINF*(d-0.5*ht)+dnssup*SigmAsSUP*(0.5*ht-d0)-Mcc

    end do

!Traitement en PIVOT A - Entièrement Tendu (ET) + Moment Negatif
!---------------------------------------------------------------
    allocate (N_P6(N_ET))
    allocate (M_P6(N_ET))

    do k = 1, N_ET

        if (k .eq. N_ET) then
            EcINF = -Esu
        else
            EcINF = -(1.e-3)*(k-1)
        end if

        EsSUP = -Esu
        EsINF = ((EsSUP-EcINF)/dneg)*(d0neg)+EcINF
        EcSUP = ((EsSUP-EcINF)/dneg)*(ht)+EcINF

        if (typdiag .eq. 1) then
        if (Abs(EsINF) .lt. Ese) then
            SigmAsINF = eys*(Abs(EsINF))
        else
            SigmAsINF = fyd+((ktys*fyd-fyd)/(Euk-Ese))*(Abs(EsINF)-Ese)
        end if
        if (Abs(EsSUP) .lt. Ese) Then
            SigmAsSUP = eys*(Abs(EsSUP))
        else
            SigmAsSUP = fyd+((ktys*fyd-fyd)/(Euk-Ese))*(Abs(EsSUP)-Ese)
        end if
        else
        if (Abs(EsINF) .lt. Ese) then
            SigmAsINF = eys*(Abs(EsINF))
        else
            SigmAsINF = fyd
        end if
        if (Abs(EsSUP) .lt. Ese) then
            SigmAsSUP = eys*(Abs(EsSUP))
        else
            SigmAsSUP = fyd
        end if
        end if

        if (EsINF .lt. 0) then
            SigmAsINF = -SigmAsINF
        end if
        if (EsSUP .lt. 0) then
            SigmAsSUP = -SigmAsSUP
        end if

        N_P6(k) = dnsinf*SigmAsINF+dnssup*SigmAsSUP
        M_P6(k) = -dnsinf*SigmAsINF*(d-0.5*ht)+dnssup*SigmAsSUP*(0.5*ht-d0)

    end do

!-----------------------------------------------------------------------
!Fin de Traitement des différents cas
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    if (N_ET+N_ET+N_PC+N_EC+N_EC+N_PCN > ntot) then
        write (6, *) "IndexError: ntot argument must be greater than", &
            N_ET+N_ET+N_PC+N_EC+N_EC+N_PCN
    else

        do k = 1, N_ET
            nrd(k) = N_P1(k)
            mrd(k) = M_P1(k)
        end do
        do k = 1, N_PC
            nrd(k+N_ET) = N_P2(k)
            mrd(k+N_ET) = M_P2(k)
        end do
        do k = 1, N_EC
            nrd(k+N_ET+N_PC) = N_P3(k)
            mrd(k+N_ET+N_PC) = M_P3(k)
        end do
        do k = 1, N_EC
            nrd(k+N_ET+N_PC+N_EC) = N_P4(k)
            mrd(k+N_ET+N_PC+N_EC) = M_P4(k)
        end do
        do k = 1, N_PCN
            nrd(k+N_ET+N_PC+N_EC+N_EC) = N_P5(k)
            mrd(k+N_ET+N_PC+N_EC+N_EC) = M_P5(k)
        end do
        do k = 1, N_ET
            nrd(k+N_ET+N_PC+N_EC+N_EC+N_PCN) = N_P6(k)
            mrd(k+N_ET+N_PC+N_EC+N_EC+N_PCN) = M_P6(k)
        end do
    end if

    deallocate (N_P1)
    deallocate (N_P2)
    deallocate (N_P3)
    deallocate (N_P4)
    deallocate (N_P5)
    deallocate (N_P6)
    deallocate (M_P1)
    deallocate (M_P2)
    deallocate (M_P3)
    deallocate (M_P4)
    deallocate (M_P5)
    deallocate (M_P6)

end subroutine
