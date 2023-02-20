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

subroutine cafeluiter(typco, alphacc, effm, effn, ht, bw, &
                      enrobi, enrobs, facier, fbeton, gammas, gammac, &
                      clacier, eys, typdiag, ferrsyme, slsyme, uc, &
                      condns, astend, ascomp, sstend, sscomp, ectend, eccomp, &
                      alpha, pivot, etat, ierr)

!_______________________________________________________________________________________________
!
!      CAFELUITER
!
!      CALCUL DES ACIERS EN FLEXION COMPOSEE A L'ELU PAR ITERATION
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
!      I FERRSYME  FERRAILLAGE SYMETRIQUE?
!                     FERRSYME = 0 (NON)
!                     FERRSYME = 1 (OUI)
!      I SLSYME    SECTION SEUIL DE TOLERANCE POUR UN FERRAILLAGE SYMETRIQUE
!      I UC        UNITE DES CONTRAINTES :
!                     UC = 0 CONTRAINTES EN Pa
!                     UC = 1 CONTRAINTES EN MPa
!
!      O ASTEND       DENSITE DE L'ACIER TENDU
!      O ASCOMP       DENSITE DE L'ACIER COMPRIMÉ
!      O SSTEND       CONTRAINTE AU NIVEAU DE L'ACIER TENDU
!      O SSCOMP       CONTRAINTE AU NIVEAU DE L'ACIER COMPRIMÉ
!      O SCTEND       CONTRAINTE AU NIVEAU DE LA FIBRE DE BETON TENDU
!      O SCCOMP       CONTRAINTE AU NIVEAU DE LA FIBRE DE BETON COMPRIMÉ
!      O ALPHA        COEFFICIENT DE PROFONDEUR DE L'AN
!      O pivot        pivot DE FONCTIONNEMENT DE LA SECTION
!      O etat         etat DE FONCTIONNEMENT DE LA SECTION
!      O IERR         CODE RETOUR (0 = OK)
!_______________________________________________________________________________________________
!
    implicit none
#include "asterfort/verifelu.h"
#include "asterfort/wkvect.h"
#include "asterfort/jedetr.h"
#include "asterfort/juveca.h"
#include "asterfort/jeveuo.h"
!
!-----------------------------------------------------------------------
!!!!TERMES PRINCIPAUX D'ENTREE
!-----------------------------------------------------------------------

    integer :: typco
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
    integer :: clacier
    real(kind=8) :: eys
    integer :: typdiag
    integer :: ferrsyme
    real(kind=8) :: slsyme
    integer :: uc
    logical :: condns
    real(kind=8) :: astend
    real(kind=8) :: ascomp
    real(kind=8) :: sstend
    real(kind=8) :: sscomp
    real(kind=8) :: ectend
    real(kind=8) :: eccomp
    real(kind=8) :: alpha
    integer :: pivot
    integer :: etat
    integer :: ierr

!-----------------------------------------------------------------------
!!!!VARIABLES DE CALCUL
!-----------------------------------------------------------------------
    real(kind=8) :: enrob, d, d0
    real(kind=8) :: unite_pa
    real(kind=8) :: fyd, fcd, nC, ktys, xC, yC, xCt, D00, m1, m2
    real(kind=8) :: Esu, Euk, Ecu, Ec2, Ese
    real(kind=8) :: Xsup, Ncc, Mcc
    real(kind=8) :: piv_a, piv_b, piv_c, alphaAB
    real(kind=8) :: COEF1, COEF2, VAR_COEF1, VAR_COEF2
    real(kind=8) :: Calc, a11, a12, a21, a22, f1, f2, x1, y1, DELTA, Beta, yE
    real(kind=8) :: AsTOT, escomp, estend
    integer :: N_ET, N_PC, N_EC, k, kFIN, q, qsy, s, N_TOT, i, verif
    real(kind=8) :: DE, X

    character(20) :: p(56)

    real(kind=8), pointer :: ectend_ET(:) => null(), eccomp_ET(:) => null()
    real(kind=8), pointer :: estend_ET(:) => null(), escomp_ET(:) => null()
    real(kind=8), pointer :: astend_ET(:) => null(), ascomp_ET(:) => null()
    real(kind=8), pointer :: AsTOT_ET(:) => null(), alpha_ET(:) => null()
    real(kind=8), pointer :: sstend_ET(:) => null(), sscomp_ET(:) => null()
    integer, pointer :: pivot_ET(:) => null(), etat_ET(:) => null()

    real(kind=8), pointer :: ectend_PC(:) => null(), eccomp_PC(:) => null()
    real(kind=8), pointer :: estend_PC(:) => null(), escomp_PC(:) => null()
    real(kind=8), pointer :: astend_PC(:) => null(), ascomp_PC(:) => null()
    real(kind=8), pointer :: AsTOT_PC(:) => null(), alpha_PC(:) => null()
    real(kind=8), pointer :: sstend_PC(:) => null(), sscomp_PC(:) => null()
    integer, pointer :: pivot_PC(:) => null(), etat_PC(:) => null()

    real(kind=8), pointer :: ectend_EC(:) => null(), eccomp_EC(:) => null()
    real(kind=8), pointer :: estend_EC(:) => null(), escomp_EC(:) => null()
    real(kind=8), pointer :: astend_EC(:) => null(), ascomp_EC(:) => null()
    real(kind=8), pointer :: AsTOT_EC(:) => null(), alpha_EC(:) => null()
    real(kind=8), pointer :: sstend_EC(:) => null(), sscomp_EC(:) => null()
    integer, pointer :: pivot_EC(:) => null(), etat_EC(:) => null()

    real(kind=8), pointer :: ectend_TOT(:) => null(), eccomp_TOT(:) => null()
    real(kind=8), pointer :: estend_TOT(:) => null(), escomp_TOT(:) => null()
    real(kind=8), pointer :: astend_TOT(:) => null(), ascomp_TOT(:) => null()
    real(kind=8), pointer :: AsTOT_TOT(:) => null(), alpha_TOT(:) => null()
    real(kind=8), pointer :: sstend_TOT(:) => null(), sscomp_TOT(:) => null()
    integer, pointer :: pivot_TOT(:) => null(), etat_TOT(:) => null()

    real(kind=8), pointer :: Pente_astend_ET(:) => null()
    real(kind=8), pointer :: Pente_ascomp_ET(:) => null()
    real(kind=8), pointer :: Pente_astend_PC(:) => null()
    real(kind=8), pointer :: Pente_ascomp_PC(:) => null()
    real(kind=8), pointer :: Pente_astend_EC(:) => null()
    real(kind=8), pointer :: Pente_ascomp_EC(:) => null()
    real(kind=8), pointer :: Pente_astend_TOT(:) => null()
    real(kind=8), pointer :: Pente_ascomp_TOT(:) => null()

!-----------------------------------------------------------------------
!!!!LANCEMENT DU CALCUL
!-----------------------------------------------------------------------

!ET = Entièrement Tendue
!PC = Partiellement Comprimée
!EC = Entièrement Comprimée

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

!   Paramètres de calcul
    Xsup = piv_b/piv_c
    xC = (1-piv_c/piv_b)*ht
    yC = ht-xC
    xCt = xC/ht
    D00 = (ht-d0)/ht
    m1 = (((1-xCt)**(nC+1))/(2.d0*(nC+1)))*(1-(2.d0*(1-xCt))/(nC+2))
    m2 = -((1-xCt)**(nC+1))/(nC+1)

    Esu = piv_a
    Euk = Esu/0.9
    Ecu = piv_b
    Ec2 = piv_c
    Ese = fyd/eys
    alphaAB = 1./(1+Esu/Ecu)

    do i = 1, 56
        write (p(i), fmt='(A18,I2)') 'POINT_ITER_CAFELU_', i
    end do

!Traitement en pivot A - Entièrement tendu (ET)
!----------------------------------------------

    N_ET = floor(Esu*1000)+1

    call wkvect(p(1), ' V V R ', N_ET, vr=ectend_ET)
    call wkvect(p(2), ' V V R ', N_ET, vr=eccomp_ET)
    call wkvect(p(3), ' V V R ', N_ET, vr=estend_ET)
    call wkvect(p(4), ' V V R ', N_ET, vr=escomp_ET)
    call wkvect(p(5), ' V V R ', N_ET, vr=astend_ET)
    call wkvect(p(6), ' V V R ', N_ET, vr=ascomp_ET)
    call wkvect(p(7), ' V V R ', N_ET, vr=Pente_astend_ET)
    call wkvect(p(8), ' V V R ', N_ET, vr=Pente_ascomp_ET)
    call wkvect(p(9), ' V V R ', N_ET, vr=AsTOT_ET)
    call wkvect(p(10), ' V V R ', N_ET, vr=sstend_ET)
    call wkvect(p(11), ' V V R ', N_ET, vr=sscomp_ET)
    call wkvect(p(12), ' V V R ', N_ET, vr=alpha_ET)
    call wkvect(p(13), ' V V I ', N_ET, vi=pivot_ET)
    call wkvect(p(14), ' V V I ', N_ET, vi=etat_ET)

    do k = 1, N_ET

        if (k .eq. 1) then
            eccomp_ET(k) = -Esu
        else
            eccomp_ET(k) = -(1.e-3)*(N_ET-k)
        end if

        eccomp = eccomp_ET(k)
        estend = -Esu
        escomp = ((estend-eccomp)/d)*(d0)+eccomp
        ectend = ((estend-eccomp)/d)*(ht)+eccomp

        ectend_ET(k) = ectend
        escomp_ET(k) = escomp
        estend_ET(k) = estend
        pivot_ET(k) = 1
        etat_ET(k) = 4
        Calc = estend-eccomp
        if (abs(eccomp) .lt. epsilon(eccomp)) then
            alpha = 0
        elseif (abs(Calc) .gt. epsilon(Calc)) then
            alpha = -(eccomp/(estend-eccomp))
        else
            alpha = -1000
        end if

        alpha_ET(k) = alpha

        if (Abs(estend) .lt. Ese) then
            sstend = eys*(Abs(estend))
        elseif (typdiag .eq. 1) then
            sstend = fyd+((ktys*fyd-fyd)/(Euk-Ese))*(Abs(estend)-Ese)
        else
            sstend = fyd
        end if
        if (Abs(escomp) .lt. Ese) then
            sscomp = eys*(Abs(escomp))
        elseif (typdiag .eq. 1) then
            sscomp = fyd+((ktys*fyd-fyd)/(Euk-Ese))*(Abs(escomp)-Ese)
        else
            sscomp = fyd
        end if
        if (estend .lt. 0) then
            sstend = -sstend
        end if
        if (escomp .lt. 0) then
            sscomp = -sscomp
        end if

        sstend_ET(k) = sstend
        sscomp_ET(k) = sscomp

        astend_ET(k) = (abs(effm)-(0.5*ht-d0)*effn)/(Abs(sstend)*(d-d0))
        ascomp_ET(k) = -(effn+astend_ET(k)*Abs(sstend))/Abs(sscomp)

        if (k .gt. 1) then
            Pente_astend_ET(k) = (astend_ET(k)-astend_ET(k-1))/(alpha_ET(k)-alpha_ET(k-1))
            Pente_ascomp_ET(k) = (ascomp_ET(k)-ascomp_ET(k-1))/(alpha_ET(k)-alpha_ET(k-1))
        else
            Pente_astend_ET(k) = 0
            Pente_ascomp_ET(k) = 0
        end if

        AsTOT_ET(k) = astend_ET(k)+ascomp_ET(k)

    end do

!Traitement en pivot A et B - Partiellement Comprimée (PC)
!---------------------------------------------------------

    N_PC = ceiling((ht/d)*1000)+1

    call wkvect(p(15), ' V V R ', N_PC, vr=ectend_PC)
    call wkvect(p(16), ' V V R ', N_PC, vr=eccomp_PC)
    call wkvect(p(17), ' V V R ', N_PC, vr=estend_PC)
    call wkvect(p(18), ' V V R ', N_PC, vr=escomp_PC)
    call wkvect(p(19), ' V V R ', N_PC, vr=astend_PC)
    call wkvect(p(20), ' V V R ', N_PC, vr=ascomp_PC)
    call wkvect(p(21), ' V V R ', N_PC, vr=Pente_astend_PC)
    call wkvect(p(22), ' V V R ', N_PC, vr=Pente_ascomp_PC)
    call wkvect(p(23), ' V V R ', N_PC, vr=AsTOT_PC)
    call wkvect(p(24), ' V V R ', N_PC, vr=sstend_PC)
    call wkvect(p(25), ' V V R ', N_PC, vr=sscomp_PC)
    call wkvect(p(26), ' V V R ', N_PC, vr=alpha_PC)
    call wkvect(p(27), ' V V I ', N_PC, vi=pivot_PC)
    call wkvect(p(28), ' V V I ', N_PC, vi=etat_PC)

    do k = 1, N_PC

        alpha_PC(k) = (k-1)*0.001
        alpha = alpha_PC(k)

        if (alpha .lt. alphaAB) then
            pivot = 1
            estend = -Esu
            eccomp = Esu*alpha/(1-alpha)
        else
            pivot = 2
            eccomp = Ecu
            estend = -Ecu*(1-alpha)/alpha
        end if

        escomp = ((estend-eccomp)/d)*(d0)+eccomp
        ectend = ((estend-eccomp)/d)*(ht)+eccomp

        ectend_PC(k) = ectend
        escomp_PC(k) = escomp
        estend_PC(k) = estend
        eccomp_PC(k) = eccomp
        pivot_PC(k) = pivot
        etat_PC(k) = 5

        if (Abs(estend) .lt. Ese) then
            sstend = eys*(Abs(estend))
        elseif (typdiag .eq. 1) then
            sstend = fyd+((ktys*fyd-fyd)/(Euk-Ese))*(Abs(estend)-Ese)
        else
            sstend = fyd
        end if
        if (Abs(escomp) .lt. Ese) then
            sscomp = eys*(Abs(escomp))
        elseif (typdiag .eq. 1) then
            sscomp = fyd+((ktys*fyd-fyd)/(Euk-Ese))*(Abs(escomp)-Ese)
        else
            sscomp = fyd
        end if
        if (estend .lt. 0) then
            sstend = -sstend
        end if
        if (escomp .lt. 0) then
            sscomp = -sscomp
        end if

        sstend_PC(k) = sstend
        sscomp_PC(k) = sscomp

        x1 = (estend-eccomp)/Ec2
        y1 = eccomp/Ec2
        DELTA = d/ht

        if (eccomp .le. Ec2) then
            Beta = 0
        else
            yE = ((Ec2-eccomp)/(estend-eccomp))*d
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

        a11 = sscomp
        a12 = sstend
        a21 = sscomp*(0.5*ht-d0)
        a22 = -sstend*(d-0.5*ht)
        f1 = effn-Ncc
        f2 = abs(effm)-Mcc

        if (k .eq. 1001) then
            ascomp = (effn-Ncc)/sscomp
            astend = 0
            Calc = abs(effm)-Mcc-ascomp*sscomp*(0.5*ht-d0)
            if (Abs(Calc/effm) .gt. 0.01) then
                ascomp = -1
                astend = -1
            end if
        else
            ascomp = (f1*a22-a12*f2)/(a11*a22-a12*a21)
            astend = (a11*f2-a21*f1)/(a11*a22-a12*a21)
        end if

        ascomp_PC(k) = ascomp
        astend_PC(k) = astend

        if (k .gt. 1) then
            Pente_astend_PC(k) = (astend_PC(k)-astend_PC(k-1))/(alpha_PC(k)-alpha_PC(k-1))
            Pente_ascomp_PC(k) = (ascomp_PC(k)-ascomp_PC(k-1))/(alpha_PC(k)-alpha_PC(k-1))
        end if

        AsTOT_PC(k) = ascomp+astend

    end do

!Traitement en pivot C - Entièrement Comprimée (EC)
!--------------------------------------------------

    N_EC = ceiling(Xsup*100)+1

    call wkvect(p(29), ' V V R ', N_EC, vr=ectend_EC)
    call wkvect(p(30), ' V V R ', N_EC, vr=eccomp_EC)
    call wkvect(p(31), ' V V R ', N_EC, vr=estend_EC)
    call wkvect(p(32), ' V V R ', N_EC, vr=escomp_EC)
    call wkvect(p(33), ' V V R ', N_EC, vr=astend_EC)
    call wkvect(p(34), ' V V R ', N_EC, vr=ascomp_EC)
    call wkvect(p(35), ' V V R ', N_EC, vr=Pente_astend_EC)
    call wkvect(p(36), ' V V R ', N_EC, vr=Pente_ascomp_EC)
    call wkvect(p(37), ' V V R ', N_EC, vr=AsTOT_EC)
    call wkvect(p(38), ' V V R ', N_EC, vr=sstend_EC)
    call wkvect(p(39), ' V V R ', N_EC, vr=sscomp_EC)
    call wkvect(p(40), ' V V R ', N_EC, vr=alpha_EC)
    call wkvect(p(41), ' V V I ', N_EC, vi=pivot_EC)
    call wkvect(p(42), ' V V I ', N_EC, vi=etat_EC)

    do k = 1, N_EC

        X = (N_EC-k)/100.0
        if (k .eq. 1) then
            X = Xsup
        end if

        DE = X*Ec2
        ectend = Ec2-DE*(1-xCt)
        eccomp = DE+ectend
        estend = ectend+(DE/ht)*(ht-d)
        escomp = ectend+(DE/ht)*(ht-d0)

        eccomp_EC(k) = eccomp
        ectend_EC(k) = ectend
        escomp_EC(k) = escomp
        estend_EC(k) = estend
        pivot_EC(k) = 3
        etat_EC(k) = 6

        Ncc = bw*ht*fcd*(1+m2*(X**(nC)))
        Mcc = bw*ht*ht*fcd*m1*(X**(nC))
        if (abs(eccomp) .gt. epsilon(eccomp)) then
            Calc = 1-ectend/eccomp
        else
            Calc = 0
        end if
        if (abs(Calc) .gt. epsilon(Calc)) then
            alpha = (1/(1-ectend/eccomp))*(ht/d)
        else
            alpha = -1000
        end if

        alpha_EC(k) = alpha

        if (Abs(estend) .lt. Ese) then
            sstend = eys*(Abs(estend))
        elseif (typdiag .eq. 1) then
            sstend = fyd+((ktys*fyd-fyd)/(Euk-Ese))*(Abs(estend)-Ese)
        else
            sstend = fyd
        end if
        if (Abs(escomp) .lt. Ese) then
            sscomp = eys*(Abs(escomp))
        elseif (typdiag .eq. 1) then
            sscomp = fyd+((ktys*fyd-fyd)/(Euk-Ese))*(Abs(escomp)-Ese)
        else
            sscomp = fyd
        end if
        if (estend .lt. 0) then
            sstend = -sstend
        end if
        if (escomp .lt. 0) then
            sscomp = -sscomp
        end if

        sstend_EC(k) = sstend
        sscomp_EC(k) = sscomp

        a11 = sscomp
        a12 = sstend
        a21 = sscomp*(0.5*ht-d0)
        a22 = -sstend*(d-0.5*ht)
        f1 = effn-Ncc
        f2 = abs(effm)-Mcc

        ascomp = (f1*a22-a12*f2)/(a11*a22-a12*a21)
        astend = (a11*f2-a21*f1)/(a11*a22-a12*a21)

        ascomp_EC(k) = ascomp
        astend_EC(k) = astend

        if (k .gt. 1) then
            Pente_astend_EC(k) = (astend_EC(k)-astend_EC(k-1))/(alpha_EC(k)-alpha_EC(k-1))
            Pente_ascomp_EC(k) = (ascomp_EC(k)-ascomp_EC(k-1))/(alpha_EC(k)-alpha_EC(k-1))
        end if

        AsTOT_EC(k) = ascomp+astend

    end do

!-----------------------------------------------------------------------
!Fin de Traitement des différents cas
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

    N_TOT = N_ET+N_PC+N_EC

    call wkvect(p(43), ' V V R ', N_TOT, vr=ectend_TOT)
    call wkvect(p(44), ' V V R ', N_TOT, vr=eccomp_TOT)
    call wkvect(p(45), ' V V R ', N_TOT, vr=estend_TOT)
    call wkvect(p(46), ' V V R ', N_TOT, vr=escomp_TOT)
    call wkvect(p(47), ' V V R ', N_TOT, vr=astend_TOT)
    call wkvect(p(48), ' V V R ', N_TOT, vr=ascomp_TOT)
    call wkvect(p(49), ' V V R ', N_TOT, vr=Pente_astend_TOT)
    call wkvect(p(50), ' V V R ', N_TOT, vr=Pente_ascomp_TOT)
    call wkvect(p(51), ' V V R ', N_TOT, vr=AsTOT_TOT)
    call wkvect(p(52), ' V V R ', N_TOT, vr=sstend_TOT)
    call wkvect(p(53), ' V V R ', N_TOT, vr=sscomp_TOT)
    call wkvect(p(54), ' V V R ', N_TOT, vr=alpha_TOT)
    call wkvect(p(55), ' V V I ', N_TOT, vi=pivot_TOT)
    call wkvect(p(56), ' V V I ', N_TOT, vi=etat_TOT)

    do k = 1, N_ET
        ectend_TOT(k) = ectend_ET(k)
        eccomp_TOT(k) = eccomp_ET(k)
        estend_TOT(k) = estend_ET(k)
        escomp_TOT(k) = escomp_ET(k)
        astend_TOT(k) = astend_ET(k)
        ascomp_TOT(k) = ascomp_ET(k)
        Pente_astend_TOT(k) = Pente_astend_ET(k)
        Pente_ascomp_TOT(k) = Pente_ascomp_ET(k)
        AsTOT_TOT(k) = AsTOT_ET(k)
        sstend_TOT(k) = sstend_ET(k)
        sscomp_TOT(k) = sscomp_ET(k)
        alpha_TOT(k) = alpha_ET(k)
        pivot_TOT(k) = pivot_ET(k)
        etat_TOT(k) = etat_ET(k)
    end do

    do k = 1, N_PC
        ectend_TOT(N_ET+k) = ectend_PC(k)
        eccomp_TOT(N_ET+k) = eccomp_PC(k)
        estend_TOT(N_ET+k) = estend_PC(k)
        escomp_TOT(N_ET+k) = escomp_PC(k)
        astend_TOT(N_ET+k) = astend_PC(k)
        ascomp_TOT(N_ET+k) = ascomp_PC(k)
        Pente_astend_TOT(N_ET+k) = Pente_astend_PC(k)
        Pente_ascomp_TOT(N_ET+k) = Pente_ascomp_PC(k)
        AsTOT_TOT(N_ET+k) = AsTOT_PC(k)
        sstend_TOT(N_ET+k) = sstend_PC(k)
        sscomp_TOT(N_ET+k) = sscomp_PC(k)
        alpha_TOT(N_ET+k) = alpha_PC(k)
        pivot_TOT(N_ET+k) = pivot_PC(k)
        etat_TOT(N_ET+k) = etat_PC(k)
    end do

    do k = 1, N_EC
        ectend_TOT(N_ET+N_PC+k) = ectend_EC(k)
        eccomp_TOT(N_ET+N_PC+k) = eccomp_EC(k)
        estend_TOT(N_ET+N_PC+k) = estend_EC(k)
        escomp_TOT(N_ET+N_PC+k) = escomp_EC(k)
        astend_TOT(N_ET+N_PC+k) = astend_EC(k)
        ascomp_TOT(N_ET+N_PC+k) = ascomp_EC(k)
        Pente_astend_TOT(N_ET+N_PC+k) = Pente_astend_EC(k)
        Pente_ascomp_TOT(N_ET+N_PC+k) = Pente_ascomp_EC(k)
        AsTOT_TOT(N_ET+N_PC+k) = AsTOT_EC(k)
        sstend_TOT(N_ET+N_PC+k) = sstend_EC(k)
        sscomp_TOT(N_ET+N_PC+k) = sscomp_EC(k)
        alpha_TOT(N_ET+N_PC+k) = alpha_EC(k)
        pivot_TOT(N_ET+N_PC+k) = pivot_EC(k)
        etat_TOT(N_ET+N_PC+k) = etat_EC(k)
    end do

   !!Recherche des racines éventuelles

    k = 1
    do while (k .le. (N_TOT-1))
        if (((astend_TOT(k)*astend_TOT(k+1)) .le. 0) &
              & .and. ((Pente_astend_TOT(k)*Pente_astend_TOT(k+1)) .gt. 0) &
              & .and. (ascomp_TOT(k) .ge. 0) .and. (ascomp_TOT(k+1) .ge. 0)) then

            N_TOT = N_TOT+1

            do i = 43, 56
                call juveca(p(i), N_TOT)
            end do

            call jeveuo(p(43), 'E', vr=ectend_TOT)
            call jeveuo(p(44), 'E', vr=eccomp_TOT)
            call jeveuo(p(45), 'E', vr=estend_TOT)
            call jeveuo(p(46), 'E', vr=escomp_TOT)
            call jeveuo(p(47), 'E', vr=astend_TOT)
            call jeveuo(p(48), 'E', vr=ascomp_TOT)
            call jeveuo(p(49), 'E', vr=Pente_astend_TOT)
            call jeveuo(p(50), 'E', vr=Pente_ascomp_TOT)
            call jeveuo(p(51), 'E', vr=AsTOT_TOT)
            call jeveuo(p(52), 'E', vr=sstend_TOT)
            call jeveuo(p(53), 'E', vr=sscomp_TOT)
            call jeveuo(p(54), 'E', vr=alpha_TOT)
            call jeveuo(p(55), 'E', vi=pivot_TOT)
            call jeveuo(p(56), 'E', vi=etat_TOT)

            do s = 1, (N_TOT-k-1)
                ectend_TOT(N_TOT-s+1) = ectend_TOT(N_TOT-s)
                eccomp_TOT(N_TOT-s+1) = eccomp_TOT(N_TOT-s)
                estend_TOT(N_TOT-s+1) = estend_TOT(N_TOT-s)
                escomp_TOT(N_TOT-s+1) = escomp_TOT(N_TOT-s)
                astend_TOT(N_TOT-s+1) = astend_TOT(N_TOT-s)
                ascomp_TOT(N_TOT-s+1) = ascomp_TOT(N_TOT-s)
                Pente_astend_TOT(N_TOT-s+1) = Pente_astend_TOT(N_TOT-s)
                Pente_ascomp_TOT(N_TOT-s+1) = Pente_ascomp_TOT(N_TOT-s)
                AsTOT_TOT(N_TOT-s+1) = AsTOT_TOT(N_TOT-s)
                sstend_TOT(N_TOT-s+1) = sstend_TOT(N_TOT-s)
                sscomp_TOT(N_TOT-s+1) = sscomp_TOT(N_TOT-s)
                alpha_TOT(N_TOT-s+1) = alpha_TOT(N_TOT-s)
                pivot_TOT(N_TOT-s+1) = pivot_TOT(N_TOT-s)
                etat_TOT(N_TOT-s+1) = etat_TOT(N_TOT-s)
            end do

            ectend_TOT(k+1) = 0.5*(ectend_TOT(k)+ectend_TOT(k+2))
            eccomp_TOT(k+1) = 0.5*(eccomp_TOT(k)+eccomp_TOT(k+2))
            estend_TOT(k+1) = 0.5*(estend_TOT(k)+estend_TOT(k+2))
            escomp_TOT(k+1) = 0.5*(escomp_TOT(k)+escomp_TOT(k+2))
            astend_TOT(k+1) = 0
            ascomp_TOT(k+1) = 0.5*(ascomp_TOT(k)+ascomp_TOT(k+2))
            AsTOT_TOT(k+1) = ascomp_TOT(k+1)+astend_TOT(k+1)
            sstend_TOT(k+1) = 0.5*(sstend_TOT(k)+sstend_TOT(k+2))
            sscomp_TOT(k+1) = 0.5*(sscomp_TOT(k)+sscomp_TOT(k+2))
            alpha_TOT(k+1) = 0.5*(alpha_TOT(k)+alpha_TOT(k+2))
            pivot_TOT(k+1) = pivot_TOT(k)
            etat_TOT(k+1) = etat_TOT(k)
            Pente_astend_TOT(k+1) = (astend_TOT(k+1)-astend_TOT(k))/(alpha_TOT(k+1)-alpha_TOT(k))
            Pente_ascomp_TOT(k+1) = (ascomp_TOT(k+1)-ascomp_TOT(k))/(alpha_TOT(k+1)-alpha_TOT(k))

            k = k+2

        elseif (((ascomp_TOT(k)*ascomp_TOT(k+1)) .le. 0) &
                  & .and. ((Pente_ascomp_TOT(k)*Pente_ascomp_TOT(k+1)) .gt. 0) &
                  & .and. (astend_TOT(k) .ge. 0) .and. (astend_TOT(k+1) .ge. 0)) then

            N_TOT = N_TOT+1

            do i = 43, 56
                call juveca(p(i), N_TOT)
            end do

            call jeveuo(p(43), 'E', vr=ectend_TOT)
            call jeveuo(p(44), 'E', vr=eccomp_TOT)
            call jeveuo(p(45), 'E', vr=estend_TOT)
            call jeveuo(p(46), 'E', vr=escomp_TOT)
            call jeveuo(p(47), 'E', vr=astend_TOT)
            call jeveuo(p(48), 'E', vr=ascomp_TOT)
            call jeveuo(p(49), 'E', vr=Pente_astend_TOT)
            call jeveuo(p(50), 'E', vr=Pente_ascomp_TOT)
            call jeveuo(p(51), 'E', vr=AsTOT_TOT)
            call jeveuo(p(52), 'E', vr=sstend_TOT)
            call jeveuo(p(53), 'E', vr=sscomp_TOT)
            call jeveuo(p(54), 'E', vr=alpha_TOT)
            call jeveuo(p(55), 'E', vi=pivot_TOT)
            call jeveuo(p(56), 'E', vi=etat_TOT)

            do s = 1, (N_TOT-k-1)
                ectend_TOT(N_TOT-s+1) = ectend_TOT(N_TOT-s)
                eccomp_TOT(N_TOT-s+1) = eccomp_TOT(N_TOT-s)
                estend_TOT(N_TOT-s+1) = estend_TOT(N_TOT-s)
                escomp_TOT(N_TOT-s+1) = escomp_TOT(N_TOT-s)
                astend_TOT(N_TOT-s+1) = astend_TOT(N_TOT-s)
                ascomp_TOT(N_TOT-s+1) = ascomp_TOT(N_TOT-s)
                Pente_astend_TOT(N_TOT-s+1) = Pente_astend_TOT(N_TOT-s)
                Pente_ascomp_TOT(N_TOT-s+1) = Pente_ascomp_TOT(N_TOT-s)
                AsTOT_TOT(N_TOT-s+1) = AsTOT_TOT(N_TOT-s)
                sstend_TOT(N_TOT-s+1) = sstend_TOT(N_TOT-s)
                sscomp_TOT(N_TOT-s+1) = sscomp_TOT(N_TOT-s)
                alpha_TOT(N_TOT-s+1) = alpha_TOT(N_TOT-s)
                pivot_TOT(N_TOT-s+1) = pivot_TOT(N_TOT-s)
                etat_TOT(N_TOT-s+1) = etat_TOT(N_TOT-s)
            end do

            ectend_TOT(k+1) = 0.5*(ectend_TOT(k)+ectend_TOT(k+2))
            eccomp_TOT(k+1) = 0.5*(eccomp_TOT(k)+eccomp_TOT(k+2))
            estend_TOT(k+1) = 0.5*(estend_TOT(k)+estend_TOT(k+2))
            escomp_TOT(k+1) = 0.5*(escomp_TOT(k)+escomp_TOT(k+2))
            astend_TOT(k+1) = 0.5*(astend_TOT(k)+astend_TOT(k+2))
            ascomp_TOT(k+1) = 0
            AsTOT_TOT(k+1) = ascomp_TOT(k+1)+astend_TOT(k+1)
            sstend_TOT(k+1) = 0.5*(sstend_TOT(k)+sstend_TOT(k+2))
            sscomp_TOT(k+1) = 0.5*(sscomp_TOT(k)+sscomp_TOT(k+2))
            alpha_TOT(k+1) = 0.5*(alpha_TOT(k)+alpha_TOT(k+2))
            pivot_TOT(k+1) = pivot_TOT(k)
            etat_TOT(k+1) = etat_TOT(k)
            Pente_astend_TOT(k+1) = (astend_TOT(k+1)-astend_TOT(k))/(alpha_TOT(k+1)-alpha_TOT(k))
            Pente_ascomp_TOT(k+1) = (ascomp_TOT(k+1)-ascomp_TOT(k))/(alpha_TOT(k+1)-alpha_TOT(k))

            k = k+2

        else

            k = k+1

        end if

    end do

!! DETERMINATION DE LA SOLUTION DIMENSIONNANTE

    q = 0
    qsy = 0
    kFIN = 1

    do k = 1, N_TOT
        if ((ascomp_TOT(k) .ge. 0) .and. (astend_TOT(k) .ge. 0)) then
            q = q+1
            if (ferrsyme .eq. 1) then
                Calc = abs(ascomp_TOT(k)-astend_TOT(k))
                if (Calc .le. slsyme) then
                    qsy = qsy+1
                    if (qsy .eq. 1) then
                        AsTOT = AsTOT_TOT(k)
                        kFIN = k
                    else
                        if (AsTOT_TOT(k) .lt. AsTOT) Then
                            AsTOT = AsTOT_TOT(k)
                            kFIN = k
                        end if
                    end if
                end if
            else
                if (q .eq. 1) then
                    AsTOT = AsTOT_TOT(k)
                    kFIN = k
                else
                    if (AsTOT_TOT(k) .lt. AsTOT) Then
                        AsTOT = AsTOT_TOT(k)
                        kFIN = k
                    end if
                end if
            end if
        end if
    end do

    if (((q .gt. 0) .and. (ferrsyme .eq. 0)) &
        & .or. ((qsy .gt. 0) .and. (ferrsyme .eq. 1))) then

        alpha = alpha_TOT(kFIN)
        ascomp = ascomp_TOT(kFIN)
        astend = astend_TOT(kFIN)
        eccomp = eccomp_TOT(kFIN)
        ectend = ectend_TOT(kFIN)
        estend = estend_TOT(kFIN)
        escomp = escomp_TOT(kFIN)
        pivot = pivot_TOT(kFIN)
        etat = etat_TOT(kFIN)
        sstend = sstend_TOT(kFIN)
        sscomp = sscomp_TOT(kFIN)

    elseif ((ferrsyme .eq. 1) .and. ((q .gt. 0) .or. (condns .eqv. (.true.)))) then

        etat = 0
        pivot = 0
        alpha = -1000
        eccomp = 0
        ectend = -1
        estend = -1
        escomp = -1
        sscomp = -1
        sstend = -1
        ascomp = -1
        astend = -1
        ierr = 2

    else

    !! Vérification par le diagramme d'interaction avec dnsyi=dnsys=0.0
        call verifelu(typco, alphacc, ht, bw, enrobi, enrobs, facier, fbeton, &
                      gammas, gammac, clacier, eys, typdiag, uc, &
                      0.0, 0.0, effm, effn, verif)

        if (verif .eq. 0) then
           !!OK
            etat = 1
            pivot = 0
            alpha = -1000
            eccomp = -1
            ectend = -1
            estend = -1
            escomp = -1
            sscomp = -1
            sstend = -1
            ascomp = 0
            astend = 0
        else
            etat = 0
            pivot = 0
            alpha = -1000
            eccomp = -1
            ectend = -1
            estend = -1
            escomp = -1
            sscomp = -1
            sstend = -1
            ascomp = -1
            astend = -1
            ierr = 1
        end if

    end if

    do i = 1, 56
        call jedetr(p(i))
    end do

end subroutine
