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

subroutine dintelu(typco, alphacc, ht, bw, enrobi, enrobs, facier, fbeton,&
                   gammas, gammac, clacier, eys, typdiag, uc,&
                   dnsinf, dnssup, ntot, nrd, mrd)
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
!      O NTOT      DIMENSIONS DES VECTEURS
!      O NRD       VECTEUR DES EFFORTS NORMAUX RESISTANTS (DIAG INTER)
!      O MRD       VECTEUR DES MOMENTS RESISTANTS (DIAG INTER)
!
!______________________________________________________________________
!
!
    implicit none
!
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/jedetr.h"
!
    integer :: typco
    real(kind=8) :: alphacc
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
    integer :: uc
    real(kind=8) :: dnsinf
    real(kind=8) :: dnssup
    integer :: ntot
    real(kind=8) :: nrd(1:ntot)
    real(kind=8) :: mrd(1:ntot)

!-----------------------------------------------------------------------
!!!!VARIABLES DE CALCUL
!-----------------------------------------------------------------------

    real(kind=8) :: d,d0,dneg,d0neg
    real(kind=8) :: unite_pa,DE,X,Calc
    real(kind=8) :: fyd,fcd,nC,ktys,xC,yC,xCt,m1,m2
    real(kind=8) :: Esu,Euk,Ecu,Ec2,Ese,Xsup
    real(kind=8) :: piv_a,piv_b,piv_c,alpha,alphaAB,alphaR,alphaBC
    real(kind=8) :: COEF1,COEF2,VAR_COEF1,VAR_COEF2
    real(kind=8) :: DELTA,x1,y1,Beta,yE
    real(kind=8) :: EcINF,EcSUP,EsSUP,EsINF
    real(kind=8) :: SigmAsSUP,SigmAsINF,Ncc,Mcc
    integer :: N_ET,N_PC,N_EC,N_PCN,k

    character(24) :: p01, p02
    character(24) :: p03, p04
    character(24) :: p05, p06
    character(24) :: p07, p08
    character(24) :: p09, p10
    character(24) :: p11, p12

    real(kind=8), pointer :: N_P1(:) => null(), M_P1(:) => null()
    real(kind=8), pointer :: N_P2(:) => null(), M_P2(:) => null()
    real(kind=8), pointer :: N_P3(:) => null(), M_P3(:) => null()
    real(kind=8), pointer :: N_P4(:) => null(), M_P4(:) => null()
    real(kind=8), pointer :: N_P5(:) => null(), M_P5(:) => null()
    real(kind=8), pointer :: N_P6(:) => null(), M_P6(:) => null()
    
!   Initialisation des entiers
    
    N_ET = 100
    N_PC = 100
    N_EC = 100
    N_PCN = 100
    k = 1

    if (typco.eq.1) then
!       CALCUL DES PARAMETRES POUR CODIFICATION = 'BAEL91'
       
        piv_a = 10.0E-3
        piv_b = 3.5E-3
        piv_c = 2.0E-3
        nC = 2
        fyd = facier/gammas
        fcd = fbeton*alphacc/gammac

    else if (typco.eq.2) then
!       CALCUL DES PARAMETRES POUR CODIFICATION = 'EC2'

        if (uc.eq.0) then
        unite_pa = 1.e-6
        elseif (uc.eq.1) then
        unite_pa = 1.
        endif
        if (clacier.eq.0) then
        piv_a = 0.9*2.5e-2
        ktys = 1.05
        else if (clacier.eq.1) then
        piv_a = 0.9*5.e-2
        ktys = 1.08
        else
        piv_a = 0.9*7.5e-2
        ktys = 1.15
        endif 
        piv_b = min(3.5E-3,0.26*0.01+3.5*0.01*(((90.d0-fbeton*unite_pa)/100.d0)**4))
        piv_c = 2.0E-3
        if ((fbeton*unite_pa).ge.(50.d0)) then
        piv_c = 0.2*0.01+0.0085*0.01*((fbeton*unite_pa - 50.d0)**(0.53))
        endif
        nC = min(2.0,1.4+23.4*(((90.d0-fbeton*unite_pa)/100.d0)**4))
        fyd = facier/gammas
        fcd = fbeton*alphacc/gammac
    
    endif

    Esu = piv_a
    Euk = Esu/0.9
    Ecu = piv_b
    Ec2 = piv_c
    Ese = fyd/eys
    alphaAB = 1./(1+Esu/Ecu)
    alphaR = 1./(1+Ese/Ecu)
    alphaBC = 1.
    d = ht - enrobi
    d0 = enrobs
    dneg = ht - enrobs
    d0neg = enrobi
    
!   Paramètres de calcul
    Xsup = piv_b/piv_c
    xC = (1-piv_c/piv_b)*ht
    yC = ht-xC
    xCt = xC/ht
    m1 = (((1-xCt)**(nC+1))/(2.d0*(nC+1)))*(1-(2.d0*(1-xCt))/(nC+2))
    m2 = -((1-xCt)**(nC+1))/(nC+1)

!-----------------------------------------------------------------------
!Traitement des différents cas (Pivots A / B / C)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!ET = Entièrement Tendue
!PC = Partiellement Comprimée
!EC = Entièrement Comprimée


!Traitement en PIVOT A - Entièrement Tendu (ET) + Moment Positif
!---------------------------------------------------------------

   N_ET = floor(Esu*1000)+1

   p01 = 'POINT_ITER_DINTELU_01'
   p02 = 'POINT_ITER_DINTELU_02'

   call wkvect(p01, ' V V R ', N_ET, vr=N_P1)
   call wkvect(p02, ' V V R ', N_ET, vr=M_P1)
   
   do k = 1,N_ET
   
      if (k.eq.1) then
      EcSUP = -Esu
      else
      EcSUP = -(1.e-3)*(N_ET-k)
      endif
   
      EsINF = -Esu
      EsSUP = ((EsINF-EcSUP)/d)*(d0)+EcSUP
      EcINF = ((EsINF-EcSUP)/d)*(ht)+EcSUP

      if (typdiag.eq.1) then
      if (Abs(EsINF).lt.Ese) then
      SigmAsINF = eys*(Abs(EsINF))
      else
      SigmAsINF = fyd+((ktys*fyd-fyd)/(Euk-Ese))*(Abs(EsINF)-Ese)
      endif
      if (Abs(EsSUP).lt.Ese) Then
      SigmAsSUP = eys*(Abs(EsSUP))
      else
      SigmAsSUP = fyd+((ktys*fyd-fyd)/(Euk-Ese))*(Abs(EsSUP)-Ese)
      endif
      else
      if (Abs(EsINF).lt.Ese) then
      SigmAsINF = eys*(Abs(EsINF))
      else
      SigmAsINF = fyd
      endif
      if (Abs(EsSUP).lt.Ese) then
      SigmAsSUP = eys*(Abs(EsSUP))
      else
      SigmAsSUP = fyd
      endif
      endif

      if (EsINF.lt.0) then
      SigmAsINF = -SigmAsINF
      endif
      if (EsSUP.lt.0) then
      SigmAsSUP = -SigmAsSUP
      endif

      N_P1(k) = dnsinf*SigmAsINF + dnssup*SigmAsSUP
      M_P1(k) = -dnsinf*SigmAsINF*(d-0.5*ht) + dnssup*SigmAsSUP*(0.5*ht-d0)

   end do   
   
!
!Traitement en PIVOT A et B - Partiellement Comprimée (PC) + Moment Positif
!--------------------------------------------------------------------------
   
   N_PC = ceiling((ht/d)*100)+1

   p03 = 'POINT_ITER_DINTELU_03'
   p04 = 'POINT_ITER_DINTELU_04'

   call wkvect(p03, ' V V R ', N_PC, vr=N_P2)
   call wkvect(p04, ' V V R ', N_PC, vr=M_P2)

   do k=1,N_PC
   
      if (k.lt.N_PC) then
      alpha = (k-1)*0.01
      else
      alpha = ht/d
      endif
   
      if (alpha.lt.alphaAB) then
      EsINF = -Esu
      EcSUP = Esu*alpha/(1-alpha)
      else
      EcSUP = Ecu
      EsINF = -Ecu*(1-alpha)/alpha
      endif

      EsSUP = ((EsINF-EcSUP)/d)*(d0)+EcSUP
      EcINF = ((EsINF-EcSUP)/d)*(ht)+EcSUP

      if (typdiag.eq.1) then
      if (Abs(EsINF).lt.Ese) then
      SigmAsINF = eys*(Abs(EsINF))
      else
      SigmAsINF = fyd+((ktys*fyd-fyd)/(Euk-Ese))*(Abs(EsINF)-Ese)
      endif
      if (Abs(EsSUP).lt.Ese) Then
      SigmAsSUP = eys*(Abs(EsSUP))
      else
      SigmAsSUP = fyd+((ktys*fyd-fyd)/(Euk-Ese))*(Abs(EsSUP)-Ese)
      endif
      else
      if (Abs(EsINF).lt.Ese) then
      SigmAsINF = eys*(Abs(EsINF))
      else
      SigmAsINF = fyd
      endif
      if (Abs(EsSUP).lt.Ese) then
      SigmAsSUP = eys*(Abs(EsSUP))
      else
      SigmAsSUP = fyd
      endif
      endif

      if (EsINF.lt.0) then
      SigmAsINF = -SigmAsINF
      endif
      if (EsSUP.lt.0) then
      SigmAsSUP = -SigmAsSUP
      endif
   
      x1 = (EsINF-EcSUP)/Ec2
      y1 = EcSUP/Ec2
      DELTA = d/ht
   
      if (EcSUP.le.Ec2) then
      Beta = 0
      else
      yE = ((Ec2-EcSUP)/(EsINF-EcSUP))*d
      Beta = yE/d
      endif

      COEF1 = (1-y1-alpha*x1)
      COEF2 = (1-y1-Beta*x1)
      if (abs(COEF1).gt.epsilon(COEF1)) then
      VAR_COEF1 = Abs(COEF1)/COEF1
      else
      VAR_COEF1 = 1
      endif
      if (abs(COEF2).gt.epsilon(COEF2)) then
      VAR_COEF2 = Abs(COEF2)/COEF2
      else
      VAR_COEF2 = 1
      endif

      Ncc = (fcd*bw*d)*(alpha+(1/((nC+1)*x1))*(VAR_COEF1*((Abs(COEF1))**(nC+1)) &
             & -VAR_COEF2*((Abs(COEF2))**(nC + 1))))
      Mcc = (bw*d*d*fcd)*(0.5*alpha*(1/DELTA-alpha) & 
             & +(1/(2*DELTA))*(1/((nC+1)*x1))*(VAR_COEF1*((Abs(COEF1))**(nC+1)) &
             & -VAR_COEF2*((Abs(COEF2))**(nC + 1))) &
             & -(1/((nC+1)*x1))*(alpha*VAR_COEF1*((Abs(COEF1))**(nC + 1)) &
             & -Beta*VAR_COEF2*((Abs(COEF2))**(nC + 1))) &
             & -(1/((nC+1)*(nC+2)*x1*x1))*(VAR_COEF1*((Abs(COEF1))**(nC + 2)) &
             & -VAR_COEF2*((Abs(COEF2))**(nC+2))))

      N_P2(k) = dnsinf*SigmAsINF + dnssup*SigmAsSUP + Ncc
      M_P2(k) = -dnsinf*SigmAsINF*(d-0.5*ht) + dnssup*SigmAsSUP*(0.5*ht-d0) + Mcc
   
   end do


!Traitement en PIVOT C - Entièrement Comprimée (EC) + Moment Positif
!-------------------------------------------------------------------

   Xsup = piv_b/piv_c
   N_EC = ceiling(Xsup*100)+1
   
   p05 = 'POINT_ITER_DINTELU_05'
   p06 = 'POINT_ITER_DINTELU_06'

   call wkvect(p05, ' V V R ', N_EC, vr=N_P3)
   call wkvect(p06, ' V V R ', N_EC, vr=M_P3)

   do k=1,N_EC
       
      X = (N_EC-k)/100.0
      if (k.eq.1) then
      X = Xsup
      endif
       
      DE = X*Ec2
      EcINF = Ec2-DE*(1-xCt)
      EcSUP = DE+EcINF
      EsINF = EcINF+(DE/ht)*(ht-d)
      EsSUP = EcINF+(DE/ht)*(ht-d0)
  
      Ncc = bw*ht*fcd*(1+m2*(X**(nC)))
      Mcc = bw*ht*ht*fcd*m1*(X**(nC))
      Calc = EcSUP-EcINF 
      if (abs(Calc).gt.epsilon(Calc)) then
      alpha = (1/(1-EcINF/EcSUP))*(ht/d)
      else
      alpha = -1000
      endif
   
      if (typdiag.eq.1) then
      if (Abs(EsINF).lt.Ese) then
      SigmAsINF = eys*(Abs(EsINF))
      else
      SigmAsINF = fyd+((ktys*fyd-fyd)/(Euk-Ese))*(Abs(EsINF)-Ese)
      endif
      if (Abs(EsSUP).lt.Ese) Then
      SigmAsSUP = eys*(Abs(EsSUP))
      else
      SigmAsSUP = fyd+((ktys*fyd-fyd)/(Euk-Ese))*(Abs(EsSUP)-Ese)
      endif
      else
      if (Abs(EsINF).lt.Ese) then
      SigmAsINF = eys*(Abs(EsINF))
      else
      SigmAsINF = fyd
      endif
      if (Abs(EsSUP).lt.Ese) then
      SigmAsSUP = eys*(Abs(EsSUP))
      else
      SigmAsSUP = fyd
      endif
      endif

      if (EsINF.lt.0) then
      SigmAsINF = -SigmAsINF
      endif
      if (EsSUP.lt.0) then
      SigmAsSUP = -SigmAsSUP
      endif

      N_P3(k) = dnsinf*SigmAsINF + dnssup*SigmAsSUP + Ncc
      M_P3(k) = -dnsinf*SigmAsINF*(d-0.5*ht) + dnssup*SigmAsSUP*(0.5*ht-d0) + Mcc
   
   end do


!Traitement en PIVOT C - Entièrement Comprimée (EC) + Moment Negatif
!-------------------------------------------------------------------
   
   p07 = 'POINT_ITER_DINTELU_07'
   p08 = 'POINT_ITER_DINTELU_08'

   call wkvect(p07, ' V V R ', N_EC, vr=N_P4)
   call wkvect(p08, ' V V R ', N_EC, vr=M_P4)

   do k=1,N_EC
       
      X = (1-k)/100.0
      if (k.eq.N_EC) then
      X = -Xsup
      endif
       
      DE = X*Ec2
      EcINF = Ec2-DE*(1-xCt)
      EcSUP = DE+EcINF
      EsINF = EcINF+(DE/ht)*(ht-d)
      EsSUP = EcINF+(DE/ht)*(ht-d0)
  
      Ncc = bw*ht*fcd*(1+m2*(X**(nC)))
      Mcc = bw*ht*ht*fcd*m1*(X**(nC))
      Calc = EcSUP-EcINF 
      if (abs(Calc).gt.epsilon(Calc)) then
      alpha = (1/(1-EcINF/EcSUP))*(ht/d)
      else
      alpha = -1000
      endif
   
      if (typdiag.eq.1) then
      if (Abs(EsINF).lt.Ese) then
      SigmAsINF = eys*(Abs(EsINF))
      else
      SigmAsINF = fyd+((ktys*fyd-fyd)/(Euk-Ese))*(Abs(EsINF)-Ese)
      endif
      if (Abs(EsSUP).lt.Ese) Then
      SigmAsSUP = eys*(Abs(EsSUP))
      else
      SigmAsSUP = fyd+((ktys*fyd-fyd)/(Euk-Ese))*(Abs(EsSUP)-Ese)
      endif
      else
      if (Abs(EsINF).lt.Ese) then
      SigmAsINF = eys*(Abs(EsINF))
      else
      SigmAsINF = fyd
      endif
      if (Abs(EsSUP).lt.Ese) then
      SigmAsSUP = eys*(Abs(EsSUP))
      else
      SigmAsSUP = fyd
      endif
      endif

      if (EsINF.lt.0) then
      SigmAsINF = -SigmAsINF
      endif
      if (EsSUP.lt.0) then
      SigmAsSUP = -SigmAsSUP
      endif

      N_P4(k) = dnsinf*SigmAsINF + dnssup*SigmAsSUP + Ncc
      M_P4(k) = -dnsinf*SigmAsINF*(d-0.5*ht) + dnssup*SigmAsSUP*(0.5*ht-d0) - Mcc
   
   end do
   
!Traitement en PIVOT A et B - Partiellement Comprimée (PC) + Moment Negatif
!--------------------------------------------------------------------------
   
   N_PCN = ceiling((ht/dneg)*100)+1

   p09 = 'POINT_ITER_DINTELU_09'
   p10 = 'POINT_ITER_DINTELU_10'

   call wkvect(p09, ' V V R ', N_PCN, vr=N_P5)
   call wkvect(p10, ' V V R ', N_PCN, vr=M_P5)

   do k=1,N_PCN
   
      if (k.gt.1) then
      alpha = (N_PCN-k)*0.01
      else
      alpha = ht/dneg
      endif

      if (alpha.lt.alphaAB) then
      EsSUP = -Esu
      EcINF = Esu*alpha/(1-alpha)
      else
      EcINF = Ecu
      EsSUP = -Ecu*(1-alpha)/alpha
      endif

      EsINF = ((EsSUP-EcINF)/dneg)*(d0neg)+EcINF
      EcSUP = ((EsSUP-EcINF)/dneg)*(ht)+EcINF
      
      if (typdiag.eq.1) then
      if (Abs(EsINF).lt.Ese) then
      SigmAsINF = eys*(Abs(EsINF))
      else
      SigmAsINF = fyd+((ktys*fyd-fyd)/(Euk-Ese))*(Abs(EsINF)-Ese)
      endif
      if (Abs(EsSUP).lt.Ese) Then
      SigmAsSUP = eys*(Abs(EsSUP))
      else
      SigmAsSUP = fyd+((ktys*fyd-fyd)/(Euk-Ese))*(Abs(EsSUP)-Ese)
      endif
      else
      if (Abs(EsINF).lt.Ese) then
      SigmAsINF = eys*(Abs(EsINF))
      else
      SigmAsINF = fyd
      endif
      if (Abs(EsSUP).lt.Ese) then
      SigmAsSUP = eys*(Abs(EsSUP))
      else
      SigmAsSUP = fyd
      endif
      endif

      if (EsINF.lt.0) then
      SigmAsINF = -SigmAsINF
      endif
      if (EsSUP.lt.0) then
      SigmAsSUP = -SigmAsSUP
      endif
   
      x1 = (EsSUP-EcINF)/Ec2
      y1 = EcINF/Ec2
      DELTA = dneg/ht
   
      if (EcINF.le.Ec2) then
      Beta = 0
      else
      yE = ((Ec2-EcINF)/(EsSUP-EcINF))*dneg
      Beta = yE/dneg
      endif

      COEF1 = (1-y1-alpha*x1)
      COEF2 = (1-y1-Beta*x1)
      if (abs(COEF1).gt.epsilon(COEF1)) then
      VAR_COEF1 = Abs(COEF1)/COEF1
      else
      VAR_COEF1 = 1
      endif
      if (abs(COEF2).gt.epsilon(COEF2)) then
      VAR_COEF2 = Abs(COEF2)/COEF2
      else
      VAR_COEF2 = 1
      endif

      Ncc = (fcd*bw*dneg)*(alpha+(1/((nC+1)*x1))*(VAR_COEF1*((Abs(COEF1))**(nC+1)) &
             & -VAR_COEF2*((Abs(COEF2))**(nC + 1))))
      Mcc = (bw*dneg*dneg*fcd)*(0.5*alpha*(1/DELTA-alpha) & 
             & +(1/(2*DELTA))*(1/((nC+1)*x1))*(VAR_COEF1*((Abs(COEF1))**(nC+1)) &
             & -VAR_COEF2*((Abs(COEF2))**(nC + 1))) &
             & -(1/((nC+1)*x1))*(alpha*VAR_COEF1*((Abs(COEF1))**(nC + 1)) &
             & -Beta*VAR_COEF2*((Abs(COEF2))**(nC + 1))) &
             & -(1/((nC+1)*(nC+2)*x1*x1))*(VAR_COEF1*((Abs(COEF1))**(nC + 2)) &
             & -VAR_COEF2*((Abs(COEF2))**(nC+2))))

      N_P5(k) = dnsinf*SigmAsINF + dnssup*SigmAsSUP + Ncc
      M_P5(k) = -dnsinf*SigmAsINF*(d-0.5*ht) + dnssup*SigmAsSUP*(0.5*ht-d0) - Mcc      
   
   end do

!Traitement en PIVOT A - Entièrement Tendu (ET) + Moment Negatif
!---------------------------------------------------------------

   p11 = 'POINT_ITER_DINTELU_11'
   p12 = 'POINT_ITER_DINTELU_12'

   call wkvect(p11, ' V V R ', N_ET, vr=N_P6)
   call wkvect(p12, ' V V R ', N_ET, vr=M_P6)
   
   do k = 1,N_ET
   
      if (k.eq.N_ET) then
      EcINF = -Esu
      else
      EcINF = -0.1*(k-1)
      endif
   
      EsSUP = -Esu
      EsINF = ((EsSUP-EcINF)/dneg)*(d0neg)+EcINF
      EcSUP = ((EsSUP-EcINF)/dneg)*(ht)+EcINF

      if (typdiag.eq.1) then
      if (Abs(EsINF).lt.Ese) then
      SigmAsINF = eys*(Abs(EsINF))
      else
      SigmAsINF = fyd+((ktys*fyd-fyd)/(Euk-Ese))*(Abs(EsINF)-Ese)
      endif
      if (Abs(EsSUP).lt.Ese) Then
      SigmAsSUP = eys*(Abs(EsSUP))
      else
      SigmAsSUP = fyd+((ktys*fyd-fyd)/(Euk-Ese))*(Abs(EsSUP)-Ese)
      endif
      else
      if (Abs(EsINF).lt.Ese) then
      SigmAsINF = eys*(Abs(EsINF))
      else
      SigmAsINF = fyd
      endif
      if (Abs(EsSUP).lt.Ese) then
      SigmAsSUP = eys*(Abs(EsSUP))
      else
      SigmAsSUP = fyd
      endif
      endif

      if (EsINF.lt.0) then
      SigmAsINF = -SigmAsINF
      endif
      if (EsSUP.lt.0) then
      SigmAsSUP = -SigmAsSUP
      endif

      N_P6(k) = dnsinf*SigmAsINF + dnssup*SigmAsSUP
      M_P6(k) = -dnsinf*SigmAsINF*(d-0.5*ht) + dnssup*SigmAsSUP*(0.5*ht-d0)

   end do   

!-----------------------------------------------------------------------
!Fin de Traitement des différents cas
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   
   do k=1,N_ET
      nrd(k) = N_P1(k)
      mrd(k) = M_P1(k)
   end do
   do k=1,N_PC
      nrd(k+N_ET) = N_P2(k)
      mrd(k+N_ET) = M_P2(k)
   end do
   do k=1,N_EC
      nrd(k+N_ET+N_PC) = N_P3(k)
      mrd(k+N_ET+N_PC) = M_P3(k)
   end do
   do k=1,N_EC
      nrd(k+N_ET+N_PC+N_EC) = N_P4(k)
      mrd(k+N_ET+N_PC+N_EC) = M_P4(k)
   end do
   do k=1,N_PCN
      nrd(k+N_ET+N_PC+N_EC+N_EC) = N_P5(k)
      mrd(k+N_ET+N_PC+N_EC+N_EC) = M_P5(k)
   end do
   do k=1,N_ET
      nrd(k+N_ET+N_PC+N_EC+N_EC+N_PCN) = N_P6(k)
      mrd(k+N_ET+N_PC+N_EC+N_EC+N_PCN) = M_P6(k)
   end do
   
   call jedetr(p01)
   call jedetr(p02)
   call jedetr(p03)
   call jedetr(p04)
   call jedetr(p05)
   call jedetr(p06)
   call jedetr(p07)
   call jedetr(p08)
   call jedetr(p09)
   call jedetr(p10)
   call jedetr(p11)
   call jedetr(p12)

end subroutine
