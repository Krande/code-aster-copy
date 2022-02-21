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

subroutine dintels(cequi, ht, bw, enrobi, enrobs,&
                   scmaxi, scmaxs, ssmax, uc,&
                   dnsinf, dnssup, ntot, nrd, mrd)
!______________________________________________________________________
!
!      DINTELS

!      CONSTRUCTION DU DIAGRAMME D'INTERACTION D'UNE SECTION
!      FERRAILLEE - VERIFICATION D'UN FERRAILLAGE EXISTANT
!      CRITERE = LIMITATION DES CONTRAINTES (ELS)
!
!      I CEQUI     COEFFICIENT D'EQUIVALENCE ACIER/BETON
!      I HT        HAUTEUR DE LA SECTION
!      I BW        LARGEUR DE LA SECTION
!      I ENROBI    ENROBAGE DES ARMATURES INFERIEURES
!      I ENROBS    ENROBAGE DES ARMATURES SUPERIEURES
!      I SCMAXI    CONTRAINTE DE COMPRESSION MAXI DU BETON EN FIBRE INF
!      I SCMAXS    CONTRAINTE DE COMPRESSION MAXI DU BETON EN FIBRE SUP
!      I SSMAX     CONTRAINTE MAXI DE L'ACIER DE FLEXION
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
    real(kind=8) :: cequi
    real(kind=8) :: ht
    real(kind=8) :: bw
    real(kind=8) :: enrobi
    real(kind=8) :: enrobs
    real(kind=8) :: scmaxi
    real(kind=8) :: scmaxs
    real(kind=8) :: ssmax
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
    real(kind=8) :: unite_pa
    real(kind=8) :: X,scmax,scmaxneg
    real(kind=8) :: alpha_12,alpha
    real(kind=8) :: ScSUP,ScINF
    real(kind=8) :: SsSUP,SsINF,Ncc,Mcc
    integer :: N_ET,N_PC,N_EC,N_PCAC,N_PCACN,k
    integer :: N_ECN

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

!   Paramètres de calcul

    if (uc.eq.0) then
    unite_pa = 1.e-6
    else if (uc.eq.1) then
    unite_pa = 1.
    endif

    d = ht - enrobi
    d0 = enrobs
    scmax = scmaxs
    dneg = ht - enrobs
    d0neg = enrobi
    scmaxneg = scmaxi

    alpha_12 = 1.0/(1.0+(ssmax/cequi)/scmax)
    
!   Initialisation des entiers
    
    N_ET = 100
    N_PC = 100
    N_EC = 100
    N_PCAC = 100
    N_PCACN = 100
    N_ECN = 100 
    k = 1

!-----------------------------------------------------------------------
!Traitement des différents cas (Pivots A / B / C)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!ET = Entièrement Tendue
!PC = Partiellement Comprimée
!EC = Entièrement Comprimée


!Traitement en PIVOT A - Entièrement Tendu (ET) + Moment Positif
!---------------------------------------------------------------

   N_ET = 11

   p01 = 'POINT_ITER_DINTELS_01'
   p02 = 'POINT_ITER_DINTELS_02'

   call wkvect(p01, ' V V R ', N_ET, vr=N_P1)
   call wkvect(p02, ' V V R ', N_ET, vr=M_P1)
   
   do k = 1,N_ET
   
      SsINF = -ssmax
      ScSUP = -(ssmax/cequi)*(1-0.1*(k-1))
      SsSUP = ((SsINF/cequi-ScSUP)*((ht-d0)/d)+ScSUP)*cequi
      ScINF = (SsINF/cequi-ScSUP)*(ht/d)+ScSUP
      N_P1(k) = dnsinf*SsINF + dnssup*SsSUP
      M_P1(k) = -dnsinf*SsINF*(d-0.5*ht) + dnssup*SsSUP*(0.5*ht-d0)

   end do

!Traitement en PIVOT A et B - Partiellement Comprimée (PC) + Moment Positif
!--------------------------------------------------------------------------
   
   N_PC = 101
   N_PCAC = CEILING((N_PC-1)*(ht/d))+1

   p03 = 'POINT_ITER_DINTELS_03'
   p04 = 'POINT_ITER_DINTELS_04'

   call wkvect(p03, ' V V R ', N_PCAC, vr=N_P2)
   call wkvect(p04, ' V V R ', N_PCAC, vr=M_P2)

   do k=1,N_PCAC
   
      if (k.lt.N_PCAC) then
      alpha = real(k-1)/real(N_PC-1)
      else
      alpha = ht/d
      endif
      X = alpha*d
      if (alpha.eq.0) then
      SsINF = -ssmax
      ScSUP = 0
      ScINF = -(ssmax/cequi)*(ht/d)
      SsSUP = -ssmax*(ht-d0)/d
      elseif ((alpha.gt.0) .AND. (alpha.lt.alpha_12)) then
      SsINF = -ssmax
      ScSUP = (X/(d-X))*(ssmax/cequi)
      ScINF = ScSUP*(1-ht/X)
      SsSUP = ScSUP*(1-(ht-d0)/X)*cequi
      else
      ScSUP = scmax
      SsINF = scmax*cequi*(1-d/X)
      ScINF = scmax*(1-ht/X)
      SsSUP = scmax*cequi*(1-(ht-d0)/X)
      endif
      Ncc = ScSUP*0.5*X*bw
      Mcc = ScSUP*((1./4.)*ht*X - (1./6.)*X*X)*bw
      N_P2(k) = dnsinf*SsINF + dnssup*SsSUP + Ncc
      M_P2(k) = -dnsinf*SsINF*(d-0.5*ht) + dnssup*SsSUP*(0.5*ht-d0) + Mcc
   
   end do

!Traitement en PIVOT C - Entièrement Comprimée (EC) + Moment Positif
!-------------------------------------------------------------------

   N_EC = CEILING(10*(scmax*unite_pa))+1
   
   p05 = 'POINT_ITER_DINTELS_05'
   p06 = 'POINT_ITER_DINTELS_06'

   call wkvect(p05, ' V V R ', N_EC, vr=N_P3)
   call wkvect(p06, ' V V R ', N_EC, vr=M_P3)

   do k=1,N_EC
       
      ScINF = scmax*(real(k-1)/real(N_EC-1))
      ScSUP = scmax
      if (k.lt.N_EC) then
      X = (ScSUP/(ScSUP-ScINF))*ht
      alpha = X/d
      else
      alpha = -1000
      endif
      Ncc = 0.5*(ScSUP + ScINF)*ht*bw
      Mcc = (1./12.)*(ScSUP - ScINF)*ht*ht*bw
      SsSUP = ((ScINF - ScSUP)*(ht-d0)/ht + ScSUP)*cequi
      SsINF = ((ScINF - ScSUP)*d/ht + ScSUP)*cequi
      N_P3(k) = dnsinf*SsINF + dnssup*SsSUP + Ncc
      M_P3(k) = -dnsinf*SsINF*(d-0.5*ht) + dnssup*SsSUP*(0.5*ht-d0) + Mcc
   
   end do

!Traitement en PIVOT C - Entièrement Comprimée (EC) + Moment Negatif
!-------------------------------------------------------------------

   N_ECN = CEILING(10*(scmaxneg*unite_pa))+1
   
   p07 = 'POINT_ITER_DINTELS_07'
   p08 = 'POINT_ITER_DINTELS_08'

   call wkvect(p07, ' V V R ', N_ECN, vr=N_P4)
   call wkvect(p08, ' V V R ', N_ECN, vr=M_P4)

   do k=1,N_ECN
       
      ScSUP = scmaxneg*(1-real(k-1)/real(N_ECN-1))
      ScINF = scmaxneg
      if (k.ne.1) then
      X = (ScSUP/(ScSUP-ScINF))*ht
      alpha = X/d
      else
      alpha = -1000
      endif
      Ncc = 0.5*(ScSUP + ScINF)*ht*bw
      Mcc = (1./12.)*(ScINF - ScSUP)*ht*ht*bw
      SsSUP = ((ScINF - ScSUP)*(ht-d0)/ht + ScSUP)*cequi
      SsINF = ((ScINF - ScSUP)*d/ht + ScSUP)*cequi
      N_P4(k) = dnsinf*SsINF + dnssup*SsSUP + Ncc
      M_P4(k) = -dnsinf*SsINF*(d-0.5*ht) + dnssup*SsSUP*(0.5*ht-d0) - Mcc
   
   end do

!Traitement en PIVOT A et B - Partiellement Comprimée (PC) + Moment Negatif
!--------------------------------------------------------------------------

   N_PCACN = CEILING((N_PC-1)*(ht/dneg))+1

   p09 = 'POINT_ITER_DINTELS_09'
   p10 = 'POINT_ITER_DINTELS_10'

   call wkvect(p09, ' V V R ', N_PCACN, vr=N_P5)
   call wkvect(p10, ' V V R ', N_PCACN, vr=M_P5)

   do k=1,N_PCACN
   
      if (k.gt.1) then
      alpha = 1-real(k-1)/real(N_PC-1)
      else
      alpha = ht/dneg
      endif
      X = alpha*dneg
      if (alpha.eq.0) then
      SsSUP = -ssmax
      ScINF = 0
      ScSUP = -(ssmax/cequi)*(ht/dneg)
      SsINF = -ssmax*(ht-d0neg)/dneg
      elseif ((alpha.gt.0) .AND. (alpha.lt.alpha_12)) then
      SsSUP = -ssmax
      ScINF = (X/(dneg-X))*(ssmax/cequi)
      ScSUP = ScINF*(1-ht/X)
      SsINF = ScINF*(1-(ht-d0neg)/X)*cequi
      else
      ScINF = scmaxneg
      SsSUP = scmaxneg*cequi*(1-dneg/X)
      ScSUP = scmaxneg*(1-ht/X)
      SsINF = scmaxneg*cequi*(1-(ht-d0neg)/X)
      endif
      Ncc = ScINF*0.5*X*bw
      Mcc = ScINF*((1./4.)*ht*X - (1./6.)*X*X)*bw
      N_P5(k) = dnsinf*SsINF + dnssup*SsSUP + Ncc
      M_P5(k) = -dnsinf*SsINF*(d-0.5*ht) + dnssup*SsSUP*(0.5*ht-d0) - Mcc
   
   end do

!Traitement en PIVOT A - Entièrement Tendu (ET) + Moment Negatif
!---------------------------------------------------------------

   p11 = 'POINT_ITER_DINTELS_11'
   p12 = 'POINT_ITER_DINTELS_12'

   call wkvect(p11, ' V V R ', N_ET, vr=N_P6)
   call wkvect(p12, ' V V R ', N_ET, vr=M_P6)

   do k = 1,N_ET
   
      SsSUP = -ssmax
      ScINF = -(ssmax/cequi)*0.1*(k-1)
      SsINF = ((SsSUP/cequi-ScINF)*((ht-d0neg)/dneg)+ScINF)*cequi
      ScINF = (SsSUP/cequi-ScINF)*(ht/dneg)+ScINF
      N_P6(k) = dnsinf*SsINF + dnssup*SsSUP
      M_P6(k) = -dnsinf*SsINF*(d-0.5*ht) + dnssup*SsSUP*(0.5*ht-d0)

   end do   

!-----------------------------------------------------------------------
!Fin de Traitement des différents cas
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

   do k=1,N_ET
      nrd(k) = N_P1(k)
      mrd(k) = M_P1(k)
   end do
   do k=1,N_PCAC
      nrd(k+N_ET) = N_P2(k)
      mrd(k+N_ET) = M_P2(k)
   end do
   do k=1,N_EC
      nrd(k+N_ET+N_PCAC) = N_P3(k)
      mrd(k+N_ET+N_PCAC) = M_P3(k)
   end do
   do k=1,N_ECN
      nrd(k+N_ET+N_PCAC+N_EC) = N_P4(k)
      mrd(k+N_ET+N_PCAC+N_EC) = M_P4(k)
   end do
   do k=1,N_PCACN
      nrd(k+N_ET+N_PCAC+N_EC+N_ECN) = N_P5(k)
      mrd(k+N_ET+N_PCAC+N_EC+N_ECN) = M_P5(k)
   end do
   do k=1,N_ET
      nrd(k+N_ET+N_PCAC+N_EC+N_ECN+N_PCACN) = N_P6(k)
      mrd(k+N_ET+N_PCAC+N_EC+N_ECN+N_PCACN) = M_P6(k)
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
