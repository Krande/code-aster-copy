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

subroutine bresels(cequi, effmy, effmz, effn,&
                   ht, bw, enrobyi, enrobys, enrobzi, enrobzs,&
                   scmaxyi, scmaxys, scmaxzi, scmaxzs, ssmax,&
                   ferrcomp, ferrsyme, slsyme, uc,&
                   dnsyi, dnsys, dnszi, dnszs,& 
                   sigmsyi, sigmsys, sigmcyi, sigmcys,&
                   sigmszi, sigmszs, sigmczi, sigmczs,&
                   alphay, alphaz, pivoty, pivotz, etaty, etatz, ierr)

!______________________________________________________________________
!
!      BRESELS

!      CALCUL DES ACIERS EN FLEXION COMPOSEE DEVIEE A L'ELS
!      METHODE DE VERIFICATION DE BRESLER
!      CRITERE = LIMITATION DES CONTRAINTES
!
!      I CEQUI     COEFFICIENT D'EQUIVALENCE ACIER/BETON
!      I EFFMY     MOMENT DE FLEXION SUIVANT L'AXE Y
!      I EFFMZ     MOMENT DE FLEXION SUIVANT L'AXE Z
!      I EFFN      EFFORT NORMAL
!      I HT        HAUTEUR DE LA SECTION
!      I BW        LARGEUR DE LA SECTION
!      I ENROBYI   ENROBAGE DES ARMATURES INF SUIVANT L'AXE Y
!      I ENROBYS   ENROBAGE DES ARMATURES SUP SUIVANT L'AXE Y
!      I ENROBZI   ENROBAGE DES ARMATURES INF SUIVANT L'AXE Z
!      I ENROBZS   ENROBAGE DES ARMATURES SUP SUIVANT L'AXE Z
!      I SCMAXYI   CONT DE COMP MAX DU BETON EN FIBRE INF SUIVANT L'AXE Y
!      I SCMAXYS   CONT DE COMP MAX DU BETON EN FIBRE SUP SUIVANT L'AXE Y
!      I SCMAXZI   CONT DE COMP MAX DU BETON EN FIBRE INF SUIVANT L'AXE Z
!      I SCMAXZS   CONT DE COMP MAX DU BETON EN FIBRE SUP SUIVANT L'AXE Z
!      I SSMAX     CONTRAINTE MAXI DE L'ACIER DE FLEXION
!      I FERRCOMP  PRISE EN COMPTE DU FERRAILLAGE DE COMPRESSION
!                     FERRCOMP = 1 (NON)
!                     FERRCOMP = 2 (OUI)
!      I FERRSYME   FERRAILLAGE SYMETRIQUE?
!                     FERRSYME = 0 (NON)
!                     FERRSYME = 1 (OUI)
!      I SLSYME    SECTION SEUIL DE TOLERANCE POUR UN FERRAILLAGE SYMETRIQUE 
!      I UC        UNITE DES CONTRAINTES :
!                     UC = 0 CONTRAINTES EN Pa
!                     UC = 1 CONTRAINTES EN MPa
!
!      O DNSYI     DENSITE DE L'ACIER INF SUIVANT L'AXE Y
!      O DNSYS     DENSITE DE L'ACIER SUP SUIVANT L'AXE Y
!      O DNSZI     DENSITE DE L'ACIER INF SUIVANT L'AXE Z
!      O DNSZS     DENSITE DE L'ACIER SUP SUIVANT L'AXE Z
!      O SIGMSYI   CONTRAINTE DANS L'ACIER INF SUIVANT L'AXE Y
!      O SIGMSYS   CONTRAINTE DANS L'ACIER SUP SUIVANT L'AXE Y
!      O SIGMCYI   CONTRAINTE EN FIBRE INF BETON SUIVANT L'AXE Y
!      O SIGMCYS   CONTRAINTE EN FIBRE SUP BETON SUIVANT L'AXE Y
!      O SIGMSZI   CONTRAINTE DANS L'ACIER INF SUIVANT L'AXE Z
!      O SIGMSZS   CONTRAINTE DANS L'ACIER SUP SUIVANT L'AXE Z
!      O SIGMCZI   CONTRAINTE EN FIBRE INF BETON SUIVANT L'AXE Z
!      O SIGMCZS   CONTRAINTE EN FIBRE SUP BETON SUIVANT L'AXE Z
!      O ALPHAY    COEF DE PROFONDEUR DE L'AN SUIVANT L'AXE Y
!      O ALPHAZ    COEF DE PROFONDEUR DE L'AN SUIVANT L'AXE Z
!      O PIVOTY    PIVOT EN FC SUIVANT L'AXE Y
!      O PIVOTZ    PIVOT EN FC SUIVANT L'AXE Z
!      O ETATY     ETAT EN FC SUIVANT L'AXE Y
!      O ETATZ     ETAT EN FC SUIVANT L'AXE Z
!      O IERR     CODE RETOUR (0 = OK)
!
!______________________________________________________________________
!
!
    implicit none
!
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/jedetr.h"
#include "asterfort/juveca.h"
#include "asterfort/jeveuo.h"
#include "asterfort/cafels.h"
#include "asterfort/dintels.h"
!
    real(kind=8) :: cequi
    real(kind=8) :: effmy
    real(kind=8) :: effmz
    real(kind=8) :: effn
    real(kind=8) :: ht
    real(kind=8) :: bw
    real(kind=8) :: enrobyi
    real(kind=8) :: enrobys
    real(kind=8) :: enrobzi
    real(kind=8) :: enrobzs
    real(kind=8) :: scmaxyi
    real(kind=8) :: scmaxys
    real(kind=8) :: scmaxzi
    real(kind=8) :: scmaxzs
    real(kind=8) :: ssmax
    integer :: ferrcomp
    integer :: ferrsyme
    real(kind=8) :: slsyme
    integer :: uc
    real(kind=8) :: dnsyi
    real(kind=8) :: dnsys
    real(kind=8) :: dnszi
    real(kind=8) :: dnszs
    real(kind=8) :: sigmsyi
    real(kind=8) :: sigmsys
    real(kind=8) :: sigmcyi
    real(kind=8) :: sigmcys
    real(kind=8) :: sigmszi
    real(kind=8) :: sigmszs
    real(kind=8) :: sigmczi
    real(kind=8) :: sigmczs
    real(kind=8) :: alphay
    real(kind=8) :: alphaz
    integer :: pivoty
    integer :: pivotz
    integer :: etaty
    integer :: etatz
    integer :: ierr

!-----------------------------------------------------------------------
!!!!VARIABLES DE CALCUL
!-----------------------------------------------------------------------
   real(kind=8) :: Acc,fcd,fyd,coeff,Ass,Aiter
   real(kind=8) :: rhoyinf,rhoysup,rhozinf,rhozsup
   real(kind=8) :: BRES,mrdyE,mrdy1,mrdy2,mrdzE,mrdz1,mrdz2,nrdyzE,a,nrd0
   logical :: COND
   integer :: s, COUNT_BRES
   real(kind=8), pointer :: nrdy(:) => null(), mrdy(:) => null()
   real(kind=8), pointer :: nrdz(:) => null(), mrdz(:) => null()
   character(24) :: pnrdy, pmrdy, pnrdz, pmrdz
   real(kind=8) :: unite_pa, d, d0, dneg, d0neg, scmax, scmaxneg
   integer :: N_ET, N_PC, N_PCAC, N_EC, N_ECN, N_PCACN
   integer :: ntoty, ndemiy, ntotz, ndemiz

   pnrdy = 'POINT_NRD_Y'
   pmrdy = 'POINT_MRD_Y'
   pnrdz = 'POINT_NRD_Z'
   pmrdz = 'POINT_MRD_Z'
   
   Acc = bw*ht
   fcd = (scmaxyi+scmaxys+scmaxzi+scmaxzs)/4
   fyd = ssmax
   
   !Initialisation
   ntoty = 1
   ndemiy = 1
   ntotz = 1
   ndemiz = 1   
   mrdy1 = -1.0
   mrdy2 = -1.0
   mrdz1 = -1.0
   mrdz2 = -1.0
   mrdyE = -1.0
   mrdzE = -1.0
   nrdyzE = -1.0
   nrd0 = -1.0
   s = 1
   
   !Effort Axial uniquement
   !if ((effmy.eq.0) .and. (effmz.eq.0) .and. (effn.ne.0)) then
   if ((abs(effmy).lt.epsilon(effmy)) .and. (abs(effmz).lt.epsilon(effmz))) then
      call cafels(cequi, effmy, 0.5*effn, ht, bw,&
                  enrobyi, enrobys, scmaxys, scmaxyi, ssmax,& 
                  ferrcomp, ferrsyme, slsyme, uc,&
                  dnsyi, dnsys, sigmsyi, sigmsys,&
                  sigmcyi, sigmcys,&
                  alphay, pivoty, etaty, ierr)
      if (ierr.ne.0) then 
          goto 998
      endif
      call cafels(cequi, effmz, 0.5*effn, bw, ht,&
                  enrobzi, enrobzs, scmaxzs, scmaxzi, ssmax,& 
                  ferrcomp, ferrsyme, slsyme, uc,&
                  dnszi, dnszs, sigmszi, sigmszs,&
                  sigmczi, sigmczs,&
                  alphaz, pivotz, etatz, ierr)
      if (ierr.ne.0) then
          goto 998
      endif
      
   else
   
      !Calcul suivant "y"
      !if (effmy.ne.0) then
      if (abs(effmy).gt.epsilon(effmy)) then
         call cafels(cequi, effmy, effn, ht, bw,&
                     enrobyi, enrobys, scmaxys, scmaxyi, ssmax,& 
                     ferrcomp, ferrsyme, slsyme, uc,&
                     dnsyi, dnsys, sigmsyi, sigmsys,&
                     sigmcyi, sigmcys,&
                     alphay, pivoty, etaty, ierr)
         if (ierr.ne.0) then 
             goto 998
         endif
      endif
      
      !Calcul suivant "z"
      !if (effmz.ne.0) then
      if (abs(effmz).gt.epsilon(effmz)) then
         call cafels(cequi, effmz, effn, bw, ht,&
                     enrobzi, enrobzs, scmaxzs, scmaxzi, ssmax,& 
                     ferrcomp, ferrsyme, slsyme, uc,&
                     dnszi, dnszs, sigmszi, sigmszs,&
                     sigmczi, sigmczs,&
                     alphaz, pivotz, etatz, ierr)
         if (ierr.ne.0) then
             goto 998
         endif
      endif
   
   endif
   
   !if ((effmy.ne.0) .and. (effmz.ne.0)) then
   if ((abs(effmy).gt.epsilon(effmy)) .and. (abs(effmz).gt.epsilon(effmz))) then

   !Iteration Bresler
        COND = .false.
        COUNT_BRES = 0
        BRES = 1.5
        
   !Dimensionnement des vecteurs

        if (uc.eq.0) then
        unite_pa = 1.e-6
        elseif (uc.eq.1) then
        unite_pa = 1.
        endif
        
        N_ET = 11
        N_PC = 101
        
        d = ht - enrobyi
        d0 = enrobys
        scmax = scmaxys
        dneg = ht - enrobys
        d0neg = enrobyi
        scmaxneg = scmaxyi

        N_PCAC = CEILING((N_PC-1)*(ht/d))+1
        N_EC = CEILING(10*(scmax*unite_pa))+1
        N_ECN = CEILING(10*(scmaxneg*unite_pa))+1
        N_PCACN = CEILING((N_PC-1)*(ht/dneg))+1
        
        ntoty = N_ET + N_PCAC + N_EC + N_ECN + N_PCACN + N_ET
        ndemiy = N_ET + N_PCac + N_EC
        call wkvect(pnrdy, ' V V R ', ntoty, vr=nrdy)
        call wkvect(pmrdy, ' V V R ', ntoty, vr=mrdy)
        
        d = bw - enrobzi
        d0 = enrobzs
        scmax = scmaxzs
        dneg = bw - enrobzs
        d0neg = enrobzi
        scmaxneg = scmaxzi

        N_PCAC = CEILING((N_PC-1)*(bw/d))+1
        N_EC = CEILING(10*(scmax*unite_pa))+1
        N_ECN = CEILING(10*(scmaxneg*unite_pa))+1
        N_PCACN = CEILING((N_PC-1)*(bw/dneg))+1
        
        ntotz = N_ET + N_PCAC + N_EC + N_ECN + N_PCACN + N_ET
        ndemiz = N_ET + N_PCac + N_EC
        call wkvect(pnrdz, ' V V R ', ntotz, vr=nrdz)
        call wkvect(pmrdz, ' V V R ', ntotz, vr=mrdz)
        
        
        do while (COND.eqv.(.false.))
        
           Ass = dnsyi+dnsys+dnszi+dnszs
           nrdyzE= Acc*fcd + Ass*fyd

   !Determiner MRd,y
           
           do s=1,ntoty
           nrdy(s) = -1.0
           mrdy(s) = -1.0
           end do

           call dintels(cequi, ht, bw, enrobyi, enrobys,&
                        scmaxyi, scmaxys, ssmax, uc,&
                        dnsyi, dnsys, ntoty, nrdy, mrdy)

           s = 1
           nrd0 = nrdy(s)
           do while ((nrd0.lt.effn) .and. (s.lt.ndemiy))
                      s = s + 1
                      nrd0 = nrdy(s)
           end do
           if ((s.eq.1) .or. (s.eq.ndemiy)) then
                BRES = 1.5
                goto 999
           else
                mrdy1 = ((mrdy(s)-mrdy(s-1))/(nrdy(s)-nrdy(s-1)))*(effn-nrdy(s-1))+mrdy(s-1)
           endif
           s = ndemiy+1
           nrd0 = nrdy(s)
           do while ((nrd0.gt.effn) .and. (s.lt.ntoty))
                      s = s + 1
                      nrd0 = nrdy(s)
           end do
           if ((s.eq.(ndemiy+1)) .or. (s.eq.ntoty)) then
                BRES = 1.5
                goto 999
           else
                mrdy2 = ((mrdy(s)-mrdy(s-1))/(nrdy(s)-nrdy(s-1)))*(effn-nrdy(s-1))+mrdy(s-1)
           endif
           if (effmy.gt.0.0) then
               mrdy1 = max(mrdy1,0.0)
               mrdy2 = max(mrdy2,0.0)
               mrdyE = max(mrdy1,mrdy2)
           elseif (effmy.lt.0.0) then
               mrdy1 = min(mrdy1,0.0)
               mrdy2 = min(mrdy2,0.0)
               mrdyE = min(mrdy1,mrdy2)
           endif

   !Determiner MRd,z
        
           call dintels(cequi, bw, ht, enrobzi, enrobzs,&
                        scmaxzi, scmaxzs, ssmax, uc,&
                        dnszi, dnszs, ntotz, nrdz, mrdz)

           s = 1
           nrd0 = nrdz(s)
           do while ((nrd0.lt.effn) .and. (s.lt.ndemiz))
                      s = s + 1
                      nrd0 = nrdz(s)
           end do
           if ((s.eq.1) .or. (s.eq.ndemiz)) then
                BRES = 1.5
                goto 999
           else
                mrdz1 = ((mrdz(s)-mrdz(s-1))/(nrdz(s)-nrdz(s-1)))*(effn-nrdz(s-1))+mrdz(s-1)
           endif
           s = ndemiz+1
           nrd0 = nrdz(s)
           do while ((nrd0.gt.effn) .and. (s.lt.ntotz))
                      s = s + 1
                      nrd0 = nrdz(s)
           end do
           if ((s.eq.(ndemiz+1)) .or. (s.eq.ntotz)) then
                BRES = 1.5
                goto 999
           else
                mrdz2 = ((mrdz(s)-mrdz(s-1))/(nrdz(s)-nrdz(s-1)))*(effn-nrdz(s-1))+mrdz(s-1)
           endif
           if (effmz.gt.0.0) then
               mrdz1 = max(mrdz1,0.0)
               mrdz2 = max(mrdz2,0.0)
               mrdzE = max(mrdz1,mrdz2)
           elseif (effmz.lt.0.0) then
               mrdz1 = min(mrdz1,0.0)
               mrdz2 = min(mrdz2,0.0)
               mrdzE = min(mrdz1,mrdz2)
           endif

   !Calcul de 'a'

           if (abs(nrdyzE).gt.epsilon(nrdyzE)) then
               coeff = effn/nrdyzE
           else
               coeff = 0.0
           endif
           if (coeff.le.0.1) then
               a = 1.0
           elseif (coeff.le.0.7) then
               a = ((1.5-1.0)/(0.7-0.1))*(coeff-0.1) + 1.0
           elseif (coeff.le.1.0) then
               a = ((2.0-1.5)/(1.0-0.7))*(coeff-0.7) + 1.5
           else
               a = 2.0
           endif

   !Calcul de 'BRES'
           if ((abs(mrdyE).gt.epsilon(mrdyE)) .and. (abs(mrdzE).gt.epsilon(mrdzE))) then
           BRES = (effmy/mrdyE)**(a) + (effmz/mrdzE)**(a)
           endif
      
   !Verif de 'BRES' et iteration
           999 continue
           
           COUNT_BRES = COUNT_BRES + 1
           if (BRES.gt.1) then
           Aiter = 0.1*Ass
           rhoyinf = dnsyi/Ass
           rhoysup = dnsys/Ass
           rhozinf = dnszi/Ass
           rhozsup = dnszs/Ass
           dnsyi = dnsyi + rhoyinf*Aiter
           dnsys = dnsys + rhoysup*Aiter
           dnszi = dnszi + rhozinf*Aiter
           dnszs = dnszs + rhozsup*Aiter
              if (COUNT_BRES.eq.100) then
              ierr = 4
              dnsyi = -1
              dnsys = -1
              dnszi = -1
              dnszs = -1
              COND = .true.
              endif
           else
           COND = .true.
           endif

    end do
    !do while (COND.eqv.(.false.))
        
        call jedetr(pnrdy)
        call jedetr(pmrdy)
        call jedetr(pnrdz)
        call jedetr(pmrdz)
        
   endif
   !if ((effmy.ne.0) .and. (effmz.ne.0)) then
   
998 continue

end subroutine
