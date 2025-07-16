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

  subroutine vrfplq(typcmb, typco, nb, cequi, enrobi, enrobs, sigs, sigci, &
                    sigcs, alphacc, gammas, gammac, facier, eys, typdiag, fbeton, clacier, &
                    uc, ht, effrts, effref, dnsxi, dnsxs, dnsyi, dnsys, marge, theta_crit, &
                    dnsinf_crit, dnssup_crit, effn0_crit, effm0_crit, effn_crit, effm_crit, &
                    myNrd_crit, myMrd_crit, bw, tau_crit, c0c_crit, c0crd_crit)
!_____________________________________________________________________
!
!      VRFPLQ
!
!      VERIFICATION DE FERRAILLAGE D UN ELEMENT DE PLAQUE
!
!      I TYPCMB        TYPE DE COMBINAISON (0 = ELU, 1 = ELS)
!      I TYPCO         CODIFICATION UTILISEE (1 = BAEL91, 2 = EC2)
!      I NB            NOMBRE DE FACETTES
!      I CEQUI         COEFFICIENT D'EQUIVALENCE ACIER/BETON
!      I ENROBI        ENROBAGE DES ARMATURES INFERIEURES
!      I ENROBS        ENROBAGE DES ARMATURES SUPERIEURES
!      I SIGS          CONTRAINTE ADMISSIBLE DANS L'ACIER
!      I SIGCI         CONTRAINTE DE COMPRESSION MAXI DU BETON EN FACE INF
!      I SIGCS         CONTRAINTE DE COMPRESSION MAXI DU BETON EN FACE SUP
!      I ALPHACC       COEFFICIENT DE SECURITE SUR LA RESISTANCE
!                      DE CALCUL DU BETON EN COMPRESSION
!      I GAMMAS        COEFFICIENT DE SECURITE SUR LA RESISTANCE
!                      DE CALCUL DES ACIERS
!      I GAMMAC        COEFFICIENT DE SECURITE SUR LA RESISTANCE
!                      DE CALCUL DU BETON
!      I FACIER        LIMITE D'ELASTICITE DES ACIERS (CONTRAINTE)
!      I EYS           MODULE D'YOUNG DE L'ACIER
!      I TYPDIAG       TYPE DE DIAGRAMME UTILISÉ POUR L'ACIER
!                            TYPDIAG = 1 ("B1" ==> PALIER INCLINÉ)
!                            TYPDIAG = 2 ("B2" ==> PALIER HORIZONTAL)
!      I FBETON        RESISTANCE EN COMPRESSION DU BETON (CONTRAINTE)
!      I CLACIER       CLASSE DE DUCTILITE DES ACIERS (UTILISE POUR EC2) :
!                            CLACIER = 0 ACIER PEU DUCTILE (CLASSE A)
!                            CLACIER = 1 ACIER MOYENNEMENT DUCTILE (CLASSE B)
!                            CLACIER = 2 ACIER FORTEMENT DUCTILE (CLASSE C)
!      I UC            UNITE DES CONTRAINTES :
!                            UC = 0 CONTRAINTES EN Pa
!                            UC = 1 CONTRAINTES EN MPa
!      I UM            UNITE DES DIMENSIONS :
!                            UM = 0 DIMENSIONS EN m
!                            UM = 1 DIMENSIONS EN mm
!      I HT            EPAISSEUR DE LA COQUE
!      DNSXI           DENSITE DE FERRAILLAGE SELON X FACE INF
!      DNSXI           DENSITE DE FERRAILLAGE SELON X FACE SUP
!      DNSXI           DENSITE DE FERRAILLAGE SELON Y FACE INF
!      DNSXI           DENSITE DE FERRAILLAGE SELON Y FACE SUP
!      I EFFRTS        (DIM 8) TORSEUR DES EFFORTS D ETUDE, MOMENTS, ...
!      I EFFREF        (DIM 8) TORSEUR DES EFFORTS DE REFERENCE, MOMENTS, ...
!      O MARGE         LA MARGE MECANIQUE de la facette critique
!      O Nrd,Mrd       Efforts N,M : projection de la droite des efforts sur le diagramme NM
!      O tau_crit      Taux d utilisation sur la facette critique : 1-marge
!      O c0c_crit      Distance C0C entre le point de ref. et le point d etude
!      O c0crd_crit    Distance C0CRD entre le point de ref. et le point du diagramme
!      O theta_crit    L'ANGLE CORRESPONDANT A LA FACETTE CRITIQUE
!      O dnsinf_crit   DENSITE DE FERRAILLAGE FACE INF SUR LA FACETTE CRITIQUE
!      O dnssup_crit   DENSITE DE FERRAILLAGE FACE INF SUR LA FACETTE CRITIQUE
!      O effn0_crit    LES EFFORTS DE REFERENCE (N) PROJETES SUR LA FACETTE CRITIQUE
!      O effm0_crit    LES EFFORTS DE REFERENCE (M) PROJETES SUR LA FACETTE CRITIQUE
!      O effn_crit     LES EFFORTS DE CALCUL (N) PROJETES SUR LA FACETTE CRITIQUE
!      O effm_crit     LES EFFORTS DE CALCUL (M) PROJETES SUR LA FACETTE CRITIQUE
!      O myNrd_crit    N DU POINT D INTERSECTION DE LA DROITE DES EFFORTS AVEC LE
!                      DIAGRAMME D INTERACTION
!      O myMrd_crit    M DU POINT D INTERSECTION DE LA DROITE DES EFFORTS AVEC LE
!                      DIAGRAMME D INTERACTION
!_____________________________________________________________________

      implicit none
!
!
#include "asterfort/cmarge.h"
#include "asterfort/jedetr.h"
#include "asterfort/trgfct.h"
#include "asterfort/wkvect.h"
#include "extern/dintels.h"
#include "extern/dintelu.h"
#include "asterc/r8prem.h"
!
!
      integer(kind=8) :: typcmb
      integer(kind=8) :: typco
      integer(kind=8) :: nb
      real(kind=8) :: cequi
      real(kind=8) :: enrobi
      real(kind=8) :: enrobs
      real(kind=8) :: sigs
      real(kind=8) :: sigci
      real(kind=8) :: sigcs
      real(kind=8) :: alphacc
      real(kind=8) :: gammas
      real(kind=8) :: gammac
      real(kind=8) :: facier
      real(kind=8) :: eys
      integer(kind=8) :: typdiag
      real(kind=8) :: fbeton
      integer(kind=8) :: clacier
      integer(kind=8) :: uc
      real(kind=8) :: ht
      real(kind=8) :: effrts(8), effref(8)
      real(kind=8) :: dnsxi, dnsxs, dnsyi, dnsys
      real(kind=8) :: marge
      real(kind=8), allocatable :: margefct(:)
      real(kind=8) :: bw
      real(kind=8) :: dnsinf
      real(kind=8) :: dnssup

! !
!       NOMBRE DE DIVISIONS ENTRE -PI/2 ET +PI/2
      real(kind=8), allocatable :: fcttab(:, :)
!       EFFORT NORMAL DANS CETTE DIRECTION
      real(kind=8) :: effn, effn0
!       MOMENT DE FLEXION DANS CETTE DIRECTION
      real(kind=8) :: effm, effm0
!       VARIABLE D'ITERATION
      integer(kind=8) :: i
!       AUTRES VARIABLES
      character(24) :: pnrd, pmrd
      integer(kind=8) :: icrit
      real(kind=8), pointer :: mrd(:) => null()
      real(kind=8), pointer :: nrd(:) => null()
      logical :: is_facette_crit
      real(kind=8) :: myNrd, myMrd, taux_util, tau_crit, c0c, c0_crd, c0c_crit, c0crd_crit
      real(kind=8) :: dnsinf_crit, dnssup_crit, myNrd_crit, myMrd_crit, theta
      real(kind=8) :: theta_crit, effn_crit, effm_crit, effn0_crit, effm0_crit
      integer(kind=8) :: s, ntot, ndemi

      pnrd = 'POINT_NRD'
      pmrd = 'POINT_MRD'

!Initialisation
      ntot = 1
      ndemi = 1
      s = 1
!
!   INITIALISATION DES FACETTES
      allocate (fcttab(nb, 6))
      call trgfct(nb, fcttab)
!   Largeur de la section
      bw = 1.0

!   BOUCLE SUR LES FACETTES DE CAPRA ET MAURY
!   DETERMINATION DE LA MARGE POUR CHACUNE DES FACETTES
      allocate (margefct(nb))
      do 10 i = 1, nb

          theta = fcttab(i, 6)

          effn = fcttab(i, 1)*effrts(1)+fcttab(i, 2)*effrts(2)+fcttab(i, 3)*effrts(3)
          effm = fcttab(i, 1)*effrts(4)+fcttab(i, 2)*effrts(5)+fcttab(i, 3)*effrts(6)
          effn = -effn
          effm = -effm
          dnsinf = fcttab(i, 1)*dnsxi+fcttab(i, 2)*dnsyi
          dnssup = fcttab(i, 1)*dnsxs+fcttab(i, 2)*dnsys
!efforts de reference
          effn0 = fcttab(i, 1)*effref(1)+fcttab(i, 2)*effref(2)+fcttab(i, 3)*effref(3)
          effm0 = fcttab(i, 1)*effref(4)+fcttab(i, 2)*effref(5)+fcttab(i, 3)*effref(6)
          effn0 = -effn0
          effm0 = -effm0

!       CALCUL DE LA MARGE A L'ELU

          if (typcmb .eq. 0) then

!           Diagramme d'interaction A L'ELU
!

! call wkvect(pnrd, ' V V R ', ntot, vr=nrd)
! call wkvect(pmrd, ' V V R ', ntot, vr=mrd)
! call dintelu(typco, alphacc, ht, bw, enrobi, enrobs, facier, fbeton, &
!     gammas, gammac, clacier, eys, typdiag, uc, &
!     dnsinf, dnssup, ntot, nrd, mrd)

! compute vectors sizes
              ntot = -1
              call dintelu(typco, alphacc, ht, bw, enrobi, enrobs, facier, fbeton, &
                           gammas, gammac, clacier, eys, typdiag, uc, &
                           ntot, ndemi=ndemi)
              call wkvect(pnrd, ' V V R ', ntot, vr=nrd)
              call wkvect(pmrd, ' V V R ', ntot, vr=mrd)

              do s = 1, ntot
                  nrd(s) = -1.0
                  mrd(s) = -1.0
              end do
! compute nrd & mrd
              call dintelu(typco, alphacc, ht, bw, enrobi, enrobs, facier, fbeton, &
                           gammas, gammac, clacier, eys, typdiag, uc, &
                           ntot, dnsinf, dnssup, nrd, mrd)

              call cmarge(nrd, mrd, effn, effm, effn0, effm0, margefct(i), &
                          myNrd, myMrd, taux_util, c0c, c0_crd)
              if (abs(margefct(i)-2.d0) .le. r8prem()*max(abs(margefct(i)), 2.d0)) then
                  is_facette_crit = .TRUE.
                  icrit = i
                  marge = margefct(icrit)
                  theta_crit = theta
                  effn0_crit = effn0
                  effm0_crit = effm0
                  effn_crit = effn
                  effm_crit = effm
                  dnsinf_crit = dnsinf
                  dnssup_crit = dnssup
                  myNrd_crit = myNrd
                  myMrd_crit = myMrd
                  tau_crit = taux_util
                  c0c_crit = c0c
                  c0crd_crit = c0_crd
                  call jedetr(pnrd)
                  call jedetr(pmrd)
                  goto 111
              end if

              if (i == 1) then
                  is_facette_crit = .TRUE.
                  icrit = i
                  marge = margefct(icrit)
                  theta_crit = theta
                  effn0_crit = effn0
                  effm0_crit = effm0
                  effn_crit = effn
                  effm_crit = effm
                  dnsinf_crit = dnsinf
                  dnssup_crit = dnssup
                  myNrd_crit = myNrd
                  myMrd_crit = myMrd
                  tau_crit = taux_util
                  c0c_crit = c0c
                  c0crd_crit = c0_crd
              else if (i > 1) then
                  if (margefct(i) < margefct(icrit)) then
                      is_facette_crit = .TRUE.
                      icrit = i
                      marge = margefct(icrit)
                      theta_crit = theta
                      effn0_crit = effn0
                      effm0_crit = effm0
                      effn_crit = effn
                      effm_crit = effm
                      dnsinf_crit = dnsinf
                      dnssup_crit = dnssup
                      myNrd_crit = myNrd
                      myMrd_crit = myMrd
                      tau_crit = taux_util
                      c0c_crit = c0c
                      c0crd_crit = c0_crd
                  end if
              end if

              call jedetr(pnrd)
              call jedetr(pmrd)

!       CALCUL DE LA MARGE A L'ELS
          elseif (typcmb .eq. 1) then

!           Diagramme d'interaction A L'ELS
!

! compute vectors sizes
              ntot = -1
              call dintels(cequi, ht, bw, enrobi, enrobs, &
                           sigci, sigcs, sigs, uc, &
                           ntot, ndemi=ndemi)

              call wkvect(pnrd, ' V V R ', ntot, vr=nrd)
              call wkvect(pmrd, ' V V R ', ntot, vr=mrd)

              do s = 1, ntot
                  nrd(s) = -1.0
                  mrd(s) = -1.0
              end do
! compute nrd & mrd
              call dintels(cequi, ht, bw, enrobi, enrobs, &
                           sigci, sigcs, sigs, uc, &
                           ntot, dnsinf, dnssup, nrd, mrd)
              call cmarge(nrd, mrd, effn, effm, effn0, effm0, margefct(i), myNrd, &
                          myMrd, taux_util, c0c, c0_crd)

              if (abs(margefct(i)-2.d0) .le. r8prem()*max(abs(margefct(i)), 2.d0)) then
                  is_facette_crit = .TRUE.
                  icrit = i
                  marge = margefct(icrit)
                  theta_crit = theta
                  effn0_crit = effn0
                  effm0_crit = effm0
                  effn_crit = effn
                  effm_crit = effm
                  dnsinf_crit = dnsinf
                  dnssup_crit = dnssup
                  myNrd_crit = myNrd
                  myMrd_crit = myMrd
                  tau_crit = taux_util
                  c0c_crit = c0c
                  c0crd_crit = c0_crd
                  call jedetr(pnrd)
                  call jedetr(pmrd)
                  goto 111
              end if

              if (i == 1) then
                  is_facette_crit = .TRUE.
                  icrit = i
                  marge = margefct(icrit)
                  theta_crit = theta
                  effn0_crit = effn0
                  effm0_crit = effm0
                  effn_crit = effn
                  effm_crit = effm
                  dnsinf_crit = dnsinf
                  dnssup_crit = dnssup
                  myNrd_crit = myNrd
                  myMrd_crit = myMrd
                  tau_crit = taux_util
                  c0c_crit = c0c
                  c0crd_crit = c0_crd
              else if (i > 1) then
                  if (margefct(i) < margefct(icrit)) then
                      is_facette_crit = .TRUE.
                      icrit = i
                      marge = margefct(icrit)
                      theta_crit = theta
                      effn0_crit = effn0
                      effm0_crit = effm0
                      effn_crit = effn
                      effm_crit = effm
                      dnsinf_crit = dnsinf
                      dnssup_crit = dnssup
                      myNrd_crit = myNrd
                      myMrd_crit = myMrd
                      tau_crit = taux_util
                      c0c_crit = c0c
                      c0crd_crit = c0_crd
                  end if
              end if

              call jedetr(pnrd)
              call jedetr(pmrd)

          end if

10        continue

111       continue
          deallocate (fcttab)
          deallocate (margefct)

          end subroutine
