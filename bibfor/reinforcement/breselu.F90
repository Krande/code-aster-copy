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

subroutine breselu(typco, alphacc, effmy, effmz, effn, &
                   ht, bw, enrobyi, enrobys, enrobzi, enrobzs, &
                   facier, fbeton, gammas, gammac, &
                   clacier, eys, typdiag, ferrcomp, precs, ferrsyme, slsyme, &
                   uc, um, &
                   dnsyi, dnsys, dnszi, dnszs, &
                   sigmsyi, sigmsys, ecyi, ecys, &
                   sigmszi, sigmszs, eczi, eczs, &
                   alphay, alphaz, pivoty, pivotz, etaty, etatz, ierr)
!______________________________________________________________________
!
!      BRESELU

!      CALCUL DES ACIERS EN FLEXION COMPOSEE DEVIEE A L'ELU
!      METHODE DE VERIFICATION DE BRESLER
!      CRITERE = LIMITATION DES DEFORMATIONS
!
!      I TYPCO     CODIFICATION UTILISEE (1 = BAEL91, 2 = EC2)
!      I ALPHACC   COEFFICIENT DE SECURITE SUR LA RESISTANCE
!                  DE CALCUL DU BETON EN COMPRESSION
!      I EFFMY     MOMENT DE FLEXION SUIVANT L'AXE Y
!      I EFFMZ     MOMENT DE FLEXION SUIVANT L'AXE Z
!      I EFFN      EFFORT NORMAL
!      I HT        HAUTEUR DE LA SECTION
!      I BW        LARGEUR DE LA SECTION
!      I ENROBYI   ENROBAGE DES ARMATURES INF SUIVANT L'AXE Y
!      I ENROBYS   ENROBAGE DES ARMATURES SUP SUIVANT L'AXE Y
!      I ENROBZI   ENROBAGE DES ARMATURES INF SUIVANT L'AXE Z
!      I ENROBZS   ENROBAGE DES ARMATURES SUP SUIVANT L'AXE Z
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
!                     FERRCOMP = 1 (NON)
!                     FERRCOMP = 2 (OUI)
!      I PRECS     PRECISION SUPPLEMENTAIRE DANS LA RECHERCHE DE L'OPTIMUM
!                   POUR LA METHODE DES 3 PIVOTS (Intervention du 03/2023)
!                     PRECS = 0 (NON)
!                     PRECS = 1 (OUI)
!      I FERRSYME  FERRAILLAGE SYMETRIQUE?
!                     FERRSYME = 0 (NON)
!                     FERRSYME = 1 (OUI)
!      I SLSYME    SECTION SEUIL DE TOLERANCE POUR
!                     UN FERRAILLAGE SYMETRIQUE
!      I UC        UNITE DES CONTRAINTES :
!                     UC = 0 CONTRAINTES EN Pa
!                     UC = 1 CONTRAINTES EN MPa
!      I UM        UNITE DES DIMENSIONS :
!                     UM = 0 DIMENSIONS EN m
!                     UM = 1 DIMENSIONS EN mm
!
!      O DNSYI     DENSITE DE L'ACIER INF SUIVANT L'AXE Y
!      O DNSYS     DENSITE DE L'ACIER SUP SUIVANT L'AXE Y
!      O DNSZI     DENSITE DE L'ACIER INF SUIVANT L'AXE Z
!      O DNSZS     DENSITE DE L'ACIER SUP SUIVANT L'AXE Z
!      O SIGMSYI   CONTRAINTE DANS L'ACIER INF SUIVANT L'AXE Y
!      O SIGMSYS   CONTRAINTE DANS L'ACIER SUP SUIVANT L'AXE Y
!      O ECYI      DEFORMATION EN FIBRE INF SUIVANT L'AXE Y
!      O ECYS      DEFORMATION EN FIBRE SUP SUIVANT L'AXE Y
!      O SIGMSZI   CONTRAINTE DANS L'ACIER INF SUIVANT L'AXE Z
!      O SIGMSZS   CONTRAINTE DANS L'ACIER SUP SUIVANT L'AXE Z
!      O ECZI      DEFORMATION EN FIBRE INF SUIVANT L'AXE Z
!      O ECZS      DEFORMATION EN FIBRE SUP SUIVANT L'AXE Z
!      O ALPHAY    COEF DE PROFONDEUR DE L'AN SUIVANT L'AXE Y
!      O ALPHAZ    COEF DE PROFONDEUR DE L'AN SUIVANT L'AXE Z
!      O PIVOTY    PIVOT EN FC SUIVANT L'AXE Y
!      O PIVOTZ    PIVOT EN FC SUIVANT L'AXE Z
!      O ETATY     ETAT EN FC SUIVANT L'AXE Y
!      O ETATZ     ETAT EN FC SUIVANT L'AXE Z
!      O IERR      CODE RETOUR (0 = OK)
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
#include "asterfort/cafelu.h"
#include "asterfort/dintelu.h"
!
    integer :: typco
    real(kind=8) :: alphacc
    real(kind=8) :: effmy
    real(kind=8) :: effmz
    real(kind=8) :: effn
    real(kind=8) :: ht
    real(kind=8) :: bw
    real(kind=8) :: enrobyi
    real(kind=8) :: enrobys
    real(kind=8) :: enrobzi
    real(kind=8) :: enrobzs
    real(kind=8) :: facier
    real(kind=8) :: fbeton
    real(kind=8) :: gammas
    real(kind=8) :: gammac
    integer :: clacier
    real(kind=8) :: eys
    integer :: typdiag
    integer :: ferrcomp
    integer :: precs
    integer :: ferrsyme
    real(kind=8) :: slsyme
    integer :: uc
    integer :: um
    real(kind=8) :: dnsyi
    real(kind=8) :: dnsys
    real(kind=8) :: dnszi
    real(kind=8) :: dnszs
    real(kind=8) :: sigmsyi
    real(kind=8) :: sigmsys
    real(kind=8) :: ecyi
    real(kind=8) :: ecys
    real(kind=8) :: sigmszi
    real(kind=8) :: sigmszs
    real(kind=8) :: eczi
    real(kind=8) :: eczs
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
    real(kind=8) :: Acc, fcd, fyd, coeff, Ass, Aiter, Calc
    real(kind=8) :: rhoyinf, rhoysup, rhozinf, rhozsup
    real(kind=8) :: BRES, mrdyE, mrdy1, mrdy2, mrdzE, mrdz1, mrdz2, nrdyzE, a, nrd0
    logical :: COND
    integer :: s, COUNT_BRES
    real(kind=8), pointer :: nrdy(:) => null(), mrdy(:) => null()
    real(kind=8), pointer :: nrdz(:) => null(), mrdz(:) => null()
    character(24) :: pnrdy, pmrdy, pnrdz, pmrdz
    real(kind=8) :: unite_pa, unite_m
    real(kind=8) :: piv_a, piv_b, piv_c, Esu, d, d0, dneg, d0neg, Xsup
    integer :: N_ET, N_PC, N_PCN, N_EC, ntoty, ndemiy, ntotz, ndemiz

    pnrdy = 'POINT_NRD_Y'
    pmrdy = 'POINT_MRD_Y'
    pnrdz = 'POINT_NRD_Z'
    pmrdz = 'POINT_MRD_Z'

    Acc = bw*ht
    fcd = fbeton/gammac
    fyd = facier/gammas

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
    if ((abs(effmy) .lt. epsilon(effmy)) .and. (abs(effmz) .lt. epsilon(effmz))) then
        call cafelu(typco, alphacc, effmy, 0.5*effn, ht, bw, &
                    enrobzi, enrobzs, facier, fbeton, gammas, gammac, &
                    clacier, eys, typdiag, ferrcomp, precs, ferrsyme, slsyme, uc, um, &
                    dnszi, dnszs, sigmszi, sigmszs, eczi, eczs, &
                    alphaz, pivotz, etatz, ierr)
        if (ierr .ne. 0) then
            goto 998
        end if
        call cafelu(typco, alphacc, effmz, 0.5*effn, bw, ht, &
                    enrobyi, enrobys, facier, fbeton, gammas, gammac, &
                    clacier, eys, typdiag, ferrcomp, precs, ferrsyme, slsyme, uc, um, &
                    dnsyi, dnsys, sigmsyi, sigmsys, ecyi, ecys, &
                    alphay, pivoty, etaty, ierr)
        if (ierr .ne. 0) then
            goto 998
        end if

    else

        !Calcul suivant "y"
        !if (effmy.ne.0) then
        if (abs(effmy) .gt. epsilon(effmy)) then
            call cafelu(typco, alphacc, effmy, effn, ht, bw, &
                        enrobzi, enrobzs, facier, fbeton, gammas, gammac, &
                        clacier, eys, typdiag, ferrcomp, precs, ferrsyme, slsyme, uc, um, &
                        dnszi, dnszs, sigmszi, sigmszs, eczi, eczs, &
                        alphaz, pivotz, etatz, ierr)
            if (ierr .ne. 0) then
                goto 998
            end if
        end if

        !Calcul suivant "z"
        !if (effmz.ne.0) then
        if (abs(effmz) .gt. epsilon(effmz)) then
            call cafelu(typco, alphacc, effmz, effn, bw, ht, &
                        enrobyi, enrobys, facier, fbeton, gammas, gammac, &
                        clacier, eys, typdiag, ferrcomp, precs, ferrsyme, slsyme, uc, um, &
                        dnsyi, dnsys, sigmsyi, sigmsys, ecyi, ecys, &
                        alphay, pivoty, etaty, ierr)
            if (ierr .ne. 0) then
                goto 998
            end if
        end if

    end if

    !if ((effmy.ne.0) .and. (effmz.ne.0)) then
    if ((abs(effmy) .gt. epsilon(effmy)) .and. (abs(effmz) .gt. epsilon(effmz))) then

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

        if (clacier .eq. 0) then
            piv_a = 0.9*2.5e-2
        elseif (clacier .eq. 1) then
            piv_a = 0.9*5.e-2
        else
            piv_a = 0.9*7.5e-2
        end if

        Esu = piv_a
        piv_b = min(3.5E-3, 0.26*0.01+3.5*0.01*(((90.d0-fbeton*unite_pa)/100.d0)**4))
        piv_c = 2.0E-3
        if ((fbeton*unite_pa) .ge. (50.d0)) then
            piv_c = 0.2*0.01+0.0085*0.01*((fbeton*unite_pa-50.d0)**(0.53))
        end if
        Xsup = piv_b/piv_c

        !Iteration Bresler

        COND = .false.
        COUNT_BRES = 0
        BRES = 1.5

        !Dimensionnement des vecteurs

        N_ET = floor(Esu*1000)+1
        N_EC = ceiling(Xsup*100)+1

        !Pour MFY
        d = ht-enrobzi
        d0 = enrobzs
        dneg = ht-enrobzs
        d0neg = enrobzi

        N_PC = ceiling((ht/d)*100)+1
        N_PCN = ceiling((ht/dneg)*100)+1
        ntoty = N_ET+N_PC+N_EC+N_EC+N_PCN+N_ET
        ndemiy = N_ET+N_PC+N_EC

        call wkvect(pnrdy, ' V V R ', ntoty, vr=nrdy)
        call wkvect(pmrdy, ' V V R ', ntoty, vr=mrdy)

        !Pour MFZ
        d = bw-enrobyi
        d0 = enrobys
        dneg = bw-enrobys
        d0neg = enrobyi

        N_PC = ceiling((bw/d)*100)+1
        N_PCN = ceiling((bw/dneg)*100)+1
        ntotz = N_ET+N_PC+N_EC+N_EC+N_PCN+N_ET
        ndemiz = N_ET+N_PC+N_EC

        call wkvect(pnrdz, ' V V R ', ntotz, vr=nrdz)
        call wkvect(pmrdz, ' V V R ', ntotz, vr=mrdz)

        do while (COND .eqv. (.false.))

            Ass = dnsyi+dnsys+dnszi+dnszs
            nrdyzE = Acc*fcd+Ass*fyd

            !Determiner MRd,y

            call dintelu(typco, alphacc, ht, bw, enrobzi, enrobzs, facier, fbeton, &
                         gammas, gammac, clacier, eys, typdiag, uc, &
                         dnszi, dnszs, ntoty, nrdy, mrdy)

            s = 1
            nrd0 = nrdy(s)
            do while ((nrd0 .lt. effn) .and. (s .lt. ndemiy))
                s = s+1
                nrd0 = nrdy(s)
            end do
            if ((s .eq. 1) .or. (s .eq. ndemiy)) then
                BRES = 1.5
                goto 999
            else
                Calc = nrdy(s)-nrdy(s-1)
                if (abs(Calc) .gt. epsilon(Calc)) then
                    mrdy1 = ((mrdy(s)-mrdy(s-1))/(nrdy(s)-nrdy(s-1)))*(effn-nrdy(s-1))+mrdy(s-1)
                else
                    mrdy1 = 0.5*(mrdy(s-1)+mrdy(s))
                end if
            end if
            s = ndemiy+1
            nrd0 = nrdy(s)
            do while ((nrd0 .gt. effn) .and. (s .lt. ntoty))
                s = s+1
                nrd0 = nrdy(s)
            end do
            if ((s .eq. (ndemiy+1)) .or. (s .eq. ntoty)) then
                BRES = 1.5
                goto 999
            else
                Calc = nrdy(s)-nrdy(s-1)
                if (abs(Calc) .gt. epsilon(Calc)) then
                    mrdy2 = ((mrdy(s)-mrdy(s-1))/(nrdy(s)-nrdy(s-1)))*(effn-nrdy(s-1))+mrdy(s-1)
                else
                    mrdy2 = 0.5*(mrdy(s-1)+mrdy(s))
                end if
            end if
            if (effmy .ge. 0.0) then
                mrdy1 = max(mrdy1, 0.0)
                mrdy2 = max(mrdy2, 0.0)
                mrdyE = max(mrdy1, mrdy2)
            elseif (effmy .lt. 0.0) then
                mrdy1 = min(mrdy1, 0.0)
                mrdy2 = min(mrdy2, 0.0)
                mrdyE = min(mrdy1, mrdy2)
            end if

            !Determiner MRd,z

            do s = 1, ntotz
                nrdz(s) = -1.0
                mrdz(s) = -1.0
            end do

            call dintelu(typco, alphacc, bw, ht, enrobyi, enrobys, facier, fbeton, &
                         gammas, gammac, clacier, eys, typdiag, uc, &
                         dnsyi, dnsys, ntotz, nrdz, mrdz)

            s = 1
            nrd0 = nrdz(s)
            do while ((nrd0 .lt. effn) .and. (s .lt. ndemiz))
                s = s+1
                nrd0 = nrdz(s)
            end do
            if ((s .eq. 1) .or. (s .eq. ndemiz)) then
                BRES = 1.5
                goto 999
            else
                Calc = nrdz(s)-nrdz(s-1)
                if (abs(Calc) .gt. epsilon(Calc)) then
                    mrdz1 = ((mrdz(s)-mrdz(s-1))/(nrdz(s)-nrdz(s-1)))*(effn-nrdz(s-1))+mrdz(s-1)
                else
                    mrdz1 = 0.5*(mrdz(s-1)+mrdz(s))
                end if
            end if
            s = ndemiz+1
            nrd0 = nrdz(s)
            do while ((nrd0 .gt. effn) .and. (s .lt. ntotz))
                s = s+1
                nrd0 = nrdz(s)
            end do
            if ((s .eq. (ndemiz+1)) .or. (s .eq. ntotz)) then
                BRES = 1.5
                goto 999
            else
                Calc = nrdz(s)-nrdz(s-1)
                if (abs(Calc) .gt. epsilon(Calc)) then
                    mrdz2 = ((mrdz(s)-mrdz(s-1))/(nrdz(s)-nrdz(s-1)))*(effn-nrdz(s-1))+mrdz(s-1)
                else
                    mrdz2 = 0.5*(mrdz(s-1)+mrdz(s))
                end if
            end if
            if (effmz .ge. 0.0) then
                mrdz1 = max(mrdz1, 0.0)
                mrdz2 = max(mrdz2, 0.0)
                mrdzE = max(mrdz1, mrdz2)
            elseif (effmz .lt. 0.0) then
                mrdz1 = min(mrdz1, 0.0)
                mrdz2 = min(mrdz2, 0.0)
                mrdzE = min(mrdz1, mrdz2)
            end if

            !Calcul de 'a'

            if (abs(nrdyzE) .gt. epsilon(nrdyzE)) then
                coeff = effn/nrdyzE
            else
                coeff = 0.0
            end if

            if (coeff .le. 0.1) then
                a = 1.0
            elseif (coeff .le. 0.7) then
                a = ((1.5-1.0)/(0.7-0.1))*(coeff-0.1)+1.0
            elseif (coeff .le. 1.0) then
                a = ((2.0-1.5)/(1.0-0.7))*(coeff-0.7)+1.5
            else
                a = 2.0
            end if

            !Calcul de 'BRES'
            if ((abs(mrdyE) .gt. epsilon(mrdyE)) .and. (abs(mrdzE) .gt. epsilon(mrdzE))) then
                BRES = (effmy/mrdyE)**(a)+(effmz/mrdzE)**(a)
            end if

            !Verif de 'BRES' et iteration
999         continue

            COUNT_BRES = COUNT_BRES+1

            if (BRES .gt. 1) then
                if (Ass .lt. epsilon(Ass)) then
                    Ass = (1.e2)/(unite_m*unite_m)
                    rhoyinf = 0.25
                    rhoysup = 0.25
                    rhozinf = 0.25
                    rhozsup = 0.25
                else
                    if (ferrsyme .eq. 1) then
                        rhoyinf = 0.5*(dnsyi+dnsys)/Ass
                        rhoysup = 0.5*(dnsyi+dnsys)/Ass
                        rhozinf = 0.5*(dnszi+dnszs)/Ass
                        rhozsup = 0.5*(dnszi+dnszs)/Ass
                    else
                        rhoyinf = dnsyi/Ass
                        rhoysup = dnsys/Ass
                        rhozinf = dnszi/Ass
                        rhozsup = dnszs/Ass
                    end if
                end if
                Aiter = 0.10*Ass
                dnsyi = dnsyi+rhoyinf*Aiter
                dnsys = dnsys+rhoysup*Aiter
                dnszi = dnszi+rhozinf*Aiter
                dnszs = dnszs+rhozsup*Aiter
                if (COUNT_BRES .eq. 100) then
                    ierr = 4
                    dnsyi = -1
                    dnsys = -1
                    dnszi = -1
                    dnszs = -1
                    COND = .true.
                end if
            else
                COND = .true.
            end if

        end do
        !do while (COND.eqv.(.false.))

        call jedetr(pnrdy)
        call jedetr(pmrdy)
        call jedetr(pnrdz)
        call jedetr(pmrdz)

    end if
    !if ((effmy.ne.0) .and. (effmz.ne.0)) then

998 continue

end subroutine
