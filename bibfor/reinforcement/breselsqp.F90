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

subroutine breselsqp(cequi, effmy, effmz, effn, ht, bw, &
                     enrobyi, enrobys, enrobzi, enrobzs, &
                     wmaxyi, wmaxys, wmaxzi, wmaxzs, &
                     ferrcomp, precs, ferrsyme, slsyme, uc, um, &
                     kt, eys, facier, fbeton, sigelsqp, &
                     phiyi, phiys, phizi, phizs, &
                     dnsyi, dnsys, dnszi, dnszs, &
                     sigmsyi, sigmsys, sigmcyi, sigmcys, &
                     sigmszi, sigmszs, sigmczi, sigmczs, &
                     alphay, alphaz, pivoty, pivotz, etaty, etatz, &
                     wfinyi, wfinys, wfinzi, wfinzs, kvarfy, kvarfz, ierr)

!______________________________________________________________________
!
!      BRESELSQP

!      CALCUL DES ACIERS EN FLEXION COMPOSEE DEVIEE A L'ELS QP
!      METHODE DE VERIFICATION DE BRESLER
!      CRITERE = LIMITATION DES OUVERTURES DES FISSURES
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
!      I WMAXYI    OUVERTURE MAX DE FISS ADMIS EN FIBRE INF SUIVANT L'AXE Y
!      I WMAXYS    OUVERTURE MAX DE FISS ADMIS EN FIBRE SUP SUIVANT L'AXE Y
!      I WMAXZI    OUVERTURE MAX DE FISS ADMIS EN FIBRE INF SUIVANT L'AXE Z
!      I WMAXZS    OUVERTURE MAX DE FISS ADMIS EN FIBRE SUP SUIVANT L'AXE Z
!      I FERRCOMP  PRISE EN COMPTE DU FERRAILLAGE DE COMPRESSION
!                     FERRCOMP = 1 (NON)
!                     FERRCOMP = 2 (OUI)
!      I PRECS     PRECISION ITERATION
!      I FERRSYME   FERRAILLAGE SYMETRIQUE?
!                     FERRSYME = 0 (NON)
!                     FERRSYME = 1 (OUI)
!      I SLSYME    SECTION SEUIL DE TOLERANCE POUR UN FERRAILLAGE SYMETRIQUE
!      I UC        UNITE DES CONTRAINTES :
!                     UC = 0 CONTRAINTES EN Pa
!                     UC = 1 CONTRAINTES EN MPa
!      I UM        UNITE DES DIMENSIONS :
!                     UM = 0 DIMENSIONS EN m
!                     UM = 1 DIMENSIONS EN mm
!      I KT        COEFFICIENT DE DUREE DE CHARGEMENT
!      I FACIER    LIMITE D'ELASTICITÉ DES ACIERS (CONTRAINTE)
!      I FBETON    RESISTANCE EN COMPRESSION DU BÉTON (CONTRAINTE)
!      I SIGELSQP  CONTRAINTE ADMISSIBLE DANS LE BETON À L'ELS QP
!      I PHIYINF   DIAMETRE ESTIMATIF DES ARMA EN FIBRE INF SUIVANT L'AXE Y
!      I PHIYSUP   DIAMETRE ESTIMATIF DES ARMA EN FIBRE SUP SUIVANT L'AXE Y
!      I PHIZINF   DIAMETRE ESTIMATIF DES ARMA EN FIBRE INF SUIVANT L'AXE Z
!      I PHIZSUP   DIAMETRE ESTIMATIF DES ARMA EN FIBRE SUP SUIVANT L'AXE Z
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
!      O WFINYI    OUVERTURE DES FISSURES EN FIBRE INF SUIVANT L'AXE Y
!      O WFINYS    OUVERTURE DES FISSURES EN FIBRE SUP SUIVANT L'AXE Y
!      O WFINZI    OUVERTURE DES FISSURES EN FIBRE INF SUIVANT L'AXE Z
!      O WFINZS    OUVERTURE DES FISSURES EN FIBRE SUP SUIVANT L'AXE Z
!      O KVARFY    TAUX DE CONTRT DE TRACT LIM DS L'ACIER SUIVANT L'AXE Y
!      O KVARFZ    TAUX DE CONTRT DE TRACT LIM DS L'ACIER SUIVANT L'AXE Z
!      O IERR      CODE RETOUR (0 = OK)
!
!_______________________________________________________________________
!
!
    implicit none
!
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/jedetr.h"
#include "asterfort/juveca.h"
#include "asterfort/jeveuo.h"
#include "asterfort/cafelsqp.h"
#include "asterfort/dintels.h"
#include "asterc/r8prem.h"
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
    real(kind=8) :: wmaxyi
    real(kind=8) :: wmaxys
    real(kind=8) :: wmaxzi
    real(kind=8) :: wmaxzs
    integer :: ferrcomp
    integer :: precs
    integer :: ferrsyme
    real(kind=8) :: slsyme
    integer :: uc
    integer :: um
    real(kind=8) :: kt
    real(kind=8) :: eys
    real(kind=8) :: facier
    real(kind=8) :: fbeton
    real(kind=8) :: sigelsqp
    real(kind=8) :: phiyi
    real(kind=8) :: phiys
    real(kind=8) :: phizi
    real(kind=8) :: phizs
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
    real(kind=8) :: wfinyi
    real(kind=8) :: wfinys
    real(kind=8) :: wfinzi
    real(kind=8) :: wfinzs
    real(kind=8) :: kvarfy
    real(kind=8) :: kvarfz
    integer :: ierr

!-----------------------------------------------------------------------
!!!!VARIABLES DE CALCUL
!-----------------------------------------------------------------------
    real(kind=8) :: Acc, fcd, fyd, coeff, Ass, Aiter
    real(kind=8) :: rhoyi, rhoys, rhozi, rhozs
    real(kind=8) :: BRES, mrdyE, mrdy1, mrdy2, mrdzE, mrdz1, mrdz2, nrdyzE, a, nrd0
    logical :: COND
    integer :: s, COUNT_BRES
    real(kind=8), pointer :: nrdy(:) => null(), mrdy(:) => null()
    real(kind=8), pointer :: nrdz(:) => null(), mrdz(:) => null()
    character(24) :: pnrdy, pmrdy, pnrdz, pmrdz
    real(kind=8) :: unite_pa, unite_m, Calc, seuil_moment
    real(kind=8) :: d, d0, dneg, d0neg, scmax, scmaxneg, ssmaxy, ssmaxz
    integer :: N_ET, N_PC, N_PCAC, N_EC, N_ECN, N_PCACN
    integer :: ntoty, ndemiy, ntotz, ndemiz

    pnrdy = 'POINT_NRD_Y'
    pmrdy = 'POINT_MRD_Y'
    pnrdz = 'POINT_NRD_Z'
    pmrdz = 'POINT_MRD_Z'

    Acc = bw*ht
    fcd = sigelsqp

    !Initialisation
    kvarfy = 1.0
    kvarfz = 1.0
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
    seuil_moment = sqrt(r8prem())

    !Effort Axial uniquement
    !if ((effmy.eq.0) .and. (effmz.eq.0) .and. (effn.ne.0)) then
    if ((abs(effmy) .lt. seuil_moment) .and. (abs(effmz) .lt. seuil_moment)) then
        call cafelsqp(cequi, effmy, 0.5*effn, ht, bw, &
                      enrobzi, enrobzs, wmaxzi, wmaxzs, &
                      ferrcomp, precs, ferrsyme, slsyme, uc, um, &
                      kt, facier, fbeton, eys, sigelsqp, phizi, phizs, &
                      dnszi, dnszs, sigmszi, sigmszs, sigmczi, sigmczs, &
                      alphaz, pivotz, etatz, &
                      wfinzi, wfinzs, kvarfz, ierr)
        if (ierr .ne. 0) then
            goto 998
        end if
        call cafelsqp(cequi, effmz, 0.5*effn, bw, ht, &
                      enrobyi, enrobys, wmaxyi, wmaxys, &
                      ferrcomp, precs, ferrsyme, slsyme, uc, um, &
                      kt, facier, fbeton, eys, sigelsqp, phiyi, phiys, &
                      dnsyi, dnsys, sigmsyi, sigmsys, sigmcyi, sigmcys, &
                      alphay, pivoty, etaty, &
                      wfinyi, wfinys, kvarfy, ierr)
        if (ierr .ne. 0) then
            goto 998
        end if

    else

        !Calcul suivant "y"
        !if (effmy.ne.0) then
        if (abs(effmy) .gt. seuil_moment) then
            call cafelsqp(cequi, effmy, effn, ht, bw, &
                          enrobzi, enrobzs, wmaxzi, wmaxzs, &
                          ferrcomp, precs, ferrsyme, slsyme, uc, um, &
                          kt, facier, fbeton, eys, sigelsqp, phizi, phizs, &
                          dnszi, dnszs, sigmszi, sigmszs, sigmczi, sigmczs, &
                          alphaz, pivotz, etatz, &
                          wfinzi, wfinzs, kvarfz, ierr)
            if (ierr .ne. 0) then
                goto 998
            end if
        end if

        !Calcul suivant "z"
        !if (effmz.ne.0) then
        if (abs(effmz) .gt. seuil_moment) then
            call cafelsqp(cequi, effmz, effn, bw, ht, &
                          enrobyi, enrobys, wmaxyi, wmaxys, &
                          ferrcomp, precs, ferrsyme, slsyme, uc, um, &
                          kt, facier, fbeton, eys, sigelsqp, phiyi, phiys, &
                          dnsyi, dnsys, sigmsyi, sigmsys, sigmcyi, sigmcys, &
                          alphay, pivoty, etaty, &
                          wfinyi, wfinys, kvarfy, ierr)
            if (ierr .ne. 0) then
                goto 998
            end if
        end if

    end if

    fyd = 0.5*(kvarfy+kvarfz)*facier

    !if ((effmy.ne.0) .and. (effmz.ne.0)) then
    if ((abs(effmy) .gt. seuil_moment) .and. (abs(effmz) .gt. seuil_moment)) then

        !Iteration Bresler
        COND = .false.
        COUNT_BRES = 0
        BRES = 1.5

        !Dimensionnement des vecteurs

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

        scmax = sigelsqp
        scmaxneg = sigelsqp

        ssmaxy = kvarfy*facier
        ssmaxz = kvarfz*facier

        !Pour MFY
        d = ht-enrobzi
        d0 = enrobzs
        dneg = ht-enrobzs
        d0neg = enrobzi

        N_ET = 11
        N_PC = 101
        N_EC = CEILING(10*(scmax*unite_pa))+1
        N_ECN = CEILING(10*(scmaxneg*unite_pa))+1

        N_PCAC = CEILING((N_PC-1)*(ht/d))+1
        N_PCACN = CEILING((N_PC-1)*(ht/dneg))+1
        ntoty = N_ET+N_PCAC+N_EC+N_ECN+N_PCACN+N_ET
        ndemiy = N_ET+N_PCac+N_EC
        call wkvect(pnrdy, ' V V R ', ntoty, vr=nrdy)
        call wkvect(pmrdy, ' V V R ', ntoty, vr=mrdy)

        !Pour MFZ
        d = bw-enrobyi
        d0 = enrobys
        dneg = bw-enrobys
        d0neg = enrobyi

        N_PCAC = CEILING((N_PC-1)*(bw/d))+1
        N_PCACN = CEILING((N_PC-1)*(bw/dneg))+1
        ntotz = N_ET+N_PCAC+N_EC+N_ECN+N_PCACN+N_ET
        ndemiz = N_ET+N_PCac+N_EC
        call wkvect(pnrdz, ' V V R ', ntotz, vr=nrdz)
        call wkvect(pmrdz, ' V V R ', ntotz, vr=mrdz)

        do while (COND .eqv. (.false.))

            Ass = dnsyi+dnsys+dnszi+dnszs
            nrdyzE = Acc*fcd+Ass*fyd

            !Determiner MRd,y

            do s = 1, ntoty
                nrdy(s) = -1.0
                mrdy(s) = -1.0
            end do

            call dintels(cequi, ht, bw, enrobzi, enrobzs, &
                         sigelsqp, sigelsqp, ssmaxz, uc, &
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
            if (effmy .gt. 0.0) then
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

            call dintels(cequi, bw, ht, enrobyi, enrobys, &
                         sigelsqp, sigelsqp, ssmaxy, uc, &
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
            if (effmz .gt. 0.0) then
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
                    rhoyi = 0.25
                    rhoys = 0.25
                    rhozi = 0.25
                    rhozs = 0.25
                else
                    if (ferrsyme .eq. 1) then
                        rhoyi = 0.5*(dnsyi+dnsys)/Ass
                        rhoys = 0.5*(dnsyi+dnsys)/Ass
                        rhozi = 0.5*(dnszi+dnszs)/Ass
                        rhozs = 0.5*(dnszi+dnszs)/Ass
                    else
                        rhoyi = dnsyi/Ass
                        rhoys = dnsys/Ass
                        rhozi = dnszi/Ass
                        rhozs = dnszs/Ass
                    end if
                end if
                Aiter = 0.10*Ass
                dnsyi = dnsyi+rhoyi*Aiter
                dnsys = dnsys+rhoys*Aiter
                dnszi = dnszi+rhozi*Aiter
                dnszs = dnszs+rhozs*Aiter
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
