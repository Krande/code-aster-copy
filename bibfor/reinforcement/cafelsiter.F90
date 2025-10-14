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

subroutine cafelsiter(cequi, effm, effn, ht, bw, &
                      enrobi, enrobs, scmaxi, scmaxs, ssmax, ferrcomp, &
                      precs, ferrsyme, slsyme, uc, um, condns, &
                      astend, ascomp, sstend, sscomp, &
                      sctend, sccomp, &
                      alpha, pivot, etat, ierr)

!_______________________________________________________________________________________________
!
!      CAFELSITER
!
!      CALCUL DES ACIERS EN FLEXION COMPOSEE A L'ELS CARACTERISTIQUE PAR ITERATION
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
!                        FERRSYME = 0 (NON)
!                        FERRSYME = 1 (OUI)
!      I SLSYME       SECTION SEUIL DE TOLERANCE POUR UN FERRAILLAGE SYMETRIQUE
!      I UC           UNITE DES CONTRAINTES :
!                        UC = 0 CONTRAINTES EN Pa
!                        UC = 1 CONTRAINTES EN MPa
!      I UM           UNITE DES DIMENSIONS :
!                        UM = 0 DIMENSIONS EN m
!                        UM = 1 DIMENSIONS EN mm
!      I CONDNS       COND_NS DE CAFELS
!
!      O ASTEND       DENSITE DE L'ACIER TENDU
!      O ASCOMP       DENSITE DE L'ACIER COMPRIMÉ
!      O SSTEND       CONTRAINTE AU NIVEAU DE L'ACIER TENDU
!      O SSCOMP       CONTRAINTE AU NIVEAU DE L'ACIER COMPRIMÉ
!      O SCTEND       CONTRAINTE AU NIVEAU DE LA FIBRE DE BETON TENDU
!      O SCCOMP       CONTRAINTE AU NIVEAU DE LA FIBRE DE BETON COMPRIMÉ
!      O ALPHA        COEFFICIENT DE PROFONDEUR DE L'AN
!      O PIVOT        PIVOT DE FONCTIONNEMENT DE LA SECTION
!      O ETAT         ETAT DE FONCTIONNEMENT DE LA SECTION
!      O IERR         CODE RETOUR (0 = OK)
!_______________________________________________________________________________________________
!
    implicit none

#include "asterfort/jedetr.h"
#include "asterfort/jeveuo.h"
#include "asterfort/juveca.h"
#include "asterfort/verifels.h"
#include "asterfort/wkvect.h"
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
    logical :: condns
    real(kind=8) :: astend
    real(kind=8) :: ascomp
    real(kind=8) :: sstend
    real(kind=8) :: sscomp
    real(kind=8) :: sctend
    real(kind=8) :: sccomp
    real(kind=8) :: alpha
    integer(kind=8) :: pivot
    integer(kind=8) :: etat
    integer(kind=8) :: ierr

!-----------------------------------------------------------------------
!!!!VARIABLES DE CALCUL
!-----------------------------------------------------------------------
    real(kind=8) :: d, d0, NUM, DENUM, seuil_as
    real(kind=8) :: Mcalc, Ncalc, scmax, scmaxc, alpha_12, unite_pa, unite_m
    real(kind=8) :: AsTEND_F, AsCOMP_F, AsTOT_F, AsTEND_FOUND, AsCOMP_FOUND, AsTOT_FOUND
    real(kind=8) :: X, Ncc, Mcc, Dterm, Calc, yc
    real(kind=8) :: alphaI, m1, m2, residu, aG, aD, fG, fD
    logical :: COND_AJOUT_RESIDU, COND_AJOUT_AsCOMP, COND_AJOUT_AsTEND, COND_COUNT
    logical :: COND_COUNT_SYME, COND_COUNT_NOCOMP, COND_AJOUT_Proceed
    logical :: COND_COUNT_SYMENOCOMP
    integer(kind=8) :: COUNT_CARA, verif, indx_ch
    integer(kind=8) :: COUNT_CARA_SYME, COUNT_CARA_NOCOMP, COUNT_CARA_SYMENOCOMP
    logical :: COND_F
    integer(kind=8) :: COUNT_F
    integer(kind=8) :: N_ET, N_PC, N_PCAC, N_EC, N_TOT

    real(kind=8), pointer :: SsTEND_ET(:) => null(), SsCOMP_ET(:) => null()
    real(kind=8), pointer :: ScTEND_ET(:) => null(), ScCOMP_ET(:) => null()
    real(kind=8), pointer :: AsTEND_ET(:) => null(), AsCOMP_ET(:) => null()
    real(kind=8), pointer :: alpha_ET(:) => null(), RESIDU_ET(:) => null()
    integer(kind=8), pointer :: PIVOT_ET(:) => null(), ETAT_ET(:) => null()

    real(kind=8), pointer :: SsTEND_PCAC(:) => null(), SsCOMP_PCAC(:) => null()
    real(kind=8), pointer :: ScTEND_PCAC(:) => null(), ScCOMP_PCAC(:) => null()
    real(kind=8), pointer :: AsTEND_PCAC(:) => null(), AsCOMP_PCAC(:) => null()
    real(kind=8), pointer :: alpha_PCAC(:) => null(), RESIDU_PCAC(:) => null()
    integer(kind=8), pointer :: PIVOT_PCAC(:) => null(), ETAT_PCAC(:) => null()

    real(kind=8), pointer :: SsTEND_EC(:) => null(), SsCOMP_EC(:) => null()
    real(kind=8), pointer :: ScTEND_EC(:) => null(), ScCOMP_EC(:) => null()
    real(kind=8), pointer :: AsTEND_EC(:) => null(), AsCOMP_EC(:) => null()
    real(kind=8), pointer :: alpha_EC(:) => null(), RESIDU_EC(:) => null()
    integer(kind=8), pointer :: PIVOT_EC(:) => null(), ETAT_EC(:) => null()

    real(kind=8), pointer :: SsTEND_TOT(:) => null(), SsCOMP_TOT(:) => null()
    real(kind=8), pointer :: ScTEND_TOT(:) => null(), ScCOMP_TOT(:) => null()
    real(kind=8), pointer :: AsTEND_TOT(:) => null(), AsCOMP_TOT(:) => null()
    real(kind=8), pointer :: alpha_TOT(:) => null(), RESIDU_TOT(:) => null()
    integer(kind=8), pointer :: PIVOT_TOT(:) => null(), ETAT_TOT(:) => null()

    real(kind=8), dimension(2, 2) :: Dsys
    real(kind=8), dimension(2) :: SOL
    integer(kind=8) :: i, j, s, INDICE_F

    character(20) :: p(40)

!-----------------------------------------------------------------------
!!!!LANCEMENT DU CALCUL
!-----------------------------------------------------------------------

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

    seuil_as = 0.0

    !!Intervention 05/2023 - Section non ferraillée auto-équilibrée
    call verifels(cequi, ht, bw, enrobi, enrobs, &
                  scmaxi, scmaxs, ssmax, uc, &
                  seuil_as, seuil_as, effm, effn, verif)

    if (verif .eq. 0) then
        !OK
        ascomp = 0
        astend = 0
        sscomp = -1
        sstend = -1
        sccomp = -1
        sctend = -1
        alpha = -1000
        etat = 1
        pivot = 0
        goto 998
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

    Mcalc = abs(effm)
    Ncalc = effn
    N_ET = 11

    N_PC = precs+1
    N_PCAC = CEILING((N_PC-1)*(ht/d))+1

    !Determination Pivot C 'ELS'
    yc = (scmaxc/scmax)*ht
    N_EC = CEILING(10*(scmaxc*unite_pa))+1

    N_TOT = N_ET+N_PCAC+N_EC

    do i = 1, 40
        write (p(i), fmt='(A18,I2)') 'POINT_ITER_CAFELS_', i
    end do

    call wkvect(p(1), ' V V R ', N_ET, vr=SsTEND_ET)
    call wkvect(p(2), ' V V R ', N_ET, vr=SsCOMP_ET)
    call wkvect(p(3), ' V V R ', N_ET, vr=ScTEND_ET)
    call wkvect(p(4), ' V V R ', N_ET, vr=ScCOMP_ET)
    call wkvect(p(5), ' V V R ', N_ET, vr=AsTEND_ET)
    call wkvect(p(6), ' V V R ', N_ET, vr=AsCOMP_ET)
    call wkvect(p(7), ' V V R ', N_ET, vr=alpha_ET)
    call wkvect(p(8), ' V V R ', N_ET, vr=RESIDU_ET)
    call wkvect(p(9), ' V V I ', N_ET, vi=PIVOT_ET)
    call wkvect(p(10), ' V V I ', N_ET, vi=ETAT_ET)

    call wkvect(p(11), ' V V R ', N_PCAC, vr=SsTEND_PCAC)
    call wkvect(p(12), ' V V R ', N_PCAC, vr=SsCOMP_PCAC)
    call wkvect(p(13), ' V V R ', N_PCAC, vr=ScTEND_PCAC)
    call wkvect(p(14), ' V V R ', N_PCAC, vr=ScCOMP_PCAC)
    call wkvect(p(15), ' V V R ', N_PCAC, vr=AsTEND_PCAC)
    call wkvect(p(16), ' V V R ', N_PCAC, vr=AsCOMP_PCAC)
    call wkvect(p(17), ' V V R ', N_PCAC, vr=alpha_PCAC)
    call wkvect(p(18), ' V V R ', N_PCAC, vr=RESIDU_PCAC)
    call wkvect(p(19), ' V V I ', N_PCAC, vi=PIVOT_PCAC)
    call wkvect(p(20), ' V V I ', N_PCAC, vi=ETAT_PCAC)

    call wkvect(p(21), ' V V R ', N_EC, vr=SsTEND_EC)
    call wkvect(p(22), ' V V R ', N_EC, vr=SsCOMP_EC)
    call wkvect(p(23), ' V V R ', N_EC, vr=ScTEND_EC)
    call wkvect(p(24), ' V V R ', N_EC, vr=ScCOMP_EC)
    call wkvect(p(25), ' V V R ', N_EC, vr=AsTEND_EC)
    call wkvect(p(26), ' V V R ', N_EC, vr=AsCOMP_EC)
    call wkvect(p(27), ' V V R ', N_EC, vr=alpha_EC)
    call wkvect(p(28), ' V V R ', N_EC, vr=RESIDU_EC)
    call wkvect(p(29), ' V V I ', N_EC, vi=PIVOT_EC)
    call wkvect(p(30), ' V V I ', N_EC, vi=ETAT_EC)

    call wkvect(p(31), ' V V R ', N_TOT, vr=SsTEND_TOT)
    call wkvect(p(32), ' V V R ', N_TOT, vr=SsCOMP_TOT)
    call wkvect(p(33), ' V V R ', N_TOT, vr=ScTEND_TOT)
    call wkvect(p(34), ' V V R ', N_TOT, vr=ScCOMP_TOT)
    call wkvect(p(35), ' V V R ', N_TOT, vr=AsTEND_TOT)
    call wkvect(p(36), ' V V R ', N_TOT, vr=AsCOMP_TOT)
    call wkvect(p(37), ' V V R ', N_TOT, vr=alpha_TOT)
    call wkvect(p(38), ' V V R ', N_TOT, vr=RESIDU_TOT)
    call wkvect(p(39), ' V V I ', N_TOT, vi=PIVOT_TOT)
    call wkvect(p(40), ' V V I ', N_TOT, vi=ETAT_TOT)

    do i = 1, N_ET
        SsTEND_ET(i) = -ssmax
        ScCOMP_ET(i) = -(ssmax/cequi)*(1-0.1*(i-1))
        SsCOMP_ET(i) = ((SsTEND_ET(i)/cequi-ScCOMP_ET(i))*((ht-d0)/d)+ScCOMP_ET(i))*cequi
        ScTEND_ET(i) = (SsTEND_ET(i)/cequi-ScCOMP_ET(i))*(ht/d)+ScCOMP_ET(i)
        AsTEND_ET(i) = (Mcalc-Ncalc*(d0-0.5*ht))/(SsTEND_ET(i)*(-d-d0+ht))
        AsCOMP_ET(i) = (Ncalc-AsTEND_ET(i)*SsTEND_ET(i))/SsCOMP_ET(i)
        if (i .eq. 1) then
            alpha_ET(i) = -1000
        else
            alpha_ET(i) = -ScCOMP_ET(i)/(SsTEND_ET(i)/cequi-ScCOMP_ET(i))
        end if
        PIVOT_ET(i) = 1
        ETAT_ET(i) = 4
        RESIDU_ET(i) = 0
    end do

    do i = 1, N_PCAC
        if (i .lt. N_PCAC) then
            alpha_PCAC(i) = real(i-1)/real(N_PC-1)
        else
            alpha_PCAC(i) = ht/d
        end if
        alpha = alpha_PCAC(i)
        X = alpha*d
        ETAT_PCAC(i) = 5
        if (alpha .lt. alpha_12) then
            PIVOT_PCAC(i) = 1
            SsTEND_PCAC(i) = -ssmax
            ScCOMP_PCAC(i) = (X/(d-X))*(ssmax/cequi)
        else
            PIVOT_PCAC(i) = 2
            ScCOMP_PCAC(i) = scmax
            SsTEND_PCAC(i) = scmax*cequi*(1-d/X)
        end if
        m2 = ScCOMP_PCAC(i)
        m1 = (SsTEND_PCAC(i)/cequi-m2)/d
        SsCOMP_PCAC(i) = (m1*(ht-d0)+m2)*cequi
        ScTEND_PCAC(i) = (m1*ht+m2)
        Ncc = ScCOMP_PCAC(i)*0.5*X*bw
        Mcc = ScCOMP_PCAC(i)*((1./4.)*ht*X-(1./6.)*X*X)*bw
        Dsys(1, 1) = SsCOMP_PCAC(i)
        Dsys(1, 2) = SsTEND_PCAC(i)
        Dsys(2, 1) = (d0-ht/2.)*SsCOMP_PCAC(i)
        Dsys(2, 2) = -(d-ht/2.)*SsTEND_PCAC(i)
        SOL(1) = Ncalc-Ncc
        SOL(2) = Mcalc-Mcc
        Dterm = Dsys(1, 1)*Dsys(2, 2)-Dsys(1, 2)*Dsys(2, 1)
        if (abs(Dterm) .gt. epsilon(Dterm)) then
            NUM = SOL(2)-Dsys(2, 1)*SOL(1)/Dsys(1, 1)
            DENUM = Dsys(2, 2)-Dsys(1, 2)*Dsys(2, 1)/Dsys(1, 1)
            AsTEND_PCAC(i) = NUM/DENUM
            AsCOMP_PCAC(i) = (SOL(1)-Dsys(1, 2)*AsTEND_PCAC(i))/Dsys(1, 1)
            RESIDU_PCAC(i) = 0
        elseif (abs(SsCOMP_PCAC(i)) .gt. epsilon(SsCOMP_PCAC(i))) then
            AsTEND_PCAC(i) = 0
            AsCOMP_PCAC(i) = SOL(2)/Dsys(2, 1)
            RESIDU_PCAC(i) = SOL(1)-Dsys(1, 1)*AsCOMP_PCAC(i)
        elseif (abs(SsTEND_PCAC(i)) .gt. epsilon(SsTEND_PCAC(i))) then
            AsCOMP_PCAC(i) = 0
            AsTEND_PCAC(i) = SOL(2)/Dsys(2, 2)
            RESIDU_PCAC(i) = SOL(1)-Dsys(1, 2)*AsTEND_PCAC(i)
        else
            AsCOMP_PCAC(i) = -1
            AsTEND_PCAC(i) = -1
            RESIDU_PCAC(i) = -1
        end if
    end do

    do i = 1, N_EC
        ScTEND_EC(i) = scmax*(real(i-1)/real(N_EC-1))
        ScCOMP_EC(i) = (scmaxc-ScTEND_EC(i))*(ht/yc-1)+scmaxc
        if (i .lt. N_EC) then
            X = (ScCOMP_EC(i)/(ScCOMP_EC(i)-ScTEND_EC(i)))*ht
            alpha_EC(i) = X/d
        else
            alpha_EC(i) = 1000
        end if
        ETAT_EC(i) = 6
        PIVOT_EC(i) = 3
        Ncc = 0.5*(ScCOMP_EC(i)+ScTEND_EC(i))*ht*bw
        Mcc = (1./12.)*(ScCOMP_EC(i)-ScTEND_EC(i))*ht*ht*bw
        SsCOMP_EC(i) = ((ScTEND_EC(i)-ScCOMP_EC(i))*(ht-d0)/ht+ScCOMP_EC(i))*cequi
        SsTEND_EC(i) = ((ScTEND_EC(i)-ScCOMP_EC(i))*d/ht+ScCOMP_EC(i))*cequi
        Dsys(1, 1) = SsCOMP_EC(i)
        Dsys(1, 2) = SsTEND_EC(i)
        Dsys(2, 1) = (d0-ht/2.)*SsCOMP_EC(i)
        Dsys(2, 2) = -(d-ht/2.)*SsTEND_EC(i)
        SOL(1) = Ncalc-Ncc
        SOL(2) = Mcalc-Mcc
        Dterm = Dsys(1, 1)*Dsys(2, 2)-Dsys(1, 2)*Dsys(2, 1)
        if (abs(Dterm) .gt. epsilon(Dterm)) then
            NUM = SOL(2)-Dsys(2, 1)*SOL(1)/Dsys(1, 1)
            DENUM = Dsys(2, 2)-Dsys(1, 2)*Dsys(2, 1)/Dsys(1, 1)
            AsTEND_EC(i) = NUM/DENUM
            AsCOMP_EC(i) = (SOL(1)-Dsys(1, 2)*AsTEND_EC(i))/Dsys(1, 1)
            RESIDU_EC(i) = 0
        elseif (abs(SsCOMP_EC(i)) .gt. epsilon(SsCOMP_EC(i))) then
            AsTEND_EC(i) = 0
            AsCOMP_EC(i) = SOL(2)/Dsys(2, 1)
            RESIDU_EC(i) = SOL(1)-Dsys(1, 1)*AsCOMP_EC(i)
        elseif (abs(SsTEND_EC(i)) .gt. epsilon(SsTEND_EC(i))) then
            AsCOMP_EC(i) = 0
            AsTEND_EC(i) = SOL(2)/Dsys(2, 2)
            RESIDU_EC(i) = SOL(1)-Dsys(1, 2)*AsTEND_EC(i)
        else
            AsCOMP_EC(i) = -1
            AsTEND_EC(i) = -1
            RESIDU_EC(i) = -1
        end if
    end do

    do i = 1, N_ET
        alpha_TOT(i) = alpha_ET(i)
        RESIDU_TOT(i) = RESIDU_ET(i)
        AsCOMP_TOT(i) = AsCOMP_ET(i)
        AsTEND_TOT(i) = AsTEND_ET(i)
        SsCOMP_TOT(i) = SsCOMP_ET(i)
        SsTEND_TOT(i) = SsTEND_ET(i)
        ScCOMP_TOT(i) = ScCOMP_ET(i)
        ScTEND_TOT(i) = ScTEND_ET(i)
        PIVOT_TOT(i) = PIVOT_ET(i)
        ETAT_TOT(i) = ETAT_ET(i)
    end do

    do i = 1, N_PCAC
        alpha_TOT(i+N_ET) = alpha_PCAC(i)
        RESIDU_TOT(i+N_ET) = RESIDU_PCAC(i)
        AsCOMP_TOT(i+N_ET) = AsCOMP_PCAC(i)
        AsTEND_TOT(i+N_ET) = AsTEND_PCAC(i)
        SsCOMP_TOT(i+N_ET) = SsCOMP_PCAC(i)
        SsTEND_TOT(i+N_ET) = SsTEND_PCAC(i)
        ScCOMP_TOT(i+N_ET) = ScCOMP_PCAC(i)
        ScTEND_TOT(i+N_ET) = ScTEND_PCAC(i)
        PIVOT_TOT(i+N_ET) = PIVOT_PCAC(i)
        ETAT_TOT(i+N_ET) = ETAT_PCAC(i)
    end do

    do i = 1, N_EC
        alpha_TOT(i+N_ET+N_PCAC) = alpha_EC(i)
        RESIDU_TOT(i+N_ET+N_PCAC) = RESIDU_EC(i)
        AsCOMP_TOT(i+N_ET+N_PCAC) = AsCOMP_EC(i)
        AsTEND_TOT(i+N_ET+N_PCAC) = AsTEND_EC(i)
        SsCOMP_TOT(i+N_ET+N_PCAC) = SsCOMP_EC(i)
        SsTEND_TOT(i+N_ET+N_PCAC) = SsTEND_EC(i)
        ScCOMP_TOT(i+N_ET+N_PCAC) = ScCOMP_EC(i)
        ScTEND_TOT(i+N_ET+N_PCAC) = ScTEND_EC(i)
        PIVOT_TOT(i+N_ET+N_PCAC) = PIVOT_EC(i)
        ETAT_TOT(i+N_ET+N_PCAC) = ETAT_EC(i)
    end do

    !Looking for EVENTUAL Roots
    i = 1
    indx_ch = 0

    do while (i .le. (N_TOT-1))

        COND_AJOUT_RESIDU = .false.
        COND_AJOUT_AsCOMP = .false.
        COND_AJOUT_AsTEND = .false.
        COND_AJOUT_Proceed = .false.

        aG = alpha_TOT(i)
        aD = alpha_TOT(i+1)
        fG = 0
        fD = 0

        if ((RESIDU_TOT(i)*RESIDU_TOT(i+1)) .lt. 0) then
            COND_AJOUT_RESIDU = .true.
            fG = fG+RESIDU_TOT(i)
            fD = fD+RESIDU_TOT(i+1)
        end if

        if ((AsCOMP_TOT(i)*AsCOMP_TOT(i+1)) .lt. 0) then
            COND_AJOUT_AsCOMP = .true.
            fG = fG+AsCOMP_TOT(i)
            fD = fD+AsCOMP_TOT(i+1)
        end if

        if ((AsTEND_TOT(i)*AsTEND_TOT(i+1)) .lt. 0) then
            COND_AJOUT_AsTEND = .true.
            fG = fG+AsTEND_TOT(i)
            fD = fD+AsTEND_TOT(i+1)
        end if

        if ((COND_AJOUT_RESIDU .eqv. (.true.)) &
            & .or. (COND_AJOUT_AsCOMP .eqv. (.true.)) &
            & .or. (COND_AJOUT_AsTEND .eqv. (.true.))) then

            alphaI = aG+(aD-aG)/(1+abs(fD/fG))
            X = alphaI*d

            if (alphaI .le. 0) then
                !ET

                sstend = -ssmax
                m1 = (ssmax/cequi)/(X-d)
                m2 = -m1*X
                sscomp = (m1*(ht-d0)+m2)*cequi
                sctend = m1*ht+m2
                sccomp = m1*0+m2
                astend = (Mcalc-Ncalc*(d0-0.5*ht))/(sstend*(-d-d0+ht))
                ascomp = (Ncalc-astend*SsTEND_ET(i))/sscomp
                pivot = 1
                etat = 4
                residu = 0

            elseif (alphaI .le. (ht/d)) then
                !PC

                etat = 5
                if (alphaI .lt. alpha_12) then
                    pivot = 1
                    sstend = -ssmax
                    sccomp = (X/(d-X))*(ssmax/cequi)
                else
                    pivot = 2
                    sccomp = scmax
                    sstend = scmax*cequi*(1-d/X)
                end if
                m2 = sccomp
                m1 = (sstend/cequi-m2)/d
                sscomp = (m1*(ht-d0)+m2)*cequi
                sctend = (m1*ht+m2)
                Ncc = sccomp*0.5*X*bw
                Mcc = sccomp*((1./4.)*ht*X-(1./6.)*X*X)*bw
                Dsys(1, 1) = sscomp
                Dsys(1, 2) = sstend
                Dsys(2, 1) = (d0-ht/2.)*sscomp
                Dsys(2, 2) = -(d-ht/2.)*sstend
                SOL(1) = Ncalc-Ncc
                SOL(2) = Mcalc-Mcc
                Dterm = Dsys(1, 1)*Dsys(2, 2)-Dsys(1, 2)*Dsys(2, 1)
                if (abs(Dterm) .gt. epsilon(Dterm)) then
                    NUM = SOL(2)-Dsys(2, 1)*SOL(1)/Dsys(1, 1)
                    DENUM = Dsys(2, 2)-Dsys(1, 2)*Dsys(2, 1)/Dsys(1, 1)
                    astend = NUM/DENUM
                    ascomp = (SOL(1)-Dsys(1, 2)*astend)/Dsys(1, 1)
                    residu = 0
                elseif (abs(sscomp) .gt. epsilon(sscomp)) then
                    astend = 0
                    ascomp = SOL(2)/Dsys(2, 1)
                    residu = SOL(1)-Dsys(1, 1)*ascomp
                elseif (abs(sstend) .gt. epsilon(sstend)) then
                    ascomp = 0
                    astend = SOL(2)/Dsys(2, 2)
                    residu = SOL(1)-Dsys(1, 2)*astend
                else
                    ascomp = -1
                    astend = -1
                    residu = -1
                end if

            else
                !EC

                m1 = scmaxc/(ht-yc-X)
                m2 = -m1*X
                sccomp = m1*0+m2
                sscomp = (m1*(ht-d0)+m2)*cequi
                sstend = (m1*d+m2)*cequi
                sctend = m1*ht+m2
                etat = 6
                pivot = 3
                Ncc = 0.5*(sccomp+sctend)*ht*bw
                Mcc = (1./12.)*(sccomp-sctend)*ht*ht*bw
                Dsys(1, 1) = sscomp
                Dsys(1, 2) = sstend
                Dsys(2, 1) = (d0-ht/2.)*sscomp
                Dsys(2, 2) = -(d-ht/2.)*sstend
                SOL(1) = Ncalc-Ncc
                SOL(2) = Mcalc-Mcc
                Dterm = Dsys(1, 1)*Dsys(2, 2)-Dsys(1, 2)*Dsys(2, 1)
                if (abs(Dterm) .gt. epsilon(Dterm)) then
                    NUM = SOL(2)-Dsys(2, 1)*SOL(1)/Dsys(1, 1)
                    DENUM = Dsys(2, 2)-Dsys(1, 2)*Dsys(2, 1)/Dsys(1, 1)
                    astend = NUM/DENUM
                    ascomp = (SOL(1)-Dsys(1, 2)*astend)/Dsys(1, 1)
                    residu = 0
                elseif (abs(sscomp) .gt. epsilon(sscomp)) then
                    astend = 0
                    ascomp = SOL(2)/Dsys(2, 1)
                    residu = SOL(1)-Dsys(1, 1)*ascomp
                elseif (abs(sstend) .gt. epsilon(sstend)) then
                    ascomp = 0
                    astend = SOL(2)/Dsys(2, 2)
                    residu = SOL(1)-Dsys(1, 2)*astend
                else
                    ascomp = -1
                    astend = -1
                    residu = -1
                end if

            end if

            !Dinstinguishing between ROOT and ASYMPTOT
            ! & Making sure that values are ACCEPTABLE

            if (COND_AJOUT_AsCOMP .eqv. (.true.)) then
                if (ascomp*AsCOMP_TOT(i) .gt. 0) then
                    if (.not. (abs(ascomp) .lt. abs(AsCOMP_TOT(i)))) then
                        goto 11
                    end if
                elseif (ascomp*AsCOMP_TOT(i+1) .gt. 0) then
                    if (.not. (abs(ascomp) .lt. abs(AsCOMP_TOT(i+1)))) then
                        goto 11
                    end if
                end if
            elseif ((AsCOMP_TOT(i) .lt. 0) .or. (AsCOMP_TOT(i+1) .lt. 0)) then
                goto 11
            end if

            if (COND_AJOUT_AsTEND .eqv. (.true.)) then
                if (astend*AsTEND_TOT(i) .gt. 0) then
                    if (.not. (abs(astend) .lt. abs(AsTEND_TOT(i)))) then
                        goto 11
                    end if
                elseif (astend*AsTEND_TOT(i+1) .gt. 0) then
                    if (.not. (abs(astend) .lt. abs(AsTEND_TOT(i+1)))) then
                        goto 11
                    end if
                end if
            elseif ((AsTEND_TOT(i) .lt. 0) .or. (AsTEND_TOT(i+1) .lt. 0)) then
                goto 11
            end if

            if (COND_AJOUT_RESIDU .eqv. (.true.)) then
                if (residu*RESIDU_TOT(i) .gt. 0) then
                    if (.not. (abs(residu) .lt. abs(RESIDU_TOT(i)))) then
                        goto 11
                    end if
                elseif (residu*RESIDU_TOT(i+1) .gt. 0) then
                    if (.not. (abs(residu) .lt. abs(RESIDU_TOT(i+1)))) then
                        goto 11
                    end if
                end if
            elseif ((abs(RESIDU_TOT(i)) .gt. 0) .or. (abs(RESIDU_TOT(i+1)) .gt. 0)) then
                goto 11
            end if

            COND_AJOUT_Proceed = .true.

11          continue

            if (COND_AJOUT_Proceed .eqv. (.true.)) then

                indx_ch = indx_ch+1
                N_TOT = N_TOT+1

                do j = 31, 40
                    call juveca(p(j), N_TOT)
                end do

                call jeveuo(p(31), 'E', vr=SsTEND_TOT)
                call jeveuo(p(32), 'E', vr=SsCOMP_TOT)
                call jeveuo(p(33), 'E', vr=ScTEND_TOT)
                call jeveuo(p(34), 'E', vr=ScCOMP_TOT)
                call jeveuo(p(35), 'E', vr=AsTEND_TOT)
                call jeveuo(p(36), 'E', vr=AsCOMP_TOT)
                call jeveuo(p(37), 'E', vr=alpha_TOT)
                call jeveuo(p(38), 'E', vr=RESIDU_TOT)
                call jeveuo(p(39), 'E', vi=PIVOT_TOT)
                call jeveuo(p(40), 'E', vi=ETAT_TOT)

                do s = 1, (N_TOT-i-1)
                    SsTEND_TOT(N_TOT-s+1) = SsTEND_TOT(N_TOT-s)
                    SsCOMP_TOT(N_TOT-s+1) = SsCOMP_TOT(N_TOT-s)
                    ScTEND_TOT(N_TOT-s+1) = ScTEND_TOT(N_TOT-s)
                    ScCOMP_TOT(N_TOT-s+1) = ScCOMP_TOT(N_TOT-s)
                    AsTEND_TOT(N_TOT-s+1) = AsTEND_TOT(N_TOT-s)
                    AsCOMP_TOT(N_TOT-s+1) = AsCOMP_TOT(N_TOT-s)
                    alpha_TOT(N_TOT-s+1) = alpha_TOT(N_TOT-s)
                    RESIDU_TOT(N_TOT-s+1) = RESIDU_TOT(N_TOT-s)
                    PIVOT_TOT(N_TOT-s+1) = PIVOT_TOT(N_TOT-s)
                    ETAT_TOT(N_TOT-s+1) = ETAT_TOT(N_TOT-s)
                end do

                if (COND_AJOUT_AsCOMP .eqv. (.true.)) then
                    ascomp = 0.0
                end if
                if (COND_AJOUT_AsTEND .eqv. (.true.)) then
                    astend = 0.0
                end if
                if (COND_AJOUT_RESIDU .eqv. (.true.)) then
                    residu = 0.0
                end if

                SsTEND_TOT(i+1) = sstend
                SsCOMP_TOT(i+1) = sscomp
                ScTEND_TOT(i+1) = sctend
                ScCOMP_TOT(i+1) = sccomp
                AsTEND_TOT(i+1) = astend
                AsCOMP_TOT(i+1) = ascomp
                alpha_TOT(i+1) = alphaI
                RESIDU_TOT(i+1) = residu
                PIVOT_TOT(i+1) = pivot
                ETAT_TOT(i+1) = etat

                i = i+2

            else

                i = i+1

            end if

        else

            i = i+1

        end if

    end do

    !Searching for OPTIMUM Solution
    COUNT_F = 0
    COUNT_CARA = 0
    COUNT_CARA_SYME = 0
    COUNT_CARA_NOCOMP = 0
    COUNT_CARA_SYMENOCOMP = 0
    AsTOT_F = 0.0

    do i = 1, N_TOT

        COND_COUNT = .false.
        COND_COUNT_SYME = .false.
        COND_COUNT_NOCOMP = .false.
        COND_COUNT_SYMENOCOMP = .false.

        if ((abs(RESIDU_TOT(i)) .lt. epsilon(RESIDU_TOT(i))) &
            & .AND. (AsCOMP_TOT(i) .ge. 0) .AND. (AsTEND_TOT(i) .ge. 0)) then

            COUNT_CARA = COUNT_CARA+1
            COND_COUNT = .true.
            !Checking NOCOMP
            if (((SsCOMP_TOT(i) .le. 0) .or. &
                & (AsCOMP_TOT(i) .lt. epsilon(AsCOMP_TOT(i)))) &
                & .AND. ((SsTEND_TOT(i) .le. 0) .or. &
                & (AsTEND_TOT(i) .lt. epsilon(AsTEND_TOT(i))))) then
                COUNT_CARA_NOCOMP = COUNT_CARA_NOCOMP+1
                COND_COUNT_NOCOMP = .true.
            end if
            !Checking SYME
            Calc = abs(AsCOMP_TOT(i)-AsTEND_TOT(i))
            if (Calc .le. slsyme) then
                COUNT_CARA_SYME = COUNT_CARA_SYME+1
                COND_COUNT_SYME = .true.
            end if
            !Checking SYME NOCOMP
            if ((COND_COUNT_NOCOMP .eqv. (.true.)) .and. &
                & (COND_COUNT_SYME .eqv. (.true.))) then
                COND_COUNT_SYMENOCOMP = .true.
                COUNT_CARA_SYMENOCOMP = COUNT_CARA_SYMENOCOMP+1
            end if
        end if

        if ((ferrsyme .eq. 0) .and. (ferrcomp .eq. 1)) then
            COND_F = COND_COUNT
            COUNT_F = COUNT_CARA
        elseif ((ferrsyme .eq. 1) .and. (ferrcomp .eq. 1)) then
            COND_F = COND_COUNT_SYME
            COUNT_F = COUNT_CARA_SYME
        elseif ((ferrsyme .eq. 0) .and. (ferrcomp .eq. 0)) then
            COND_F = COND_COUNT_NOCOMP
            COUNT_F = COUNT_CARA_NOCOMP
        elseif ((ferrsyme .eq. 1) .and. (ferrcomp .eq. 0)) then
            COND_F = COND_COUNT_SYMENOCOMP
            COUNT_F = COUNT_CARA_SYMENOCOMP
        end if

        if (COND_F .eqv. (.true.)) then

            AsCOMP_FOUND = AsCOMP_TOT(i)
            AsTEND_FOUND = AsTEND_TOT(i)
            AsTOT_FOUND = AsCOMP_FOUND+AsTEND_FOUND

            if ((COUNT_F .eq. 1) .or. (AsTOT_FOUND .lt. AsTOT_F)) then

                INDICE_F = i
                AsCOMP_F = AsCOMP_FOUND
                AsTEND_F = AsTEND_FOUND
                AsTOT_F = AsTOT_FOUND

            end if

        end if

    end do

    if (COUNT_F .gt. 0) then

        ascomp = AsCOMP_F
        astend = AsTEND_F
        sscomp = SsCOMP_TOT(INDICE_F)
        sstend = SsTEND_TOT(INDICE_F)
        sccomp = ScCOMP_TOT(INDICE_F)
        sctend = ScTEND_TOT(INDICE_F)
        alpha = alpha_TOT(INDICE_F)
        etat = ETAT_TOT(INDICE_F)
        pivot = PIVOT_TOT(INDICE_F)

    elseif ((ferrcomp .eq. 0) .and. (COUNT_CARA .gt. 0)) then

        etat = 0
        pivot = 0
        ascomp = -1
        astend = -1
        sscomp = -1
        sstend = -1
        sccomp = -1
        sctend = -1
        alpha = -1000
        ierr = 1

    elseif ((ferrsyme .eq. 1) .and. ((COUNT_CARA .gt. 0) .or. condns .eqv. (.true.))) then

        ascomp = -1
        astend = -1
        sscomp = -1
        sstend = -1
        sccomp = -1
        sctend = -1
        alpha = -1000
        etat = 0
        pivot = 0
        ierr = 2

    else

        ascomp = 0
        astend = 0
        sscomp = -1
        sstend = -1
        sccomp = -1
        sctend = -1
        alpha = -1000
        etat = 1
        pivot = 0

    end if

    do i = 1, 40
        call jedetr(p(i))
    end do

998 continue

end subroutine
