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

subroutine cafelsiter(cequi, effm, effn, ht, bw, &
                      enrobi, enrobs, scmaxi, scmaxs, ssmax, &
                      ferrsyme, slsyme, uc, condns, &
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
!      I FERRSYME     FERRAILLAGE SYMETRIQUE?
!                        FERRSYME = 0 (NON)
!                        FERRSYME = 1 (OUI)
!      I SLSYME       SECTION SEUIL DE TOLERANCE POUR UN FERRAILLAGE SYMETRIQUE
!      I UC           UNITE DES CONTRAINTES :
!                        UC = 0 CONTRAINTES EN Pa
!                        UC = 1 CONTRAINTES EN MPa
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
#include "asterfort/verifels.h"
#include "asterfort/wkvect.h"
#include "asterfort/jedetr.h"
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
    integer :: ferrsyme
    real(kind=8) :: slsyme
    integer :: uc
    logical :: condns
    real(kind=8) :: astend
    real(kind=8) :: ascomp
    real(kind=8) :: sstend
    real(kind=8) :: sscomp
    real(kind=8) :: sctend
    real(kind=8) :: sccomp
    real(kind=8) :: alpha
    integer :: pivot
    integer :: etat
    integer :: ierr

!-----------------------------------------------------------------------
!!!!VARIABLES DE CALCUL
!-----------------------------------------------------------------------
    real(kind=8) :: d, d0, NUM, DENUM, Diffa
    real(kind=8) :: Mcalc, Ncalc, scmax, scmaxc, alpha_12, unite_pa
    real(kind=8) :: AsTEND_F, AsCOMP_F, AsTOT_F, AsTEND_FOUND, AsCOMP_FOUND, AsTOT_FOUND
    real(kind=8) :: X, Ncc, Mcc, Dterm, Calc, yc
    logical :: COND_AJOUT_RESIDU, COND_AJOUT_AsCOMP, COND_AJOUT_AsTEND, COND_COUNT
    logical :: COND_COUNT_SYME
    integer :: COUNT_CARA, verif
    integer :: COUNT_CARA_SYME
    logical :: COND_F, PIV_C
    integer :: COUNT_F
    integer :: N_ET, N_PC, N_PCAC, N_EC, N_TOT

    real(kind=8), pointer :: SsTEND_ET(:) => null(), SsCOMP_ET(:) => null()
    real(kind=8), pointer :: ScTEND_ET(:) => null(), ScCOMP_ET(:) => null()
    real(kind=8), pointer :: AsTEND_ET(:) => null(), AsCOMP_ET(:) => null()
    real(kind=8), pointer :: alpha_ET(:) => null(), RESIDU_ET(:) => null()
    integer, pointer :: PIVOT_ET(:) => null(), ETAT_ET(:) => null()

    real(kind=8), pointer :: SsTEND_PCAC(:) => null(), SsCOMP_PCAC(:) => null()
    real(kind=8), pointer :: ScTEND_PCAC(:) => null(), ScCOMP_PCAC(:) => null()
    real(kind=8), pointer :: AsTEND_PCAC(:) => null(), AsCOMP_PCAC(:) => null()
    real(kind=8), pointer :: alpha_PCAC(:) => null(), RESIDU_PCAC(:) => null()
    integer, pointer :: PIVOT_PCAC(:) => null(), ETAT_PCAC(:) => null()

    real(kind=8), pointer :: SsTEND_EC(:) => null(), SsCOMP_EC(:) => null()
    real(kind=8), pointer :: ScTEND_EC(:) => null(), ScCOMP_EC(:) => null()
    real(kind=8), pointer :: AsTEND_EC(:) => null(), AsCOMP_EC(:) => null()
    real(kind=8), pointer :: alpha_EC(:) => null(), RESIDU_EC(:) => null()
    integer, pointer :: PIVOT_EC(:) => null(), ETAT_EC(:) => null()

    real(kind=8), pointer :: SsTEND_TOT(:) => null(), SsCOMP_TOT(:) => null()
    real(kind=8), pointer :: ScTEND_TOT(:) => null(), ScCOMP_TOT(:) => null()
    real(kind=8), pointer :: AsTEND_TOT(:) => null(), AsCOMP_TOT(:) => null()
    real(kind=8), pointer :: alpha_TOT(:) => null(), RESIDU_TOT(:) => null()
    integer, pointer :: PIVOT_TOT(:) => null(), ETAT_TOT(:) => null()
    real(kind=8), pointer :: Pente_AsTEND_TOT(:) => null(), Pente_AsCOMP_TOT(:) => null()
    real(kind=8), pointer :: PENTE_RESIDU_TOT(:) => null()

    real(kind=8), dimension(2, 2) :: Dsys
    real(kind=8), dimension(2) :: SOL
    integer :: i, INDICE_F
    real(kind=8) :: INDICE_F_AVATAR

    character(8), dimension(43) :: p

!-----------------------------------------------------------------------
!!!!LANCEMENT DU CALCUL
!-----------------------------------------------------------------------

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

    if (uc .eq. 0) then
        unite_pa = 1.e-6
    elseif (uc .eq. 1) then
        unite_pa = 1.
    end if

    alpha_12 = 1.0/(1.0+(ssmax/cequi)/scmax)

    Mcalc = abs(effm)
    Ncalc = effn
    N_ET = 11
    N_PC = 101
    N_PCAC = CEILING((N_PC-1)*(ht/d))+1

    !Determination Pivot C 'ELS'
    PIV_C = .false.
    if (scmaxc .lt. scmax) then
        PIV_C = .true.
        yc = (scmaxc/scmax)*ht
    end if
    N_EC = CEILING(10*(scmaxc*unite_pa))+1

    N_TOT = N_ET+N_PCAC+N_EC

    do i = 1, 43
        write (p(i), '(A6,I2)') 'POINT_', i
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

    call wkvect(p(41), ' V V R ', N_TOT, vr=PENTE_RESIDU_TOT)
    call wkvect(p(42), ' V V R ', N_TOT, vr=PENTE_AsCOMP_TOT)
    call wkvect(p(43), ' V V R ', N_TOT, vr=PENTE_AsTEND_TOT)

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
        if (abs(alpha) .lt. epsilon(alpha)) then
            PIVOT_PCAC(i) = 1
            RESIDU_PCAC(i) = Mcalc+Ncalc*(d-0.5*ht)
            SsTEND_PCAC(i) = -ssmax
            ScCOMP_PCAC(i) = 0
            ScTEND_PCAC(i) = -(ssmax/cequi)*(ht/d)
            SsCOMP_PCAC(i) = -ssmax*(ht-d0)/d
        elseif ((alpha .gt. 0) .AND. (alpha .lt. alpha_12)) then
            PIVOT_PCAC(i) = 1
            SsTEND_PCAC(i) = -ssmax
            ScCOMP_PCAC(i) = (X/(d-X))*(ssmax/cequi)
            ScTEND_PCAC(i) = ScCOMP_PCAC(i)*(1-ht/X)
            SsCOMP_PCAC(i) = ScCOMP_PCAC(i)*(1-(ht-d0)/X)*cequi
        else
            PIVOT_PCAC(i) = 2
            ScCOMP_PCAC(i) = scmax
            SsTEND_PCAC(i) = scmax*cequi*(1-d/X)
            ScTEND_PCAC(i) = scmax*(1-ht/X)
            SsCOMP_PCAC(i) = scmax*cequi*(1-(ht-d0)/X)
        end if
        Ncc = ScCOMP_PCAC(i)*0.5*X*bw
        Mcc = ScCOMP_PCAC(i)*((1./4.)*ht*X-(1./6.)*X*X)*bw
        Dsys(1, 1) = SsCOMP_PCAC(i)
        Dsys(1, 2) = SsTEND_PCAC(i)
        Dsys(2, 1) = (d0-ht/2.)*SsCOMP_PCAC(i)
        Dsys(2, 2) = -(d-ht/2.)*SsTEND_PCAC(i)
        SOL(1) = Ncalc-Ncc
        SOL(2) = Mcalc-Mcc
        Dterm = Dsys(1, 1)*Dsys(2, 2)-Dsys(1, 2)*Dsys(2, 1)
        !if (Dterm.ne.0) then
        if (abs(Dterm) .gt. epsilon(Dterm)) then
            NUM = SOL(2)-Dsys(2, 1)*SOL(1)/Dsys(1, 1)
            DENUM = Dsys(2, 2)-Dsys(1, 2)*Dsys(2, 1)/Dsys(1, 1)
            AsTEND_PCAC(i) = NUM/DENUM
            AsCOMP_PCAC(i) = (SOL(1)-Dsys(1, 2)*AsTEND_PCAC(i))/Dsys(1, 1)
            RESIDU_PCAC(i) = 0
            !elseif (SsCOMP_PCAC(i).ne.0) then
        elseif (abs(SsCOMP_PCAC(i)) .gt. epsilon(SsCOMP_PCAC(i))) then
            AsTEND_PCAC(i) = 0
            AsCOMP_PCAC(i) = SOL(2)/Dsys(2, 1)
            RESIDU_PCAC(i) = SOL(1)-Dsys(1, 1)*AsCOMP_PCAC(i)
            !elseif (SsTEND_PCAC(i).ne.0) then
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
        if (PIV_C .eqv. (.true.)) then
            ScTEND_EC(i) = scmaxc*(real(i-1)/real(N_EC-1))
            ScCOMP_EC(i) = (scmaxc-ScTEND_EC(i))*(ht/yc-1)+scmaxc
        else
            ScTEND_EC(i) = scmax*(real(i-1)/real(N_EC-1))
            ScCOMP_EC(i) = scmax
        end if
        if (i .lt. N_EC) then
            X = (ScCOMP_EC(i)/(ScCOMP_EC(i)-ScTEND_EC(i)))*ht
            alpha_EC(i) = X/d
        else
            alpha_EC(i) = -1000
        end if
        ETAT_EC(i) = 6
        PIVOT_EC(i) = 2
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
        !if (Dterm.ne.0) then
        if (abs(Dterm) .gt. epsilon(Dterm)) then
            NUM = SOL(2)-Dsys(2, 1)*SOL(1)/Dsys(1, 1)
            DENUM = Dsys(2, 2)-Dsys(1, 2)*Dsys(2, 1)/Dsys(1, 1)
            AsTEND_EC(i) = NUM/DENUM
            AsCOMP_EC(i) = (SOL(1)-Dsys(1, 2)*AsTEND_EC(i))/Dsys(1, 1)
            RESIDU_EC(i) = 0
            !elseif (SsCOMP_EC(i).ne.0) then
        elseif (abs(SsCOMP_EC(i)) .gt. epsilon(SsCOMP_EC(i))) then
            AsTEND_EC(i) = 0
            AsCOMP_EC(i) = SOL(2)/Dsys(2, 1)
            RESIDU_EC(i) = SOL(1)-Dsys(1, 1)*AsCOMP_EC(i)
            !elseif (SsTEND_EC(i).ne.0) then
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

    do i = 1, N_TOT
        if (i .eq. 1) then
            Diffa = 0
        else
            Diffa = alpha_TOT(i-1)-alpha_TOT(i)
        end if
        if (i .eq. 1) then
            PENTE_RESIDU_TOT(i) = 0
            PENTE_AsCOMP_TOT(i) = 0
            PENTE_AsTEND_TOT(i) = 0
            !elseif (alpha_TOT(i-1).ne.alpha_TOT(i)) then
        elseif (abs(Diffa) .gt. (epsilon(Diffa))) then
            PENTE_RESIDU_TOT(i) = RESIDU_TOT(i)-RESIDU_TOT(i-1)
            PENTE_RESIDU_TOT(i) = (PENTE_RESIDU_TOT(i))/(alpha_TOT(i)-alpha_TOT(i-1))
            PENTE_AsCOMP_TOT(i) = (AsCOMP_TOT(i)-AsCOMP_TOT(i-1))
            PENTE_AsCOMP_TOT(i) = (PENTE_AsCOMP_TOT(i))/(alpha_TOT(i)-alpha_TOT(i-1))
            PENTE_AsTEND_TOT(i) = (AsTEND_TOT(i)-AsTEND_TOT(i-1))
            PENTE_AsTEND_TOT(i) = (PENTE_AsTEND_TOT(i))/(alpha_TOT(i)-alpha_TOT(i-1))
        else
            PENTE_RESIDU_TOT(i) = 0
            PENTE_AsCOMP_TOT(i) = 0
            PENTE_AsTEND_TOT(i) = 0
        end if
    end do

    COUNT_CARA = 0
    COUNT_CARA_SYME = 0

    do i = 1, N_TOT
        COND_AJOUT_RESIDU = .false.
        COND_AJOUT_AsCOMP = .false.
        COND_AJOUT_AsTEND = .false.
        COND_COUNT = .false.
        COND_COUNT_SYME = .false.
        if ((abs(RESIDU_TOT(i)) .lt. epsilon(RESIDU_TOT(i))) &
             & .AND. (AsCOMP_TOT(i) .ge. 0) .AND. (AsTEND_TOT(i) .ge. 0)) then
            COUNT_CARA = COUNT_CARA+1
            COND_COUNT = .true.
            Calc = abs(AsCOMP_TOT(i)-AsTEND_TOT(i))
            if (Calc .le. slsyme) then
                COUNT_CARA_SYME = COUNT_CARA_SYME+1
                COND_COUNT_SYME = .true.
            end if
        elseif (i .gt. 1) then
            if ((AsCOMP_TOT(i) .ge. 0) .AND. (AsTEND_TOT(i) .ge. 0) &
                 & .AND. (AsCOMP_TOT(i-1) .ge. 0) .AND. (AsTEND_TOT(i-1) .ge. 0)) then
                if (((RESIDU_TOT(i)*RESIDU_TOT(i-1)) .le. 0) &
                      & .AND. ((PENTE_RESIDU_TOT(i)*PENTE_RESIDU_TOT(i-1)) .ge. 0)) then
                    COUNT_CARA = COUNT_CARA+1
                    COND_COUNT = .true.
                    COND_AJOUT_RESIDU = .true.
                    Calc = abs(AsCOMP_TOT(i)-AsTEND_TOT(i))
                    if (Calc .le. slsyme) then
                        COUNT_CARA_SYME = COUNT_CARA_SYME+1
                        COND_COUNT_SYME = .true.
                    end if
                end if
            elseif ((AsCOMP_TOT(i) .ge. 0) .AND. (AsCOMP_TOT(i-1) .ge. 0)) then
                if (((RESIDU_TOT(i)*RESIDU_TOT(i-1)) .le. 0) &
                      & .AND. ((PENTE_RESIDU_TOT(i)*PENTE_RESIDU_TOT(i-1)) .ge. 0) &
                      & .AND. ((AsTEND_TOT(i)*AsTEND_TOT(i-1)) .le. 0) &
                      & .AND. ((PENTE_AsTEND_TOT(i)*PENTE_AsTEND_TOT(i-1)) .ge. 0)) then
                    COUNT_CARA = COUNT_CARA+1
                    COND_COUNT = .true.
                    COND_AJOUT_AsTEND = .true.
                    Calc = abs(AsCOMP_TOT(i)-AsTEND_TOT(i))
                    if (Calc .le. slsyme) then
                        COUNT_CARA_SYME = COUNT_CARA_SYME+1
                        COND_COUNT_SYME = .true.
                    end if
                end if
            elseif ((AsTEND_TOT(i) .ge. 0) .AND. (AsTEND_TOT(i-1) .ge. 0)) then
                if (((RESIDU_TOT(i)*RESIDU_TOT(i-1)) .le. 0) &
                      & .AND. ((PENTE_RESIDU_TOT(i)*PENTE_RESIDU_TOT(i-1)) .ge. 0) &
                      & .AND. ((AsCOMP_TOT(i)*AsCOMP_TOT(i-1)) .le. 0) &
                      & .AND. ((PENTE_AsCOMP_TOT(i)*PENTE_AsCOMP_TOT(i-1)) .ge. 0)) then
                    COUNT_CARA = COUNT_CARA+1
                    COND_COUNT = .true.
                    COND_AJOUT_AsCOMP = .true.
                    Calc = abs(AsCOMP_TOT(i)-AsTEND_TOT(i))
                    if (Calc .le. slsyme) then
                        COUNT_CARA_SYME = COUNT_CARA_SYME+1
                        COND_COUNT_SYME = .true.
                    end if
                end if
            else
                if (((RESIDU_TOT(i)*RESIDU_TOT(i-1)) .le. 0) &
                      & .AND. ((PENTE_RESIDU_TOT(i)*PENTE_RESIDU_TOT(i-1)) .ge. 0) &
                      & .AND. ((AsCOMP_TOT(i)*AsCOMP_TOT(i-1)) .le. 0) &
                      & .AND. ((PENTE_AsCOMP_TOT(i)*PENTE_AsCOMP_TOT(i-1)) .ge. 0) &
                      & .AND. ((AsTEND_TOT(i)*AsTEND_TOT(i-1)) .le. 0) &
                      & .AND. ((PENTE_AsTEND_TOT(i)*PENTE_AsTEND_TOT(i-1)) .ge. 0)) then
                    COUNT_CARA = COUNT_CARA+1
                    COND_COUNT = .true.
                    COND_AJOUT_AsCOMP = .true.
                    COND_AJOUT_AsTEND = .true.
                    Calc = abs(AsCOMP_TOT(i)-AsTEND_TOT(i))
                    if (Calc .le. slsyme) then
                        COUNT_CARA_SYME = COUNT_CARA_SYME+1
                        COND_COUNT_SYME = .true.
                    end if
                end if
            end if
        end if

        if (ferrsyme .eq. 0) then
            COND_F = COND_COUNT
            COUNT_F = COUNT_CARA
        elseif (ferrsyme .eq. 1) then
            COND_F = COND_COUNT_SYME
            COUNT_F = COUNT_CARA_SYME
        end if

        if (COND_F .eqv. (.true.)) then
            if (COUNT_F .eq. 1) then
                if ((COND_AJOUT_RESIDU .eqv. (.false.)) &
                     & .AND. (COND_AJOUT_AsCOMP .eqv. (.false.)) &
                     & .AND. (COND_AJOUT_AsTEND .eqv. (.false.))) then
                    AsCOMP_F = AsCOMP_TOT(i)
                    AsTEND_F = AsTEND_TOT(i)
                    AsTOT_F = AsCOMP_F+AsTEND_F
                    INDICE_F = i
                    INDICE_F_AVATAR = i
                else
                    INDICE_F_AVATAR = i-0.5
                    INDICE_F = i
                    if ((COND_AJOUT_AsCOMP .eqv. (.true.)) &
                         & .AND. (COND_AJOUT_AsTEND .eqv. (.true.))) then
                        AsCOMP_F = 0
                        AsTEND_F = 0
                    elseif (COND_AJOUT_AsTEND .eqv. (.true.)) then
                        AsTEND_F = 0
                        AsCOMP_F = 0.5*(AsCOMP_TOT(i-1)+AsCOMP_TOT(i))
                    elseif (COND_AJOUT_AsCOMP .eqv. (.true.)) then
                        AsTEND_F = 0.5*(AsTEND_TOT(i-1)+AsTEND_TOT(i))
                        AsCOMP_F = 0
                    else
                        AsTEND_F = 0.5*(AsTEND_TOT(i-1)+AsTEND_TOT(i))
                        AsCOMP_F = 0.5*(AsCOMP_TOT(i-1)+AsCOMP_TOT(i))
                    end if
                    AsTOT_F = AsTEND_F+AsCOMP_F
                end if
            else
                if ((COND_AJOUT_RESIDU .eqv. (.false.)) &
                     & .AND. (COND_AJOUT_AsCOMP .eqv. (.false.)) &
                     & .AND. (COND_AJOUT_AsTEND .eqv. (.false.))) then
                    AsCOMP_FOUND = AsCOMP_TOT(i)
                    AsTEND_FOUND = AsTEND_TOT(i)
                    AsTOT_FOUND = AsCOMP_FOUND+AsTEND_FOUND
                    if (AsTOT_FOUND .lt. AsTOT_F) then
                        INDICE_F = i
                        INDICE_F_AVATAR = i
                        AsCOMP_F = AsCOMP_FOUND
                        AsTEND_F = AsTEND_FOUND
                        AsTOT_F = AsTOT_FOUND
                    end if
                else
                    if ((COND_AJOUT_RESIDU .eqv. (.false.)) &
                         & .AND. (COND_AJOUT_AsCOMP .eqv. (.false.)) &
                         & .AND. (COND_AJOUT_AsTEND .eqv. (.false.))) then
                        AsCOMP_FOUND = 0
                        AsTEND_FOUND = 0
                    elseif (COND_AJOUT_AsTEND .eqv. (.true.)) then
                        AsTEND_FOUND = 0
                        AsCOMP_FOUND = 0.5*(AsCOMP_TOT(i-1)+AsCOMP_TOT(i))
                    elseif (COND_AJOUT_AsCOMP .eqv. (.true.)) then
                        AsTEND_FOUND = 0.5*(AsTEND_TOT(i-1)+AsTEND_TOT(i))
                        AsCOMP_FOUND = 0
                    elseif (COND_AJOUT_RESIDU .eqv. (.true.)) then
                        AsTEND_FOUND = 0.5*(AsTEND_TOT(i-1)+AsTEND_TOT(i))
                        AsCOMP_FOUND = 0.5*(AsCOMP_TOT(i-1)+AsCOMP_TOT(i))
                    end if
                    AsTOT_FOUND = AsCOMP_FOUND+AsTEND_FOUND
                    if (AsTOT_FOUND .lt. AsTOT_F) then
                        INDICE_F = i
                        INDICE_F_AVATAR = i-0.5
                        AsCOMP_F = AsCOMP_FOUND
                        AsTEND_F = AsTEND_FOUND
                        AsTOT_F = AsTOT_FOUND
                    end if
                end if
            end if
        end if
    end do

    if (((ferrsyme .eq. 0) .and. (COUNT_CARA .gt. 0)) &
          & .or. ((ferrsyme .eq. 1) .and. (COUNT_CARA_SYME .gt. 0))) then

        Diffa = CEILING(INDICE_F_AVATAR)-INDICE_F_AVATAR

        if (abs(Diffa) .lt. epsilon(Diffa)) then
            ascomp = AsCOMP_F
            astend = AsTEND_F
            sscomp = SsCOMP_TOT(INDICE_F)
            sstend = SsTEND_TOT(INDICE_F)
            sccomp = ScCOMP_TOT(INDICE_F)
            sctend = ScTEND_TOT(INDICE_F)
            alpha = alpha_TOT(INDICE_F)
            etat = ETAT_TOT(INDICE_F)
            pivot = PIVOT_TOT(INDICE_F)
        else
            ascomp = AsCOMP_F
            astend = AsTEND_F
            sscomp = 0.5*(SsCOMP_TOT(INDICE_F)+SsCOMP_TOT(INDICE_F-1))
            sstend = 0.5*(SsTEND_TOT(INDICE_F)+SsTEND_TOT(INDICE_F-1))
            sccomp = 0.5*(ScCOMP_TOT(INDICE_F)+ScCOMP_TOT(INDICE_F-1))
            sctend = 0.5*(ScTEND_TOT(INDICE_F)+ScTEND_TOT(INDICE_F-1))
            alpha = 0.5*(alpha_TOT(INDICE_F)+alpha_TOT(INDICE_F-1))
            etat = ETAT_TOT(INDICE_F)
            pivot = PIVOT_TOT(INDICE_F)
        end if

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

    !! Vérification par le diagramme d'interaction avec dnsyi=dnsys=0.0
        call verifels(cequi, ht, bw, enrobi, enrobs, &
                      scmaxi, scmaxs, ssmax, uc, &
                      0.0, 0.0, effm, effn, verif)

        if (verif .eq. 0) then
           !!OK
            ascomp = 0
            astend = 0
            sscomp = -1
            sstend = -1
            sccomp = -1
            sctend = -1
            alpha = -1000
            etat = 1
            pivot = 0
        else
            ascomp = -1
            astend = -1
            sscomp = -1
            sstend = -1
            sccomp = -1
            sctend = -1
            alpha = -1000
            etat = 0
            pivot = 0
            ierr = 1
        end if

    end if

    do i = 1, 43
        call jedetr(p(i))
    end do

end subroutine
