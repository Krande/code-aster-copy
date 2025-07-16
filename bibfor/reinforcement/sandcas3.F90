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

subroutine sandcas3(effrts, ht, enrobi, enrobs, facier, fbeton, gammas, gammac, &
                    thiter, cond109, ferrcomp, ferrsyme, slsyme, uc, &
                    dnsxi, dnsxs, dnsyi, dnsys, etsxi, etsxs, etsyi, etsys, &
                    snsxi, snsxs, snsyi, snsys, ncmaxi, ncmini, ncmaxs, ncmins, &
                    t_inf, t_sup, theta_inf, theta_sup, ierr)

!______________________________________________________________________
!
!      SANDCAS3
!
!      CALCUL DES ACIERS A L'ELU PAR LA METHODE SANDWICH
!      CAS 3 - FERRAILLAGE [+] REQUIS EN INF
!
!      I EFFRTS        (DIM 6) TORSEUR DES EFFORTS, MOMENTS, ...
!                         EFFRTS(1) = NXX
!                         EFFRTS(2) = NYY
!                         EFFRTS(3) = NXY
!                         EFFRTS(4) = MXX
!                         EFFRTS(5) = MYY
!                         EFFRTS(6) = MXY
!      I HT            EPAISSEUR DE LA COQUE
!      I ENROBI        ENROBAGE DES ARMATURES INFERIEURES
!      I ENROBS        ENROBAGE DES ARMATURES SUPERIEURES
!      I FACIER        LIMITE D'ELASTICITE DES ACIERS (CONTRAINTE)
!      I FBETON        RESISTANCE EN COMPRESSION DU BETON (CONTRAINTE)
!      I GAMMAS        COEFFICIENT DE SECURITE SUR LA RESISTANCE
!                      DE CALCUL DES ACIERS
!      I GAMMAC        COEFFICIENT DE SECURITE SUR LA RESISTANCE
!                      DE CALCUL DU BETON
!      I THITER        ANGLE D'ITERATION POUR LA RECHERCHE DES SOLUTIONS
!      I COND109       APPLICATION DE LA CLAUSE §109 POUR LA VALORISATION
!                      DE LA RESISTANCE DES MEMBRANES
!                         COND109 = 0 (NON)
!                         COND109 = 1 (OUI)
!      I FERRCOMP      PRISE EN COMPTE DU FERRAILLAGE DE COMPRESSION
!                         FERRCOMP = 0 (NON)
!                         FERRCOMP = 1 (OUI)
!      I FERRSYME      FERRAILLAGE SYMETRIQUE?
!                         FERRSYME = 0 (NON)
!                         FERRSYME = 1 (OUI)
!      I SLSYME        SECTION SEUIL DE TOLERANCE POUR UN
!                      FERRAILLAGE SYMETRIQUE
!      I UC            UNITE DES CONTRAINTES :
!                         UC = 0 CONTRAINTES EN Pa
!                         UC = 1 CONTRAINTES EN MPa
!
!      O DNSXI         DENSITE DE L'ACIER INFERIEUR SUIVANT X
!      O DNSXS         DENSITE DE L'ACIER SUPERIEUR SUIVANT X
!      O DNSYI         DENSITE DE L'ACIER INFERIEUR SUIVANT Y
!      O DNSYS         DENSITE DE L'ACIER SUPERIEUR SUIVANT Y
!      O ETSXI         ETAT DE L'ACIER INFERIEUR SUIVANT X
!      O ETSXS         ETAT DE L'ACIER SUPERIEUR SUIVANT X
!      O ETSYI         ETAT DE L'ACIER INFERIEUR SUIVANT Y
!      O ETSYS         ETAT DE L'ACIER SUPERIEUR SUIVANT Y
!                         ETAT = 0 (TRACTION)
!                         ETAT = 1 (COMPRESSION)
!                         ETAT = 2 (PAS D'EQUILIBRE)
!      O SNSXI         CONTRAINTE DE L'ACIER INFERIEUR SUIVANT X
!      O SNSXS         CONTRAINTE DE L'ACIER SUPERIEUR SUIVANT X
!      O SNSYI         CONTRAINTE DE L'ACIER INFERIEUR SUIVANT Y
!      O SNSYS         CONTRAINTE DE L'ACIER SUPERIEUR SUIVANT Y
!      O NCMAXI        CONTRAINTE PRINCIPALE MAJEURE
!                      DE COMPRESSION DANS LA COUCHE INF DE BETON
!      O NCMINI        CONTRAINTE PRINCIPALE MINEURE
!                      DE COMPRESSION DANS LA COUCHE INF DE BETON
!      O NCMAXS        CONTRAINTE PRINCIPALE MAJEURE
!                      DE COMPRESSION DANS LA COUCHE SUP DE BETON
!      O NCMINS        CONTRAINTE PRINCIPALE MINEURE
!                      DE COMPRESSION DANS LA COUCHE SUP DE BETON
!      O T_INF         EPAISSEUR DE LA COUCHE INFERIEURE
!      O T_SUP         EPAISSEUR DE LA COUCHE SUPERIEURE
!      O THETA_INF     ANGLE D'INCLINAISON DE LA BIELLE PRINCIPALE
!                      DE COMPRESSION DANS LA COUCHE INFERIEURE
!      O THETA_SUP     ANGLE D'INCLINAISON DE LA BIELLE PRINCIPALE
!                      DE COMPRESSION DANS LA COUCHE SUPERIEURE
!      O IERR          CODE RETOUR (0 = OK)
!
!______________________________________________________________________
!
    use sand_solvers_module
    implicit none
#include "asterc/r8pi.h"
#include "asterc/r8dgrd.h"
#include "asterfort/wkvect.h"
#include "asterfort/jedetr.h"

!VARIABLES PRINCIPALES
!------------------------------
    real(kind=8) :: effrts(6)
    real(kind=8) :: ht
    real(kind=8) :: enrobi
    real(kind=8) :: enrobs
    real(kind=8) :: facier
    real(kind=8) :: fbeton
    real(kind=8) :: gammas
    real(kind=8) :: gammac
    real(kind=8) :: thiter
    integer(kind=8) :: cond109
    integer(kind=8) :: ferrcomp
    integer(kind=8) :: ferrsyme
    real(kind=8) :: slsyme
    integer(kind=8) :: uc
    real(kind=8) :: dnsxi
    real(kind=8) :: dnsxs
    real(kind=8) :: dnsyi
    real(kind=8) :: dnsys
    integer(kind=8) :: etsxi
    integer(kind=8) :: etsxs
    integer(kind=8) :: etsyi
    integer(kind=8) :: etsys
    real(kind=8) :: snsxi
    real(kind=8) :: snsxs
    real(kind=8) :: snsyi
    real(kind=8) :: snsys
    real(kind=8) :: ncmaxi
    real(kind=8) :: ncmini
    real(kind=8) :: ncmaxs
    real(kind=8) :: ncmins
    real(kind=8) :: t_inf
    real(kind=8) :: t_sup
    real(kind=8) :: theta_inf
    real(kind=8) :: theta_sup
    integer(kind=8) :: ierr

!Variables de calcul
    real(kind=8) :: fcd, fyd, ySUP, yINF, Z, fcd1, fc, pi, vect(20), my_epsi
    real(kind=8) :: Nxx, Nyy, Nxy, Mxx, Myy, Mxy
    integer(kind=8) :: N_SUP, N_INF, N_TOT, i
    real(kind=8) :: unite_pa, Calc1, Calc2
    character(20) :: p(15)
    real(kind=8), pointer :: nSX_SUP(:) => null(), nSX_INF(:) => null()
    real(kind=8), pointer :: nSY_SUP(:) => null(), nSY_INF(:) => null()
    real(kind=8), pointer :: nS_TOT(:) => null()
    real(kind=8), pointer :: tSUP(:) => null(), tINF(:) => null()
    real(kind=8), pointer :: ncMAX_SUP(:) => null(), ncMIN_SUP(:) => null()
    real(kind=8), pointer :: ncMAX_INF(:) => null(), ncMIN_INF(:) => null()
    real(kind=8), pointer :: RESIDU(:) => null()
    real(kind=8), pointer :: AngleSUP(:) => null(), AngleINF(:) => null()
    real(kind=8) :: tSUPA, tSUPB, tSUPmed, RESIDUA, RESIDUB, RESIDUmed
    logical :: condmed

    pi = r8pi()
    ierr = 0
    fcd = fbeton/gammac
    fyd = facier/gammas
    ySUP = 0.5*ht-enrobs
    yINF = 0.5*ht-enrobi
    Z = ySUP+yINF
    Nxx = effrts(1)
    Nyy = effrts(2)
    Nxy = effrts(3)
    Mxx = effrts(4)
    Myy = effrts(5)
    Mxy = effrts(6)

!-------------------------------------------------------------------------------------------------
!On commence à itérer sur le MODELE SANDWICH A ADOPTER (Determination des epaisseurs des 3 Couches)
!-------------------------------------------------------------------------------------------------

    N_INF = 2*ceiling(90.0/thiter)+1
    N_SUP = 1
    N_TOT = N_SUP*N_INF

    if (uc .eq. 0) then
        unite_pa = 1.e-6
    elseif (uc .eq. 1) then
        unite_pa = 1.
    end if

    if (cond109 .eq. 1) then
        fcd1 = 0.6*(1-fbeton*unite_pa/250.0)*fcd
    else
        fcd1 = fcd
    end if

    my_epsi = 0.01/unite_pa

    fc = fcd1

    do i = 1, 15
        write (p(i), fmt='(A9,I2)') 'SANDCAS3_', i
    end do

    call wkvect(p(1), ' V V R ', N_TOT, vr=nSX_INF)
    call wkvect(p(2), ' V V R ', N_TOT, vr=nSY_INF)
    call wkvect(p(3), ' V V R ', N_TOT, vr=nSX_SUP)
    call wkvect(p(4), ' V V R ', N_TOT, vr=nSY_SUP)
    call wkvect(p(5), ' V V R ', N_TOT, vr=nS_TOT)
    call wkvect(p(6), ' V V R ', N_TOT, vr=tINF)
    call wkvect(p(7), ' V V R ', N_TOT, vr=tSUP)
    call wkvect(p(8), ' V V R ', N_TOT, vr=ncMAX_INF)
    call wkvect(p(9), ' V V R ', N_TOT, vr=ncMIN_INF)
    call wkvect(p(10), ' V V R ', N_TOT, vr=ncMAX_SUP)
    call wkvect(p(11), ' V V R ', N_TOT, vr=ncMIN_SUP)
    call wkvect(p(12), ' V V R ', N_TOT, vr=RESIDU)
!!!ICI ON A UN PROBLEME DE COHERENCE GLOBALE
!!EN EFFET, DANS LE solver_sandcas1, AngleSUP a N_SUP
!!COMME DIMENSION , ALORS QU'ICI, c'est N_INF!!!
!!==> RESOLVE with variable MONO
    call wkvect(p(13), ' V V R ', N_INF, vr=AngleINF)
    call wkvect(p(14), ' V V R ', N_INF, vr=AngleSUP)

    do i = 1, N_INF

        if (i .eq. 1) then
            AngleINF(i) = -90
        elseif (i .eq. N_INF) then
            AngleINF(i) = 90
        else
            AngleINF(i) = -((N_INF-1)/2-i+1)*thiter
        end if

        theta_inf = AngleINF(i)*r8dgrd()

        !Calc pour tSUPA
        tSUPA = 0
        call solver_sandcas2(2, ySUP, yINF, ht, effrts, fcd, fcd1, cond109, &
                             AngleINF(i), tSUPA, &
                             tINF(i), AngleSUP(i), ResiduA, &
                             nSX_SUP(i), nSY_SUP(i), nSX_INF(i), nSY_INF(i), nS_TOT(i), &
                             ncMAX_SUP(i), ncMIN_SUP(i), ncMAX_INF(i), ncMIN_INF(i))

        !Calc pour tSUPB
        tSUPB = 0.5*ht
        call solver_sandcas2(2, ySUP, yINF, ht, effrts, fcd, fcd1, cond109, &
                             AngleINF(i), tSUPB, &
                             tINF(i), AngleSUP(i), ResiduB, &
                             nSX_SUP(i), nSY_SUP(i), nSX_INF(i), nSY_INF(i), nS_TOT(i), &
                             ncMAX_SUP(i), ncMIN_SUP(i), ncMAX_INF(i), ncMIN_INF(i))

        !Calc pour tSUPmed
        condmed = .false.
        do while (condmed .eqv. (.false.))
            tSUPmed = 0.5*(tSUPA+tSUPB)
            call solver_sandcas2(2, ySUP, yINF, ht, effrts, fcd, fcd1, cond109, &
                                 AngleINF(i), tSUPmed, &
                                 tINF(i), AngleSUP(i), Residumed, &
                                 nSX_SUP(i), nSY_SUP(i), nSX_INF(i), nSY_INF(i), nS_TOT(i), &
                                 ncMAX_SUP(i), ncMIN_SUP(i), ncMAX_INF(i), ncMIN_INF(i))

            !Iteration sur tSUPmed
            if ((RESIDUA .ne. (-1.d0)) .and. (RESIDUB .ne. (-1.d0))) then
                if (RESIDUA*RESIDUB .gt. 0) then
                    if (abs(RESIDUA) .le. abs(RESIDUB)) then
                        tSUPB = tSUPA
                        RESIDUB = RESIDUA
                        tSUPA = max(tSUPA-(tSUPmed-tSUPA), 0.0)
                        !Calc NEW tSUPA
                        call solver_sandcas2(2, ySUP, yINF, ht, effrts, fcd, fcd1, cond109, &
                                             AngleINF(i), tSUPA, &
                                             tINF(i), AngleSUP(i), ResiduA, &
                                             nSX_SUP(i), nSY_SUP(i), nSX_INF(i), &
                                             nSY_INF(i), nS_TOT(i), &
                                             ncMAX_SUP(i), ncMIN_SUP(i), &
                                             ncMAX_INF(i), ncMIN_INF(i))
                    else
                        tSUPA = tSUPB
                        RESIDUA = RESIDUB
                        tSUPB = min(tSUPB+tSUPB-tSUPmed, 0.5/ht)
                        !Calc NEW tSUPB
                        call solver_sandcas2(2, ySUP, yINF, ht, effrts, fcd, fcd1, cond109, &
                                             AngleINF(i), tSUPB, &
                                             tINF(i), AngleSUP(i), ResiduB, &
                                             nSX_SUP(i), nSY_SUP(i), nSX_INF(i), &
                                             nSY_INF(i), nS_TOT(i), &
                                             ncMAX_SUP(i), ncMIN_SUP(i), &
                                             ncMAX_INF(i), ncMIN_INF(i))
                    end if
                else
                    if ((RESIDUA*RESIDUmed) .gt. 0) then
                        tSUPA = tSUPmed
                        RESIDUA = RESIDUmed
                    elseif ((RESIDUB*RESIDUmed) .gt. 0) then
                        tSUPB = tSUPmed
                        RESIDUB = RESIDUmed
                    else
                        ierr = 1
                        goto 99
                    end if
                end if

            elseif (RESIDUA .ne. (-1.d0)) then
                tSUPB = tSUPmed
                RESIDUB = RESIDUmed

            elseif (RESIDUB .ne. (-1.d0)) then
                tSUPA = tSUPmed
                RESIDUB = RESIDUmed
            else

                tINF(i) = -1.d0
                tSUP(i) = -1.d0
                RESIDU(i) = -1.d0
                goto 99
            end if

            if (condmed .eqv. (.false.)) then
                Calc1 = RESIDUmed
                Calc2 = tSUPA-tSUPB
                if ((abs(Calc1) .le. my_epsi) .or. (abs(Calc2) .le. epsilon(Calc2))) then
                    condmed = .true.
                    tSUP(i) = tSUPmed
                    RESIDU(i) = RESIDUmed
                end if
            end if

        end do
        !Pour iteration sur tSUPmed

99      continue
        !Verification de l'aboutissement de l'equilibre
        Calc1 = RESIDU(i)
        if ((tSUP(i) .eq. (-1.d0)) .or. (tINF(i) .eq. (-1.d0)) &
             & .or. ((tSUP(i)+tINF(i)) .gt. ht) .or. (abs(Calc1) .gt. my_epsi)) then
            tSUP(i) = -1.d0
            tINF(i) = -1.d0
            nSX_SUP(i) = -1.d0
            nSX_INF(i) = -1.d0
            nSY_SUP(i) = -1.d0
            nSY_INF(i) = -1.d0
            nS_TOT(i) = -1.d0
            RESIDU(i) = -1.d0
            AngleSUP(i) = -1.d0
            ncMAX_SUP(i) = -1.d0
            ncMIN_SUP(i) = -1.d0
            ncMAX_INF(i) = -1.d0
            ncMIN_INF(i) = -1.d0
        end if

    end do
!iteration sur i=1,N_INF

!First step ok"

!DETERMINATION DES RACINES EVENTUELLEMENT DETECTABLES - ITERATION SUR LES COLONNES, SOIT N_INF

    call solver_sandcas1(N_INF, N_INF, N_TOT, 1, .true., p, &
                         nSX_INF, nSY_INF, nSX_SUP, nSY_SUP, nS_TOT, &
                         tINF, tSUP, ncMAX_INF, ncMIN_INF, ncMAX_SUP, ncMIN_SUP, &
                         RESIDU, AngleINF, AngleSUP)

!---------------------------

!Determination des valeurs optimales

!---------------------------

    call solver_optimum(N_INF, N_SUP, .true., &
                        nSX_INF, nSY_INF, nSX_SUP, nSY_SUP, nS_TOT, &
                        tINF, tSUP, ncMAX_INF, ncMIN_INF, ncMAX_SUP, ncMIN_SUP, &
                        AngleINF, AngleSUP, &
                        fyd, ferrcomp, ferrsyme, slsyme, &
                        vect, ierr)

    dnsxi = vect(1)
    dnsxs = vect(2)
    dnsyi = vect(3)
    dnsys = vect(4)
    etsxi = to_aster_int(vect(5))
    etsxs = to_aster_int(vect(6))
    etsyi = to_aster_int(vect(7))
    etsys = to_aster_int(vect(8))
    snsxi = vect(9)
    snsxs = vect(10)
    snsyi = vect(11)
    snsys = vect(12)
    ncmaxi = vect(13)
    ncmini = vect(14)
    ncmaxs = vect(15)
    ncmins = vect(16)
    t_inf = vect(17)
    t_sup = vect(18)
    theta_inf = vect(19)
    theta_sup = vect(20)

    do i = 1, 15
        call jedetr(p(i))
    end do

end subroutine
