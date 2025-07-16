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

subroutine sandcas4(effrts, ht, enrobi, enrobs, facier, fbeton, gammas, gammac, &
                    thiter, epiter, cond109, ferrcomp, ferrsyme, slsyme, uc, um, &
                    dnsxi, dnsxs, dnsyi, dnsys, etsxi, etsxs, etsyi, etsys, &
                    snsxi, snsxs, snsyi, snsys, ncmaxi, ncmini, ncmaxs, ncmins, &
                    t_inf, t_sup, theta_inf, theta_sup, ierr)

!______________________________________________________________________
!
!      SANDCAS4
!
!      CALCUL DES ACIERS A L'ELU PAR LA METHODE SANDWICH
!      CAS 4 - FERRAILLAGE [+] NON REQUIS EN SUP ET INF
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
!      I EPITER        TAUX D'ITERATION SUR LES EPAISSEURS DES COUCHES
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
!      I UM            UNITE DES DIMENSIONS :
!                         UM = 0 DIMENSIONS EN m
!                         UM = 1 DIMENSIONS EN mm
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
#include "asterc/r8rddg.h"
#include "asterfort/wkvect.h"
#include "asterfort/jedetr.h"
#include "asterfort/mgauss.h"

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
    real(kind=8) :: epiter
    integer(kind=8) :: cond109
    integer(kind=8) :: ferrcomp
    integer(kind=8) :: ferrsyme
    real(kind=8) :: slsyme
    integer(kind=8) :: uc
    integer(kind=8) :: um
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
    real(kind=8) :: pi, fcd, fyd, ySUP, yINF, Z, fcd2, fc
    real(kind=8) :: Nxx, Nxy, Nyy, Mxx, Mxy, Myy
    real(kind=8) :: Ds(6, 6), SOL(6)
    real(kind=8) :: unite_m, det, nctot_opt, Calc, alpha
    integer(kind=8) :: N_IV, N_TOT, i, j, count_sol, iret
    integer(kind=8) :: indx, indx1, indx2, indx3, indx4, indx5
    logical :: cond_inf, cond_sup
    character(20) :: p(10)
    real(kind=8) :: ncX_SUP, ncY_SUP, ncXY_SUP, ncX_INF, ncY_INF, ncXY_INF

    real(kind=8), pointer :: tSUP(:) => null(), tINF(:) => null()
    real(kind=8), pointer :: ncMAX_SUP(:) => null(), ncMIN_SUP(:) => null()
    real(kind=8), pointer :: ncMAX_INF(:) => null(), ncMIN_INF(:) => null()
    real(kind=8), pointer :: RESIDU_INF(:) => null(), RESIDU_SUP(:) => null()
    real(kind=8), pointer :: AngleSUP(:) => null(), AngleINF(:) => null()

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

    Ds(1, 1) = 1
    Ds(1, 2) = 0
    Ds(1, 3) = 0
    Ds(1, 4) = 1
    Ds(1, 5) = 0
    Ds(1, 6) = 0

    Ds(2, 1) = 0
    Ds(2, 2) = 1
    Ds(2, 3) = 0
    Ds(2, 4) = 0
    Ds(2, 5) = 1
    Ds(2, 6) = 0

    Ds(3, 1) = 0
    Ds(3, 2) = 0
    Ds(3, 3) = 1
    Ds(3, 4) = 0
    Ds(3, 5) = 0
    Ds(3, 6) = 1

    Ds(4, 2) = 0
    Ds(4, 3) = 0
    Ds(4, 5) = 0
    Ds(4, 6) = 0

    Ds(5, 1) = 0
    Ds(5, 3) = 0
    Ds(5, 4) = 0
    Ds(5, 6) = 0

    Ds(6, 1) = 0
    Ds(6, 2) = 0
    Ds(6, 4) = 0
    Ds(6, 5) = 0

!-------------------------------------------------------------------------------------------------
!On commence à itérer sur le MODELE SANDWICH A ADOPTER (Determination des epaisseurs des 3 Couches)
!-------------------------------------------------------------------------------------------------

    if (um .eq. 0) then
        unite_m = 1.e3
    Else
        unite_m = 1.
    end if

    N_IV = ceiling(0.5/epiter)+1
    N_TOT = N_IV**2
    do i = 1, 10
        write (p(i), fmt='(A9,I2)') 'SANDCAS4_', i
    end do

    call wkvect(p(1), ' V V R ', N_IV, vr=tINF)
    call wkvect(p(2), ' V V R ', N_IV, vr=tSUP)
    call wkvect(p(3), ' V V R ', N_TOT, vr=ncMAX_INF)
    call wkvect(p(4), ' V V R ', N_TOT, vr=ncMIN_INF)
    call wkvect(p(5), ' V V R ', N_TOT, vr=ncMAX_SUP)
    call wkvect(p(6), ' V V R ', N_TOT, vr=ncMIN_SUP)
    call wkvect(p(7), ' V V R ', N_TOT, vr=RESIDU_INF)
    call wkvect(p(8), ' V V R ', N_TOT, vr=RESIDU_SUP)
    call wkvect(p(9), ' V V R ', N_TOT, vr=AngleINF)
    call wkvect(p(10), ' V V R ', N_TOT, vr=AngleSUP)

    do i = 1, N_IV
        if (i .eq. 1) then
            tSUP(i) = 0.d0
            tINF(i) = 0.d0
        elseif (i .eq. N_IV) then
            tSUP(i) = 0.5*ht
            tINF(i) = 0.5*ht
        else
            tSUP(i) = tSUP(i-1)+epiter*ht
            tINF(i) = tINF(i-1)+epiter*ht
        end if
    end do

    do i = 1, N_IV
    do j = 1, N_IV

        !(i,j) ==> indx
        indx = i+N_IV*(j-1)

        if ((tSUP(i)+tINF(j)) .le. ht) then

            Ds(4, 1) = 0.5*(ht-tSUP(i))
            Ds(4, 4) = -0.5*(ht-tINF(j))
            Ds(5, 2) = 0.5*(ht-tSUP(i))
            Ds(5, 5) = -0.5*(ht-tINF(j))
            Ds(6, 3) = 0.5*(ht-tSUP(i))
            Ds(6, 6) = -0.5*(ht-tINF(j))

            SOL(1) = Nxx
            SOL(2) = Nyy
            SOL(3) = Nxy
            SOL(4) = Mxx
            SOL(5) = Myy
            SOL(6) = Mxy

            call mgauss('NFSP', Ds, SOL, 6, 6, 1, det, iret)

            ncX_SUP = SOL(1)
            ncY_SUP = SOL(2)
            ncXY_SUP = SOL(3)
            ncX_INF = SOL(4)
            ncY_INF = SOL(5)
            ncXY_INF = SOL(6)

            !Calc RESIDU_SUP

            ncMAX_SUP(indx) = 0.5*(ncX_SUP+ncY_SUP)
            ncMAX_SUP(indx) = ncMAX_SUP(indx)-Sqrt(((0.5*ncX_SUP-0.5*ncY_SUP)**2)+(ncXY_SUP**2))
            ncMAX_SUP(indx) = -ncMAX_SUP(indx)

            ncMIN_SUP(indx) = 0.5*(ncX_SUP+ncY_SUP)
            ncMIN_SUP(indx) = ncMIN_SUP(indx)+Sqrt(((0.5*ncX_SUP-0.5*ncY_SUP)**2)+(ncXY_SUP**2))
            ncMIN_SUP(indx) = -ncMIN_SUP(indx)

            if (cond109 .eq. 1) then
                if (abs(ncMAX_SUP(indx)) .gt. epsilon(abs(ncMAX_SUP(indx)))) Then
                    alpha = abs(ncMIN_SUP(indx)/ncMAX_SUP(indx))
                    fcd2 = 0.85*fcd*(1+3.8*alpha)/((1+alpha)**2)
                else
                    fcd2 = fcd
                end if
            else
                fcd2 = fcd
            end if
            fc = fcd2
            RESIDU_SUP(indx) = 0.5*(ncX_SUP+ncY_SUP)
            RESIDU_SUP(indx) = RESIDU_SUP(indx)-Sqrt(((0.5*ncX_SUP-0.5*ncY_SUP)**2)+(ncXY_SUP**2))
            RESIDU_SUP(indx) = RESIDU_SUP(indx)+tSUP(i)*fc

            if (abs(ncXY_SUP) .gt. epsilon(abs(ncXY_SUP))) then
                AngleSUP(indx) = atan((-tSUP(i)*fc-ncX_SUP)/(ncXY_SUP))
            elseif (abs(ncX_SUP) .ge. abs(ncY_SUP)) then
                AngleSUP(indx) = pi/2.
            else
                AngleSUP(indx) = 0
            end if

            if (AngleSUP(indx) .ge. 0) then
                AngleSUP(indx) = 90-r8rddg()*AngleSUP(indx)
            else
                AngleSUP(indx) = -90-r8rddg()*AngleSUP(indx)
            end if

            !Calc RESIDU_INF

            ncMAX_INF(indx) = 0.5*(ncX_INF+ncY_INF)
            ncMAX_INF(indx) = ncMAX_INF(indx)-Sqrt(((0.5*ncX_INF-0.5*ncY_INF)**2)+(ncXY_INF**2))
            ncMAX_INF(indx) = -ncMAX_INF(indx)

            ncMIN_INF(indx) = 0.5*(ncX_INF+ncY_INF)
            ncMIN_INF(indx) = ncMIN_INF(indx)+Sqrt(((0.5*ncX_INF-0.5*ncY_INF)**2)+(ncXY_INF**2))
            ncMIN_INF(indx) = -ncMIN_INF(indx)

            if (cond109 .eq. 1) then
                if (abs(ncMAX_INF(indx)) .gt. epsilon(abs(ncMAX_INF(indx)))) Then
                    alpha = abs(ncMIN_INF(indx)/ncMAX_INF(indx))
                    fcd2 = 0.85*fcd*(1+3.8*alpha)/((1+alpha)**2)
                else
                    fcd2 = fcd
                end if
            else
                fcd2 = fcd
            end if
            fc = fcd2
            RESIDU_INF(indx) = 0.5*(ncX_INF+ncY_INF)
            RESIDU_INF(indx) = RESIDU_INF(indx)-Sqrt(((0.5*ncX_INF-0.5*ncY_INF)**2)+(ncXY_INF**2))
            RESIDU_INF(indx) = RESIDU_INF(indx)+tINF(j)*fc

            if (abs(ncXY_INF) .gt. epsilon(abs(ncXY_INF))) then
                AngleINF(indx) = atan((-tINF(j)*fc-ncX_INF)/(ncXY_INF))
            elseif (abs(ncX_INF) .ge. abs(ncY_INF)) then
                AngleINF(indx) = pi/2.
            else
                AngleINF(indx) = 0
            end if

            if (AngleINF(indx) .ge. 0) then
                AngleINF(indx) = 90-r8rddg()*AngleINF(indx)
            else
                AngleSUP(indx) = -90-r8rddg()*AngleINF(indx)
            end if

        else

            ncMAX_SUP(indx) = -1.d0
            ncMIN_SUP(indx) = -1.d0
            ncMAX_INF(indx) = -1.d0
            ncMIN_INF(indx) = -1.d0
            AngleSUP(indx) = -1.d0
            AngleINF(indx) = -1.d0
            RESIDU_SUP(indx) = -1.d0
            RESIDU_INF(indx) = -1.d0

        end if

    end do
    end do

!---------------------------

!Determination des valeurs optimales

!---------------------------

    count_sol = 0

    do i = 2, (N_IV-1)
    do j = 2, (N_IV-1)

        !(i,j) ==> indx
        indx1 = i+N_IV*(j-1)

        indx2 = (i+1)+N_IV*(j-1)
        indx3 = (i-1)+N_IV*(j-1)

        indx4 = i+N_IV*((j-1)-1)
        indx5 = i+N_IV*((j+1)-1)

        cond_sup = .false.
        cond_inf = .false.

        if ((RESIDU_SUP(indx1)*RESIDU_SUP(indx2) .le. 0) &
            & .or. (RESIDU_SUP(indx1)*RESIDU_SUP(indx3) .le. 0) &
            & .or. (RESIDU_SUP(indx1)*RESIDU_SUP(indx4) .le. 0) &
            & .or. (RESIDU_SUP(indx1)*RESIDU_SUP(indx5) .le. 0)) then
            cond_sup = .true.
        end if
        if ((RESIDU_INF(indx1)*RESIDU_INF(indx2) .le. 0) &
            & .or. (RESIDU_INF(indx1)*RESIDU_INF(indx3) .le. 0) &
            & .or. (RESIDU_INF(indx1)*RESIDU_INF(indx4) .le. 0) &
            & .or. (RESIDU_INF(indx1)*RESIDU_INF(indx5) .le. 0)) then
            cond_inf = .true.
        end if

        if ((cond_sup .eqv. (.true.)) .and. (cond_inf .eqv. (.true.))) then
            count_sol = count_sol+1
            if (count_sol .eq. 1) then
                t_sup = tSUP(i)
                t_inf = tINF(i)
                theta_sup = AngleSUP(indx1)
                theta_inf = AngleINF(indx1)
                ncmaxs = ncMAX_SUP(indx1)
                ncmins = ncMIN_SUP(indx1)
                ncmaxi = ncMAX_INF(indx1)
                ncmini = ncMIN_INF(indx1)
                nctot_opt = abs(ncmaxs)+abs(ncmins)+abs(ncmaxi)+abs(ncmini)
            else
                Calc = abs(ncmaxs)+abs(ncmins)+abs(ncmaxi)+abs(ncmini)
                if (Calc .lt. nctot_opt) then
                    t_sup = tSUP(i)
                    t_inf = tINF(i)
                    theta_sup = AngleSUP(indx1)
                    theta_inf = AngleINF(indx1)
                    ncmaxs = ncMAX_SUP(indx1)
                    ncmins = ncMIN_SUP(indx1)
                    ncmaxi = ncMAX_INF(indx1)
                    ncmini = ncMIN_INF(indx1)
                end if
            end if
        end if
    end do
    end do

    if (count_sol .ne. 0) then
        snsxi = 0.d0
        snsxs = 0.d0
        snsyi = 0.d0
        snsys = 0.d0
        dnsxi = 0.d0
        dnsxs = 0.d0
        dnsyi = 0.d0
        dnsys = 0.d0
        etsxi = 0.d0
        etsxs = 0.d0
        etsyi = 0.d0
        etsys = 0.d0
    else
        ierr = 1
        snsxi = -1.d0
        snsxs = -1.d0
        snsyi = -1.d0
        snsys = -1.d0
        dnsxi = -1.d0
        dnsxs = -1.d0
        dnsyi = -1.d0
        dnsys = -1.d0
        etsxi = -1.d0
        etsxs = -1.d0
        etsyi = -1.d0
        etsys = -1.d0
        t_sup = -1.d0
        t_inf = -1.d0
        theta_sup = -1.d0
        theta_inf = -1.d0
        ncmaxs = -1.d0
        ncmins = -1.d0
        ncmaxi = -1.d0
        ncmini = -1.d0
    end if

    do i = 1, 10
        call jedetr(p(i))
    end do

end subroutine
