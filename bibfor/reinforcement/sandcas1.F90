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

subroutine sandcas1(effrts, ht, enrobi, enrobs, facier, fbeton, gammas, gammac, &
                    thiter, epiter, cond109, ferrcomp, ferrsyme, slsyme, uc, &
                    dnsxi, dnsxs, dnsyi, dnsys, etsxi, etsxs, etsyi, etsys, &
                    snsxi, snsxs, snsyi, snsys, ncmaxi, ncmini, ncmaxs, ncmins, &
                    t_inf, t_sup, theta_inf, theta_sup, ierr)

!______________________________________________________________________
!
!      SANDCAS1
!
!      CALCUL DES ACIERS A L'ELU PAR LA METHODE SANDWICH
!      CAS 1 - FERRAILLAGE [+] REQUIS EN SUP ET INF
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
    real(kind=8) :: fcd, fyd, ySUP, yINF, Z, fcd1, fc, pi, vect(20)
    real(kind=8) :: Nxx, Nyy, Nxy, Mxx, Myy, Mxy
    integer(kind=8) :: N_SUP, N_INF, N_TOT, i, j, indx
    integer(kind=8) :: k, N1, N2, l
    integer(kind=8) :: iret
    real(kind=8) :: nS_TOT_opt, Calc, CalcX, CalcY, det
    real(kind=8) :: unite_pa, a00, b00, c00, DELTA, xA, xB
    real(kind=8) :: nC_INF, mC_INF, nC_SUP, mC_SUP, Ds(4, 4), SOL(4)
    real(kind=8) :: Calc1, Calc2, Calc3, Calc4
    logical :: COND_xA, COND_xB
    character(20) :: p(15)
    real(kind=8), pointer :: nSX_SUP(:) => null(), nSX_INF(:) => null()
    real(kind=8), pointer :: nSY_SUP(:) => null(), nSY_INF(:) => null()
    real(kind=8), pointer :: nS_TOT(:) => null()
    real(kind=8), pointer :: tSUP(:) => null(), tINF(:) => null()
    real(kind=8), pointer :: ncMAX_SUP(:) => null(), ncMIN_SUP(:) => null()
    real(kind=8), pointer :: ncMAX_INF(:) => null(), ncMIN_INF(:) => null()
    real(kind=8), pointer :: RESIDU(:) => null()
    real(kind=8), pointer :: AngleSUP(:) => null(), AngleINF(:) => null()
    integer(kind=8) :: indx_F, count_F
    real(kind=8) :: theta_sup_F, theta_inf_F, t_sup_F, t_inf_F, alpha_INF
    real(kind=8) :: ci, cs, ss, si

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
    Ds(1, 2) = 1
    Ds(1, 3) = 0
    Ds(1, 4) = 0
    Ds(2, 1) = 0
    Ds(2, 2) = 0
    Ds(2, 3) = 1
    Ds(2, 4) = 1
    Ds(3, 1) = ySUP
    Ds(3, 2) = -yINF
    Ds(3, 3) = 0
    Ds(3, 4) = 0
    Ds(4, 1) = 0
    Ds(4, 2) = 0
    Ds(4, 3) = ySUP
    Ds(4, 4) = -yINF
    count_F = 0
    indx_F = 0
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

    fc = fcd1

    do i = 1, 15
        write (p(i), fmt='(A9,I2)') 'SANDCAS1_', i
    end do

!-------------------------------------------------------------------------------------------------
!On commence à itérer sur le MODELE SANDWICH A ADOPTER (Determination des epaisseurs des 3 Couches)
!-------------------------------------------------------------------------------------------------

    if ((abs(Nxy) .le. epsilon(Nxy)) .and. (abs(Mxy) .le. epsilon(Mxy))) then

        N1 = 2
        N2 = ceiling(0.5/epiter)+1
        N_TOT = N1*N1*N2*N2

        call wkvect(p(1), ' V V R ', N1, vr=AngleINF)
        call wkvect(p(2), ' V V R ', N1, vr=AngleSUP)
        call wkvect(p(3), ' V V R ', N2, vr=tINF)
        call wkvect(p(4), ' V V R ', N2, vr=tSUP)

        call wkvect(p(5), ' V V R ', N_TOT, vr=ncMAX_INF)
        call wkvect(p(6), ' V V R ', N_TOT, vr=ncMIN_INF)
        call wkvect(p(7), ' V V R ', N_TOT, vr=ncMAX_SUP)
        call wkvect(p(8), ' V V R ', N_TOT, vr=ncMIN_SUP)

        call wkvect(p(9), ' V V R ', N_TOT, vr=nsX_INF)
        call wkvect(p(10), ' V V R ', N_TOT, vr=nsY_INF)
        call wkvect(p(11), ' V V R ', N_TOT, vr=nsX_SUP)
        call wkvect(p(12), ' V V R ', N_TOT, vr=nsY_SUP)
        call wkvect(p(13), ' V V R ', N_TOT, vr=ns_TOT)

        call wkvect(p(14), ' V V R ', N_TOT, vr=RESIDU)

        do i = 1, N1

            if (i .eq. 1) then
                AngleSUP(i) = 90.0
            else
                AngleSUP(i) = 0.0
            end if
            theta_sup = AngleSUP(i)*r8dgrd()

            do j = 1, N1

                if (j .eq. 1) then
                    AngleINF(j) = 90.0
                else
                    AngleINF(j) = 0.0
                end if
                theta_inf = AngleINF(j)*r8dgrd()

                do k = 1, N2

                    if (k .eq. 1) then
                        tSUP(k) = 0.d0
                    elseif (k .eq. N2) then
                        tSUP(k) = 0.5*ht
                    else
                        tSUP(k) = tSUP(k-1)+epiter*ht
                    end if
                    t_sup = tSUP(k)

                    do l = 1, N2

                        if (l .eq. 1) then
                            tINF(l) = 0.d0
                        elseif (l .eq. N2) then
                            tINF(l) = 0.5*ht
                        else
                            tINF(l) = tINF(l-1)+epiter*ht
                        end if
                        t_inf = tINF(l)

                        indx = i+N1*(j-1)+N1*N1*(k-1)+N1*N1*N2*(l-1)

                        nC_SUP = -t_sup*fc
                        nC_INF = -t_inf*fc
                        mC_SUP = 0.5*(ht-t_sup)*nC_SUP
                        mC_INF = 0.5*(ht-t_inf)*nC_INF

                        ncMAX_INF(indx) = -nC_INF
                        ncMIN_INF(indx) = 0.0
                        ncMAX_SUP(indx) = -nC_SUP
                        ncMIN_SUP(indx) = 0.0

                        alpha_INF = 0.d0

                        ci = cos(theta_inf)
                        cs = cos(theta_sup)
                        si = sin(theta_inf)
                        ss = sin(theta_sup)
                        SOL(1) = Nxx-nC_SUP*((ss)**2)
                        SOL(1) = SOL(1)-nC_INF*(((si)**2)+alpha_INF*((ci)**2))
                        SOL(2) = Nyy-nC_SUP*((cs)**2)
                        SOL(2) = SOL(2)-nC_INF*(((ci)**2)+alpha_INF*((si)**2))
                        SOL(3) = Mxx-mC_SUP*((ss)**2)
                        SOL(3) = SOL(3)+mC_INF*(((si)**2)+alpha_INF*((ci)**2))
                        SOL(4) = Myy-mC_SUP*((cs)**2)
                        SOL(4) = SOL(4)+mC_INF*(((ci)**2)+alpha_INF*((si)**2))

                        call mgauss('NFSP', Ds, SOL, 4, 4, 1, det, iret)

                        nSX_SUP(indx) = SOL(1)
                        nSX_INF(indx) = SOL(2)
                        nSY_SUP(indx) = SOL(3)
                        nSY_INF(indx) = SOL(4)
                        ns_TOT(indx) = abs(nSX_SUP(indx))+abs(nSX_INF(indx))
                        ns_TOT(indx) = ns_TOT(indx)+abs(nSY_SUP(indx))+abs(nSY_INF(indx))

                        if ((t_inf+t_sup) .le. ht) then
                            CalcX = abs(abs(nSX_SUP(indx))-abs(nSX_INF(indx)))
                            CalcY = abs(abs(nSY_SUP(indx))-abs(nSY_INF(indx)))
                            if (ferrsyme .eq. 1) then
                                if ((CalcX .gt. slsyme) .or. (CalcY .gt. slsyme)) then
                                    ierr = 2
                                    goto 99
                                end if
                            end if
                            count_F = count_F+1
                            if ((count_F .eq. 1) .or. (ns_TOT(indx) .lt. ns_TOT_opt)) then
                                indx_F = indx
                                ns_TOT_opt = ns_TOT(indx)
                                theta_sup_F = theta_sup
                                theta_inf_F = theta_inf
                                t_sup_F = t_sup
                                t_inf_F = t_inf
                            end if
                        end if
99                      continue
                    end do
                end do
            end do
        end do

        if (count_F .gt. 0) then

            ierr = 0
            t_sup = t_sup_F
            t_inf = t_inf_F
            theta_sup = r8rddg()*theta_sup_F
            theta_inf = r8rddg()*theta_inf_F
            dnsxi = nSX_INF(indx_F)/fyd
            dnsxs = nSX_SUP(indx_F)/fyd
            dnsyi = nSY_INF(indx_F)/fyd
            dnsys = nSY_SUP(indx_F)/fyd
            if (dnsxi .ge. 0) then
                etsxi = 0
            else
                etsxi = 1
            end if
            if (dnsxs .ge. 0) then
                etsxs = 0
            else
                etsxs = 1
            end if
            if (dnsyi .ge. 0) then
                etsyi = 0
            else
                etsyi = 1
            end if
            if (dnsys .ge. 0) then
                etsys = 0
            else
                etsys = 1
            end if

            dnsxi = abs(dnsxi)
            dnsxs = abs(dnsxs)
            dnsyi = abs(dnsyi)
            dnsys = abs(dnsys)

            if (abs(t_sup) .ge. epsilon(t_sup)) then
                ncmaxs = ncMAX_SUP(indx_F)/t_sup
                ncmins = ncMIN_SUP(indx_F)/t_sup
            else
                ncmaxs = 0.d0
                ncmins = 0.d0
            end if
            if (abs(t_inf) .ge. epsilon(t_inf)) then
                ncmaxi = ncMAX_INF(indx_F)/t_inf
                ncmini = ncMIN_INF(indx_F)/t_inf
            else
                ncmaxi = 0.d0
                ncmini = 0.d0
            end if

            if (abs(dnsxs) .ge. epsilon(dnsxs)) then
                snsxs = fyd
            else
                snsxs = 0.d0
            end if
            if (abs(dnsys) .ge. epsilon(dnsys)) then
                snsys = fyd
            else
                snsys = 0.d0
            end if
            if (abs(dnsxi) .ge. epsilon(dnsxi)) then
                snsxi = fyd
            else
                snsxi = 0.d0
            end if
            if (abs(dnsyi) .ge. epsilon(dnsyi)) then
                snsyi = fyd
            else
                snsyi = 0.d0
            end if

        else

            dnsxs = -1.d0
            dnsys = -1.d0
            dnsxi = -1.d0
            dnsyi = -1.d0

            etsxs = 2
            etsys = 2
            etsxi = 2
            etsyi = 2

            t_sup = -1.d0
            t_inf = -1.d0
            theta_sup = -1.d0
            theta_inf = -1.d0

            ncmaxs = -1.d0
            ncmins = -1.d0
            ncmaxi = -1.d0
            ncmini = -1.d0

            snsxs = -1.d0
            snsys = -1.d0
            snsxi = -1.d0
            snsyi = -1.d0

            if (ierr .eq. 0) then
                ierr = 1
            end if

        end if

    else

        N_SUP = 2*ceiling(90.0/thiter)+1
        N_INF = N_SUP
        N_TOT = N_SUP*N_INF

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
        call wkvect(p(13), ' V V R ', N_INF, vr=AngleINF)
        call wkvect(p(14), ' V V R ', N_SUP, vr=AngleSUP)

        do i = 1, N_SUP
        do j = 1, N_INF

            !(i,j) ==> indx
            indx = i+N_SUP*(j-1)
            RESIDU(indx) = 0.d0

            if (i .eq. 1) then
                AngleSUP(i) = -90
            elseif (i .eq. N_SUP) then
                AngleSUP(i) = 90
            else
                AngleSUP(i) = -((N_SUP-1)/2-i+1)*thiter
            end if
            if (j .eq. 1) then
                AngleINF(j) = -90
            elseif (j .eq. N_INF) then
                AngleINF(j) = 90
            else
                AngleINF(j) = -((N_INF-1)/2-j+1)*thiter
            end if

            theta_sup = AngleSUP(i)*r8dgrd()
            theta_inf = AngleINF(j)*r8dgrd()

            !Distinction sur ThetaSUP et ThetaINF (cas TRIVIAUX ou non)
            !1ere Configuration
            Calc1 = abs(AngleSUP(i))
            Calc2 = abs(AngleINF(j))
            Calc3 = abs(abs(AngleSUP(i))-90)
            Calc4 = abs(abs(AngleINF(j))-90)

            if ((Calc1 .gt. epsilon(Calc1)) .and. (Calc2 .gt. epsilon(Calc2)) .and. &
                 & (Calc3 .gt. epsilon(Calc3)) .and. (Calc4 .gt. epsilon(Calc4))) then

                si = Sin(theta_inf)
                ci = Cos(theta_inf)
                ss = Sin(theta_sup)
                cs = Cos(theta_sup)
                a00 = fc*ss*cs-fc*si*ci*(((ss*cs)/(si*ci))**2)
                b00 = -2*ht*fc*ss*cs-2*Nxy*((ss*cs)/(si*ci))
                c00 = -ht*Nxy-2*Mxy-Nxy*Nxy/(fc*si*ci)

                if (abs(a00) .gt. epsilon(a00)) then

                    DELTA = b00*b00-4*a00*c00

                    if (DELTA .gt. epsilon(DELTA)) then
                        xA = (-b00+sqrt(DELTA))/(2*a00)
                        xB = (-b00-sqrt(DELTA))/(2*a00)
                        if ((xA .ge. 0) .and. (xA .le. (0.5*ht))) then
                            COND_xA = .true.
                        else
                            COND_xA = .false.
                        end if
                        if ((xB .ge. 0) .and. (xB .le. (0.5*ht))) then
                            COND_xB = .true.
                        else
                            COND_xB = .false.
                        end if
                    else
                        COND_xA = .false.
                        COND_xB = .false.
                    end if

                    if (COND_xA .eqv. (.true.) .and. (COND_xB .eqv. (.true.))) then
                        tSUP(indx) = Min(xA, xB)
                    elseif (COND_xA .eqv. (.true.)) then
                        tSUP(indx) = xA
                    elseif (COND_xB .eqv. (.true.)) then
                        tSUP(indx) = xB
                    else
                        tSUP(indx) = -1.d0
                    end if

                else

                    if (abs(b00) .gt. epsilon(b00)) then
                        tSUP(indx) = -c00/b00
                        if ((tSUP(indx) .lt. 0) .or. (tSUP(indx) .gt. ht)) then
                            tSUP(indx) = -1.d0
                        end if
                    else
                        tSUP(indx) = -1.d0
                    end if

                end if

                if (tSUP(indx) .ne. (-1.d0)) then
                    tINF(indx) = -((ss*cs)/(si*ci))*tSUP(indx)
                    tINF(indx) = tINF(indx)-Nxy/(fc*si*ci)
                    if ((tINF(indx) .lt. 0) .or. (tINF(indx) .gt. ht)) then
                        tINF(indx) = -1.d0
                    end if
                else
                    tINF(indx) = -1.d0
                end if

                !2eme Configuration
            elseif (((Calc1 .lt. epsilon(Calc1)) .or. (Calc3 .lt. epsilon(Calc3))) &
                     & .and. (Calc2 .ge. epsilon(Calc2)) .and. (Calc4 .ge. epsilon(Calc4))) then

                tSUP(indx) = 0
                if (abs(Nxy) .ge. epsilon(Nxy)) then
                    tINF(indx) = -Nxy/(fc*Sin(theta_inf)*Cos(theta_inf))
                    nC_INF = -tINF(indx)*fc
                    mC_INF = 0.5*(ht-tINF(indx))*nC_INF
                    Calc = Mxy+mC_INF*Sin(theta_inf)*Cos(theta_inf)
                    if (abs(Calc) .gt. epsilon(Calc)) then
                        tSUP(indx) = -1.d0
                        tINF(indx) = -1.d0
                    elseif ((tINF(indx) .lt. 0) .or. (tINF(indx) .gt. (0.5*ht))) then
                        tSUP(indx) = -1.d0
                        tINF(indx) = -1.d0
                    end if
                else
                    if (abs(Mxy) .lt. epsilon(Mxy)) then
                        tINF(indx) = 0
                    else
                        tSUP(indx) = -1.d0
                        tINF(indx) = -1.d0
                    end if
                end if

                !3eme Configuration
            elseif (((Calc2 .lt. epsilon(Calc2)) .or. (Calc4 .lt. epsilon(Calc4))) &
                     & .and. (Calc1 .ge. epsilon(Calc1)) .and. (Calc3 .ge. epsilon(Calc3))) then

                tINF(indx) = 0
                if (abs(Nxy) .ge. epsilon(Nxy)) then
                    tSUP(indx) = -Nxy/(fc*Sin(theta_sup)*Cos(theta_sup))
                    nC_SUP = -tSUP(indx)*fc
                    mC_SUP = 0.5*(ht-tSUP(indx))*nC_SUP
                    Calc = Mxy-mC_SUP*Sin(theta_sup)*Cos(theta_sup)
                    if (abs(Calc) .ge. epsilon(Calc)) then
                        tSUP(indx) = -1.d0
                        tINF(indx) = -1.d0
                    elseif ((tSUP(indx) .lt. 0) .or. (tSUP(indx) .gt. (0.5*ht))) then
                        tSUP(indx) = -1.d0
                        tINF(indx) = -1.d0
                    end if
                else
                    if (abs(Mxy) .lt. epsilon(Mxy)) then
                        tSUP(indx) = 0
                    else
                        tSUP(indx) = -1.d0
                        tINF(indx) = -1.d0
                    end if
                end if

                !4eme Configuration
            elseif (((Calc1 .lt. epsilon(Calc1)) .or. (Calc3 .lt. epsilon(Calc3))) &
                     & .and. ((Calc2 .lt. epsilon(Calc2)) .or. (Calc4 .lt. epsilon(Calc4)))) then

                if ((abs(Nxy) .lt. epsilon(Nxy)) .and. (abs(Mxy) .lt. epsilon(Mxy))) then
                    tSUP(indx) = 0
                    tINF(indx) = 0
                else
                    tSUP(indx) = -1.d0
                    tINF(indx) = -1.d0
                end if

            end if

            !Determination des ferraillage en fonction des epaisseurs calculées
            if ((tSUP(indx) .ne. (-1.d0)) .and. (tINF(indx) .ne. (-1.d0)) &
                 & .and. (tSUP(indx)+tINF(indx) .le. ht)) then

                nC_SUP = -tSUP(indx)*fc
                nC_INF = -tINF(indx)*fc
                mC_SUP = 0.5*(ht-tSUP(indx))*nC_SUP
                mC_INF = 0.5*(ht-tINF(indx))*nC_INF

                ncMAX_SUP(indx) = -nC_SUP
                ncMIN_SUP(indx) = 0
                ncMAX_INF(indx) = -nC_INF
                ncMIN_INF(indx) = 0

                SOL(1) = Nxx-nC_SUP*((Sin(theta_sup))**2)-nC_INF*((Sin(theta_inf))**2)
                SOL(2) = Nyy-nC_SUP*((Cos(theta_sup))**2)-nC_INF*((Cos(theta_inf))**2)
                SOL(3) = Mxx-mC_SUP*((Sin(theta_sup))**2)+mC_INF*((Sin(theta_inf))**2)
                SOL(4) = Myy-mC_SUP*((Cos(theta_sup))**2)+mC_INF*((Cos(theta_inf))**2)

                call mgauss('NFSP', Ds, SOL, 4, 4, 1, det, iret)
                nSX_SUP(indx) = SOL(1)
                nSX_INF(indx) = SOL(2)
                nSY_SUP(indx) = SOL(3)
                nSY_INF(indx) = SOL(4)
                nS_TOT(indx) = Abs(nSX_SUP(indx))
                nS_TOT(indx) = nS_TOT(indx)+Abs(nSX_INF(indx))
                nS_TOT(indx) = nS_TOT(indx)+Abs(nSY_SUP(indx))
                nS_TOT(indx) = nS_TOT(indx)+Abs(nSY_INF(indx))

            else

                nSX_SUP(indx) = -1.d0
                nSX_INF(indx) = -1.d0
                nSY_SUP(indx) = -1.d0
                nSY_INF(indx) = -1.d0
                nS_TOT(indx) = -1.d0
                ncMAX_SUP(indx) = -1.d0
                ncMIN_SUP(indx) = -1.d0
                ncMAX_INF(indx) = -1.d0
                ncMIN_INF(indx) = -1.d0

            end if

        end do
        end do

!DETERMINATION DES RACINES EVENTUELLEMENT DETECTABLES - ITERATION SUR LES COLONNES, SOIT N_INF

        call solver_sandcas1(N_INF, N_SUP, N_TOT, 1, .false., p, &
                             nSX_INF, nSY_INF, nSX_SUP, nSY_SUP, nS_TOT, &
                             tINF, tSUP, ncMAX_INF, ncMIN_INF, ncMAX_SUP, ncMIN_INF, &
                             RESIDU, AngleINF, AngleSUP)

!DETERMINATION DES RACINES EVENTUELLEMENT DETECTABLES - ITERATION SUR LES LIGNES, SOIT N_SUP

        call solver_sandcas1(N_INF, N_SUP, N_TOT, 2, .false., p, &
                             nSX_INF, nSY_INF, nSX_SUP, nSY_SUP, nS_TOT, &
                             tINF, tSUP, ncMAX_INF, ncMIN_INF, ncMAX_SUP, ncMIN_INF, &
                             RESIDU, AngleINF, AngleSUP)

!---------------------------

!Determination des valeurs optimales

!---------------------------

        call solver_optimum(N_INF, N_SUP, .false., &
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

    end if

    do i = 1, 14
        call jedetr(p(i))
    end do

end subroutine
