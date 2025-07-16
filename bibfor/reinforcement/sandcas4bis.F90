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

subroutine sandcas4bis(effrts, ht, enrobi, enrobs, facier, fbeton, gammas, gammac, &
                       thiter, epiter, aphiter, cond109, ferrcomp, ferrsyme, slsyme, uc, um, &
                       dnsxi, dnsxs, dnsyi, dnsys, etsxi, etsxs, etsyi, etsys, &
                       snsxi, snsxs, snsyi, snsys, ncmaxi, ncmini, ncmaxs, ncmins, &
                       t_inf, t_sup, theta_inf, theta_sup, ierr)

!______________________________________________________________________
!
!      SANDCAS4BIS
!
!      CALCUL DES ACIERS A L'ELU PAR LA METHODE SANDWICH
!      CAS 4 BIS - FERRAILLAGE [+] NON REQUIS EN SUP ET INF
!                  FERRAILLAGE [-] REQUIS EN SUP ET/OU INF
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
!      I APHITER       TAUX D'ITERATION SUR LE RATIO DES CONTRAINTES
!                      PRINCIPALES DE COMPRESSION
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
    real(kind=8) :: aphiter
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
    real(kind=8) :: pi, fcd, fyd, ySUP, yINF, alpha_INF, alpha_SUP, Z
    integer(kind=8) :: i, j, k, l, m, n, N1, N2, N3, N_IIbis, indx, count_F, indx_F, iret
    real(kind=8) :: Nxx, Nxy, Nyy, Mxx, Mxy, Myy, CalcX, CalcY, denum
    real(kind=8) :: unite_m, fcSUP, fcINF, nC_INF, nC_SUP, mC_INF, mC_SUP
    real(kind=8) :: a11, a12, a21, a22, a00, b00, c00, Delta, X_alpha_inf, X_alpha_sup, t1, t2
    logical :: cond_t1, cond_t2
    real(kind=8) :: Ds(4, 4), SOL(4), det, cs, ci, ss, si
    real(kind=8) :: ns_TOT_opt, t_inf_F, t_sup_F, theta_inf_F, theta_sup_F
    real(kind=8), pointer :: nSX_SUP(:) => null(), nSX_INF(:) => null()
    real(kind=8), pointer :: nSY_SUP(:) => null(), nSY_INF(:) => null()
    real(kind=8), pointer :: nS_TOT(:) => null()
    real(kind=8), pointer :: tSUP(:) => null(), tINF(:) => null()
    real(kind=8), pointer :: ncMAX_SUP(:) => null(), ncMIN_SUP(:) => null()
    real(kind=8), pointer :: ncMAX_INF(:) => null(), ncMIN_INF(:) => null()
    real(kind=8), pointer :: AngleSUP(:) => null(), AngleINF(:) => null()
    character(20) :: p(13)

    do i = 1, 13
        write (p(i), fmt='(A12,I2)') 'SANDCAS4BIS_', i
    end do

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

    if (um .eq. 0) then
        unite_m = 1.e3
    else
        unite_m = 1.
    end if

!------------------------------------------------------
!On commence à itérer sur le MODELE SANDWICH A ADOPTER
!------------------------------------------------------

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

    if ((abs(Nxy) .le. epsilon(Nxy)) .and. (abs(Mxy) .le. epsilon(Nxy))) then

        !AngleSUP = 90 / 0 (i=1,N1)
        !AngleINF = 90 / 0 (j=1,N1)
        !tSUP varie de 0 à h/2 (avec un pas de h*epiter) (k=1,N2)
        !tINF varie de 0 à h/2 (avec un pas de h*epiter) (l=1,N2)
        !alpha_SUP varie de 0 à 1 (avec un pas de aphiter) (m=1,N3)
        !alpha_INF varie de 0 à 1 (avec un pas de aphiter) (n=1,N3)

        N1 = 2
        N2 = ceiling(0.5/epiter)+1
        N3 = ceiling(1/aphiter)+1

        N_IIbis = N1*N1*N2*N2*N3*N3

        call wkvect(p(1), ' V V R ', N1, vr=AngleINF)
        call wkvect(p(2), ' V V R ', N1, vr=AngleSUP)
        call wkvect(p(3), ' V V R ', N2, vr=tINF)
        call wkvect(p(4), ' V V R ', N2, vr=tSUP)

        call wkvect(p(5), ' V V R ', N_IIbis, vr=ncMAX_INF)
        call wkvect(p(6), ' V V R ', N_IIbis, vr=ncMIN_INF)
        call wkvect(p(7), ' V V R ', N_IIbis, vr=ncMAX_SUP)
        call wkvect(p(8), ' V V R ', N_IIbis, vr=ncMIN_SUP)

        call wkvect(p(9), ' V V R ', N_IIbis, vr=nsX_INF)
        call wkvect(p(10), ' V V R ', N_IIbis, vr=nsY_INF)
        call wkvect(p(11), ' V V R ', N_IIbis, vr=nsX_SUP)
        call wkvect(p(12), ' V V R ', N_IIbis, vr=nsY_SUP)
        call wkvect(p(13), ' V V R ', N_IIbis, vr=ns_TOT)

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

                        do m = 1, N3

                            if (m .eq. N3) then
                                alpha_SUP = 1.d0
                            else
                                alpha_SUP = aphiter*(m-1)
                            end if

                            do n = 1, N3

                                if (n .eq. N3) then
                                    alpha_INF = 1.d0
                                else
                                    alpha_INF = aphiter*(n-1)
                                end if

                                !(i,j,k,l,m) ==> indx
                                indx = i
                                indx = indx+N1*(j-1)
                                indx = indx+N1*N1*(k-1)
                                indx = indx+N1*N1*N2*(l-1)
                                indx = indx+N1*N1*N2*N2*(m-1)
                                indx = indx+N1*N1*N2*N2*N3*(n-1)

                                if (cond109 .eq. 1) then
                                    fcSUP = 0.85*fcd*(1+3.8*alpha_SUP)/((1+alpha_SUP)**2)
                                    fcINF = 0.85*fcd*(1+3.8*alpha_INF)/((1+alpha_INF)**2)
                                else
                                    fcSUP = fcd
                                    fcINF = fcd
                                end if

                                nC_SUP = -t_sup*fcSUP
                                nC_INF = -t_inf*fcINF
                                mC_SUP = 0.5*(ht-t_sup)*nC_SUP
                                mC_INF = 0.5*(ht-t_inf)*nC_INF

                                ncMAX_INF(indx) = -nC_INF
                                ncMIN_INF(indx) = -alpha_INF*nC_INF
                                ncMAX_SUP(indx) = -nC_SUP
                                ncMIN_SUP(indx) = -alpha_SUP*nC_SUP

                                ci = cos(theta_inf)
                                cs = cos(theta_sup)
                                si = sin(theta_inf)
                                ss = sin(theta_sup)

                                SOL(1) = Nxx-nC_INF*(((si)**2)+alpha_INF*((ci)**2))
                                SOL(1) = SOL(1)-nC_SUP*(((ss)**2)+alpha_SUP*((cs)**2))
                                SOL(2) = Nyy-nC_INF*(((si)**2)+alpha_INF*((ci)**2))
                                SOL(2) = SOL(2)-nC_SUP*(((cs)**2)+alpha_SUP*((ss)**2))
                                SOL(3) = Mxx+mC_INF*(((si)**2)+alpha_INF*((ci)**2))
                                SOL(3) = SOL(3)-mC_SUP*(((ss)**2)+alpha_SUP*((cs)**2))
                                SOL(4) = Myy+mC_INF*(((si)**2)+alpha_INF*((ci)**2))
                                SOL(4) = SOL(4)-mC_SUP*(((cs)**2)+alpha_SUP*((ss)**2))

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
99                              continue
                            end do
                        end do
                    end do
                end do
            end do
        end do

    else

        !AngleSUP varie de -90 à 90 (i=1,N1)
        !AngleINF varie de -90 à 90 (j=1,N1)
        !tSUP varie de 0 à h/2 (avec un pas de h*epiter) (k=1,N2)
        !tINF varie de 0 à h/2 (avec un pas de h*epiter) (l=1,N2)

        N1 = 2*ceiling(90.0/thiter)+1
        N2 = ceiling(0.5/epiter)+1

        N_IIbis = N1*N1*N2*N2

        call wkvect(p(1), ' V V R ', N1, vr=AngleINF)
        call wkvect(p(2), ' V V R ', N1, vr=AngleSUP)
        call wkvect(p(3), ' V V R ', N2, vr=tINF)
        call wkvect(p(4), ' V V R ', N2, vr=tSUP)

        call wkvect(p(5), ' V V R ', N_IIbis, vr=ncMAX_INF)
        call wkvect(p(6), ' V V R ', N_IIbis, vr=ncMIN_INF)
        call wkvect(p(7), ' V V R ', N_IIbis, vr=ncMAX_SUP)
        call wkvect(p(8), ' V V R ', N_IIbis, vr=ncMIN_SUP)

        call wkvect(p(9), ' V V R ', N_IIbis, vr=nsX_INF)
        call wkvect(p(10), ' V V R ', N_IIbis, vr=nsY_INF)
        call wkvect(p(11), ' V V R ', N_IIbis, vr=nsX_SUP)
        call wkvect(p(12), ' V V R ', N_IIbis, vr=nsY_SUP)
        call wkvect(p(13), ' V V R ', N_IIbis, vr=ns_TOT)

        do i = 1, N1

            if (i .eq. 1) then
                AngleSUP(i) = -90
            elseif (i .eq. N1) then
                AngleSUP(i) = 90
            else
                AngleSUP(i) = -((N1-1)/2-i+1)*thiter
            end if
            theta_sup = AngleSUP(i)*r8dgrd()

            do j = 1, N1

                if (j .eq. 1) then
                    AngleINF(j) = -90
                elseif (j .eq. N1) then
                    AngleINF(j) = 90
                else
                    AngleINF(j) = -((N1-1)/2-j+1)*thiter
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

                        !Solve 1st equation
                        a11 = -fcd*t_sup*cos(theta_sup)*sin(theta_sup)
                        a12 = -fcd*t_inf*cos(theta_inf)*sin(theta_inf)
                        a21 = 0.5*a11
                        a22 = 0.5*a12

                        if (a12 .le. epsilon(a12)) then
                            alpha_INF = -1.d0
                            alpha_SUP = -1.d0
                        else
                            denum = a21+a11*a22/a12
                            if (abs(denum) .lt. epsilon(denum)) then
                                goto 999
                            end if
                            X_alpha_sup = (Mxy+(a22/a12)*Nxy)/denum
                            X_alpha_inf = (Nxy-a11*X_alpha_sup)/a12
                            if (cond109 .eq. 1) then
                                !Solve alpha_SUP
                                a00 = 3.23+X_alpha_sup
                                b00 = 2*X_alpha_sup-2.38
                                c00 = X_alpha_sup-0.85
                                Delta = b00**2-4*a00*c00
                                if (a00 .le. epsilon(a00)) then
                                    if (b00 .gt. epsilon(b00)) then
                                        alpha_SUP = -c00/b00
                                    else
                                        alpha_SUP = -1.d0
                                    end if
                                elseif (Delta .ge. 0) then
                                    t1 = (-b00-sqrt(Delta))/(2*a00)
                                    t2 = (-b00+sqrt(Delta))/(2*a00)
                                    cond_t1 = .false.
                                    cond_t2 = .false.
                                    if ((t1 .ge. 0) .and. (t1 .le. 1.0)) then
                                        cond_t1 = .true.
                                    elseif ((t2 .ge. 0) .and. (t2 .le. 1.0)) then
                                        cond_t2 = .true.
                                    end if
                                    if ((cond_t1 .eqv. (.true.)) &
                                         &.and. (cond_t2 .eqv. (.true.))) then
                                        alpha_SUP = min(t1, t2)
                                    elseif (cond_t1 .eqv. (.true.)) then
                                        alpha_SUP = t1
                                    elseif (cond_t2 .eqv. (.true.)) then
                                        alpha_SUP = t2
                                    else
                                        alpha_SUP = -1.d0
                                    end if
                                else
                                    alpha_SUP = -1.d0
                                end if
                                !Solve alpha_INF
                                a00 = 3.23+X_alpha_inf
                                b00 = 2*X_alpha_inf-2.38
                                c00 = X_alpha_inf-0.85
                                Delta = b00**2-4*a00*c00
                                if (a00 .le. epsilon(a00)) then
                                    if (b00 .gt. epsilon(b00)) then
                                        alpha_INF = -c00/b00
                                    else
                                        alpha_INF = -1.d0
                                    end if
                                elseif (Delta .ge. 0) then
                                    t1 = (-b00-sqrt(Delta))/(2*a00)
                                    t2 = (-b00+sqrt(Delta))/(2*a00)
                                    cond_t1 = .false.
                                    cond_t2 = .false.
                                    if ((t1 .ge. 0) .and. (t1 .le. 1.0)) then
                                        cond_t1 = .true.
                                    elseif ((t2 .ge. 0) .and. (t2 .le. 1.0)) then
                                        cond_t2 = .true.
                                    end if
                                    if ((cond_t1 .eqv. (.true.)) &
                                         &.and. (cond_t2 .eqv. (.true.))) then
                                        alpha_INF = min(t1, t2)
                                    elseif (cond_t1 .eqv. (.true.)) then
                                        alpha_INF = t1
                                    elseif (cond_t2 .eqv. (.true.)) then
                                        alpha_INF = t2
                                    else
                                        alpha_INF = -1.d0
                                    end if
                                else
                                    alpha_INF = -1.d0
                                end if
                            else
                                !Solve alpha_SUP
                                alpha_SUP = 1-X_alpha_sup
                                if ((alpha_SUP .lt. 0) .or. (alpha_SUP .gt. 1.0)) then
                                    alpha_SUP = -1.d0
                                end if
                                !Solve alpha_INF
                                alpha_INF = 1-X_alpha_inf
                                if ((alpha_INF .lt. 0) .or. (alpha_INF .gt. 1.0)) then
                                    alpha_INF = -1.d0
                                end if
                            end if
                        end if

                        if ((alpha_INF .ne. (-1.d0)) .and. (alpha_SUP .ne. (-1.d0))) then
                            if (cond109 .eq. 1) then
                                fcSUP = 0.85*fcd*(1+3.8*alpha_SUP)/((1+alpha_SUP)**2)
                                fcINF = 0.85*fcd*(1+3.8*alpha_INF)/((1+alpha_INF)**2)
                            else
                                fcSUP = fcd
                                fcINF = fcd
                            end if

                            nC_SUP = -t_sup*fcSUP
                            nC_INF = -t_inf*fcINF
                            mC_SUP = 0.5*(ht-t_sup)*nC_SUP
                            mC_INF = 0.5*(ht-t_inf)*nC_INF

                            ncMAX_INF(indx) = -nC_INF
                            ncMIN_INF(indx) = -alpha_INF*nC_INF
                            ncMAX_SUP(indx) = -nC_SUP
                            ncMIN_SUP(indx) = -alpha_SUP*nC_SUP

                            ci = cos(theta_inf)
                            cs = cos(theta_sup)
                            si = sin(theta_inf)
                            ss = sin(theta_sup)

                            SOL(1) = Nxx-nC_INF*(((si)**2)+alpha_INF*((ci)**2))
                            SOL(1) = SOL(1)-nC_SUP*(((ss)**2)+alpha_SUP*((cs)**2))
                            SOL(2) = Nyy-nC_INF*(((si)**2)+alpha_INF*((ci)**2))
                            SOL(2) = SOL(2)-nC_SUP*(((cs)**2)+alpha_SUP*((ss)**2))
                            SOL(3) = Mxx+mC_INF*(((si)**2)+alpha_INF*((ci)**2))
                            SOL(3) = SOL(3)-mC_SUP*(((ss)**2)+alpha_SUP*((cs)**2))
                            SOL(4) = Myy+mC_INF*(((si)**2)+alpha_INF*((ci)**2))
                            SOL(4) = SOL(4)-mC_SUP*(((cs)**2)+alpha_SUP*((ss)**2))

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
                                        goto 999
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
                        end if
999                     continue
                    end do
                end do
            end do
        end do
    end if

!Choix de la solution finale
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

    do i = 1, 13
        call jedetr(p(i))
    end do

end subroutine
