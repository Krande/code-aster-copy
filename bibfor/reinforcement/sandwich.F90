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

subroutine sandwich(enrobi, enrobs, facier, fbeton, gammas, gammac, &
                    thiter, epiter, aphiter, cond109, &
                    ferrcomp, ferrsyme, slsyme, &
                    epucisa, ferrmin, rholmin, rhotmin, compress, &
                    alphacc, eys, typdiag, clacier, &
                    uc, um, &
                    ht, effrts, dnsits, ierrl, ierrt)

!______________________________________________________________________
!
!      SANDWICH
!
!      CALCUL DES ACIERS A L'ELU PAR LA METHODE SANDWICH
!
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
!      I EPUCISA       IMPACT DE L'EFFORT TRANCHANT ET DE LA TORSION SUR LE
!                      FERRAILLAGE LONGITUDINAL?
!                      (0 = NON, 1 = OUI)
!      I FERRMIN       PRISE EN COMPTE DU FERRA MINI (0 = NON, 1 = OUI, 2 = CODE)
!      I RHOLMIN       RATIO DE FERRAILLAGE LONGI MINI (A RENSEIGNER SI FERMIN='OUI')
!      I RHOTMIN       RATIO DE FERRAILLAGE TRNSV MINI (A RENSEIGNER SI FERMIN='OUI')
!      I COMPRESS      VALORISATION DE LA COMPRESSION POUR LES ACIERS TRANSVERSAUX
!                      (0 = NON, 1 = OUI)
!      I ALPHACC       COEFFICIENT DE SECURITE SUR LA RESISTANCE
!                      DE CALCUL DU BETON EN COMPRESSION
!      I EYS           MODULE D'YOUNG DE L'ACIER
!      I TYPDIAG       TYPE DE DIAGRAMME UTILISÉ POUR L'ACIER
!                            TYPDIAG = 1 ("B1" ==> PALIER INCLINÉ)
!                            TYPDIAG = 2 ("B2" ==> PALIER HORIZONTAL)
!      I CLACIER       CLASSE DE DUCTILITE DES ACIERS (UTILISE POUR EC2) :
!                            CLACIER = 0 ACIER PEU DUCTILE (CLASSE A)
!                            CLACIER = 1 ACIER MOYENNEMENT DUCTILE (CLASSE B)
!                            CLACIER = 2 ACIER FORTEMENT DUCTILE (CLASSE C)
!      I UC            UNITE DES CONTRAINTES :
!                         UC = 0 CONTRAINTES EN Pa
!                         UC = 1 CONTRAINTES EN MPa
!      I UM            UNITE DES DIMENSIONS :
!                         UM = 0 DIMENSIONS EN m
!                         UM = 1 DIMENSIONS EN mm
!      I HT            EPAISSEUR DE LA COQUE
!      I EFFRTS        (DIM 8) TORSEUR DES EFFORTS, MOMENTS, ...
!                         EFFRTS(1) = NXX
!                         EFFRTS(2) = NYY
!                         EFFRTS(3) = NXY
!                         EFFRTS(4) = MXX
!                         EFFRTS(5) = MYY
!                         EFFRTS(6) = MXY
!                         EFFRTS(7) = QX
!                         EFFRTS(8) = QY
!
!      O DNSITS        (DIM 6) DENSITES
!                            1..4 : SURFACES D'ACIER LONGITUDINAL
!                            5..6 : TRANSVERSAL
!      O IERRL     CODE RETOUR LONGI (0 = OK)
!      O IERRT     CODE RETOUR TRNSV (0 = OK)
!
!______________________________________________________________________
!
    implicit none
#include "asterc/r8pi.h"
#include "asterc/r8dgrd.h"
#include "asterfort/sandcas1.h"
#include "asterfort/sandcas2.h"
#include "asterfort/sandcas3.h"
#include "asterfort/sandcas4.h"
#include "asterfort/sandcas2bis.h"
#include "asterfort/sandcas3bis.h"
#include "asterfort/sandcas4bis.h"
#include "asterfort/clcplq.h"

!VARIABLES PRINCIPALES
!------------------------------
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
    integer(kind=8) :: epucisa
    integer(kind=8) :: ferrmin
    real(kind=8) :: rholmin
    real(kind=8) :: rhotmin
    integer(kind=8) :: compress
    real(kind=8) :: alphacc
    real(kind=8) :: eys
    integer(kind=8) :: typdiag
    integer(kind=8) :: clacier
    integer(kind=8) :: uc
    integer(kind=8) :: um
    real(kind=8) :: ht
    real(kind=8) :: effrts(8)
    real(kind=8) :: dnsits(6)
    integer(kind=8) :: ierrl
    integer(kind=8) :: ierrt

!Variables de calcul
    real(kind=8) :: fcd, fyd, ySUP, yINF, Z, pi, zI, zS, denom, d, fctm
    real(kind=8) :: Dnsx, Dnsxy, Dnsy, Dnsx_NEW, Dnsxy_NEW, Dnsy_NEW
    logical :: PREMIERE_ITERATION
    real(kind=8) :: Nxx, Nyy, Nxy, Mxx, Myy, Mxy, Qx, Qy
    real(kind=8) :: Nx_SUP, Nx_INF, Ny_SUP, Ny_INF, Nxy_SUP, Nxy_INF
    real(kind=8) :: cond_trac_inf, cond_trac_sup
    integer(kind=8) :: CAS_SUP, CAS_INF, COUNT_ITER, j, ierr, nb, precs
    real(kind=8) :: dnsxi, dnsxs, dnsyi, dnsys, dnsxt, dnsyt
    integer(kind=8) :: etsxi, etsxs, etsyi, etsys
    real(kind=8) :: snsxi, snsxs, snsyi, snsys
    real(kind=8) :: ncmaxi, ncmini, ncmaxs, ncmins
    real(kind=8) :: t_inf, t_sup, theta_inf, theta_sup, unite_pa, unite_m
    real(kind=8) :: rhosxi, rhosxs, rhosyi, rhosys
    real(kind=8) :: betha, rhoCALC, VEd, NEd, Scp, alphaCW, Nu, Nu1, ratio
    real(kind=8) :: tC, zMOY, CRdc, kBAR, rho, k1, vmin, vCALC, VRdc, AsT, ThetaB, VRdmax
    real(kind=8) :: thetaB_ITER(233), EQ_ITER(233), AsT_ITER(233)

    pi = r8pi()
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
    Qx = effrts(7)
    Qy = effrts(8)

    dnsits(1) = -1.d0
    dnsits(2) = -1.d0
    dnsits(3) = -1.d0
    dnsits(4) = -1.d0
    dnsits(5) = -1.d0
    dnsits(6) = -1.d0

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

    Dnsx = 0
    Dnsxy = 0
    Dnsy = 0
    Dnsx_NEW = 0
    Dnsxy_NEW = 0
    Dnsy_NEW = 0
    PREMIERE_ITERATION = .true.
    ierrl = 0
    ierrt = 0

    do while ((PREMIERE_ITERATION .eqv. (.true.)) .or. (Dnsx_NEW .gt. Dnsx) &
               & .or. (Dnsy_NEW .gt. Dnsy) .or. (Dnsxy_NEW .gt. Dnsxy))

        Dnsx = Dnsx_NEW
        Dnsy = Dnsy_NEW
        Dnsxy = Dnsxy_NEW

        Nx_SUP = Nxx*((Z-ySUP)/Z)+Mxx/Z+0.5*Dnsx
        Nx_INF = Nxx*((Z-yINF)/Z)-Mxx/Z+0.5*Dnsx
        Ny_SUP = Nyy*((Z-ySUP)/Z)+Myy/Z+0.5*Dnsy
        Ny_INF = Nyy*((Z-yINF)/Z)-Myy/Z+0.5*Dnsy
        Nxy_SUP = Nxy*((Z-ySUP)/Z)+Mxy/Z+0.5*Dnsxy
        Nxy_INF = Nxy*((Z-yINF)/Z)-Mxy/Z+0.5*Dnsxy

        Nxx = Nxx+Dnsx
        Nyy = Nyy+Dnsy
        Nxy = Nxy+Dnsxy
        effrts(1) = Nxx
        effrts(2) = Nyy
        effrts(3) = Nxy

   !!Cas de la nappe SUP
        if ((Nx_SUP .ge. (-abs(Nxy_SUP))) .and. (Ny_SUP .ge. (-abs(Nxy_SUP)))) then
            CAS_SUP = 1
        elseif (Nx_SUP .lt. (-abs(Nxy_SUP))) then
            if (Ny_SUP .gt. ((Nxy_SUP*Nxy_SUP)/Nx_SUP)) then
                CAS_SUP = 2
            else
                CAS_SUP = 4
            end if
        elseif (Ny_SUP .lt. (-abs(Nxy_SUP))) then
            if (Nx_SUP .gt. ((Nxy_SUP*Nxy_SUP)/Ny_SUP)) then
                CAS_SUP = 3
            else
                CAS_SUP = 4
            end if
        else
            CAS_SUP = 4
        end if

   !!Cas de la nappe INF
        if ((Nx_INF .ge. (-abs(Nxy_INF))) .and. (Ny_INF .ge. (-abs(Nxy_INF)))) then
            CAS_INF = 1
        elseif (Nx_INF .lt. (-abs(Nxy_INF))) then
            if (Ny_INF .gt. ((Nxy_INF*Nxy_INF)/Nx_INF)) then
                CAS_INF = 2
            else
                CAS_INF = 4
            end if
        elseif (Ny_INF .lt. (-abs(Nxy_INF))) then
            if (Nx_INF .gt. ((Nxy_INF*Nxy_INF)/Ny_INF)) then
                CAS_INF = 3
            else
                CAS_INF = 4
            end if
        else
            CAS_INF = 4
        end if

        !On commence à itérer sur le MODELE SANDWICH A ADOPTER
        !------------------------------------------------------

        !CAS I - REINF NEEDED FOR SUP AND INF
        !------------------------------------
        if ((CAS_SUP .ne. 4) .and. (CAS_INF .ne. 4)) then

            call sandcas1(effrts(1:6), ht, enrobi, enrobs, facier, fbeton, gammas, gammac, &
                          thiter, epiter, cond109, ferrcomp, ferrsyme, slsyme, uc, &
                          dnsxi, dnsxs, dnsyi, dnsys, etsxi, etsxs, etsyi, etsys, &
                          snsxi, snsxs, snsyi, snsys, ncmaxi, ncmini, ncmaxs, ncmins, &
                          t_inf, t_sup, theta_inf, theta_sup, ierr)

            !CAS II - REINF NEEDED FOR SUP ONLY
            !---------------------------------------
        elseif ((CAS_INF .eq. 4) .and. (CAS_SUP .ne. 4)) then

            call sandcas2(effrts(1:6), ht, enrobi, enrobs, facier, fbeton, gammas, gammac, &
                          thiter, cond109, ferrcomp, ferrsyme, slsyme, uc, &
                          dnsxi, dnsxs, dnsyi, dnsys, etsxi, etsxs, etsyi, etsys, &
                          snsxi, snsxs, snsyi, snsys, ncmaxi, ncmini, ncmaxs, ncmins, &
                          t_inf, t_sup, theta_inf, theta_sup, ierr)

            if ((ierr .eq. 1) .and. (ferrcomp .eq. 1)) then

                call sandcas2bis(effrts(1:6), ht, enrobi, enrobs, facier, fbeton, gammas, gammac, &
                                 thiter, epiter, aphiter, cond109, ferrcomp, ferrsyme, &
                                 slsyme, uc, um, &
                                 dnsxi, dnsxs, dnsyi, dnsys, etsxi, etsxs, etsyi, etsys, &
                                 snsxi, snsxs, snsyi, snsys, ncmaxi, ncmini, ncmaxs, ncmins, &
                                 t_inf, t_sup, theta_inf, theta_sup, ierr)

            end if

            !CAS III - REINF NEEDED FOR INF ONLY
            !---------------------------------------
        elseif ((CAS_INF .ne. 4) .and. (CAS_SUP .eq. 4)) then

            call sandcas3(effrts(1:6), ht, enrobi, enrobs, facier, fbeton, gammas, gammac, &
                          thiter, cond109, ferrcomp, ferrsyme, slsyme, uc, &
                          dnsxi, dnsxs, dnsyi, dnsys, etsxi, etsxs, etsyi, etsys, &
                          snsxi, snsxs, snsyi, snsys, ncmaxi, ncmini, ncmaxs, ncmins, &
                          t_inf, t_sup, theta_inf, theta_sup, ierr)

            if ((ierr .eq. 1) .and. (ferrcomp .eq. 1)) then

                call sandcas3bis(effrts(1:6), ht, enrobi, enrobs, facier, fbeton, gammas, gammac, &
                                 thiter, epiter, aphiter, cond109, ferrcomp, ferrsyme, &
                                 slsyme, uc, um, &
                                 dnsxi, dnsxs, dnsyi, dnsys, etsxi, etsxs, etsyi, etsys, &
                                 snsxi, snsxs, snsyi, snsys, ncmaxi, ncmini, ncmaxs, ncmins, &
                                 t_inf, t_sup, theta_inf, theta_sup, ierr)

            end if

            !CAS IV - NO REINF NEEDED
            !------------------------
        else

            call sandcas4(effrts(1:6), ht, enrobi, enrobs, facier, fbeton, gammas, gammac, &
                          thiter, epiter, cond109, ferrcomp, ferrsyme, slsyme, uc, um, &
                          dnsxi, dnsxs, dnsyi, dnsys, etsxi, etsxs, etsyi, etsys, &
                          snsxi, snsxs, snsyi, snsys, ncmaxi, ncmini, ncmaxs, ncmins, &
                          t_inf, t_sup, theta_inf, theta_sup, ierr)

            if ((ierr .eq. 1) .and. (ferrcomp .eq. 1)) then

                call sandcas4bis(effrts(1:6), ht, enrobi, enrobs, facier, fbeton, gammas, gammac, &
                                 thiter, epiter, aphiter, cond109, ferrcomp, ferrsyme, &
                                 slsyme, uc, um, &
                                 dnsxi, dnsxs, dnsyi, dnsys, etsxi, etsxs, etsyi, etsys, &
                                 snsxi, snsxs, snsyi, snsys, ncmaxi, ncmini, ncmaxs, ncmins, &
                                 t_inf, t_sup, theta_inf, theta_sup, ierr)

            end if

        end if
        !Distinction entre CAS_SUP et CAS_INF

!Ajout de l'impact du cisaillement HORS-PLAN
!-------------------------------------------------------------------------------------------------

        ierrl = ierr

        if (ierrl .gt. 0) then

            !NON CONVERGENCE DETECTED => CALL CAPRA MAURY
            nb = ceiling(180/thiter)
            precs = ceiling(1/epiter)
            call clcplq(0, 2, nb, precs, &
                        ferrsyme, slsyme, ferrcomp, epucisa, &
                        ferrmin, rholmin, rhotmin, compress, -1.d0, &
                        enrobi, enrobs, -1.d0, -1.d0, -1.d0, &
                        alphacc, gammas, gammac, facier, eys, typdiag, &
                        fbeton, clacier, uc, um, &
                        -1.d0, -1.d0, -1.d0, -1.d0, -1.d0, -1.d0, -1.d0, -1.d0, &
                        ht, effrts, dnsits, ierrl, ierrt)
            ierrl = 2001
            goto 998

        end if

        if (ierrl .gt. 0) then
            dnsxi = 0
            dnsxs = 0
            dnsyi = 0
            dnsys = 0
            etsxi = 0
            etsxs = 0
            etsyi = 0
            etsys = 0
            t_inf = 0
            t_sup = 0
        end if

        cond_trac_sup = 0.d0
        cond_trac_inf = 0.d0

        if (etsxi .eq. 0) then
            rhosxi = dnsxi
            cond_trac_inf = 1.d0
        else
            rhosxi = 0.d0
        end if
        if (etsyi .eq. 0) then
            rhosyi = dnsyi
            cond_trac_inf = 1.d0
        else
            rhosyi = 0.d0
        end if
        if (etsxs .eq. 0) then
            rhosxs = dnsxs
            cond_trac_sup = 1.d0
        else
            rhosxs = 0.d0
        end if
        if (etsys .eq. 0) then
            rhosys = dnsys
            cond_trac_sup = 1.d0
        else
            rhosys = 0.d0
        end if

!Definition de d
        if ((cond_trac_sup .eq. (1.d0)) .and. (cond_trac_inf .eq. (1.d0))) then
            d = ht-max(enrobi, enrobs)
        elseif (cond_trac_sup .eq. (1.d0)) then
            d = ht-enrobs
        elseif (cond_trac_inf .eq. (1.d0)) then
            d = ht-enrobi
        else
            d = ht
        end if

        if (Qx .ne. 0) then
            betha = atan(Qy/Qx)
        else
            betha = pi/2
        end if

        rhoCALC = (rhosxi+rhosxs)*Cos(betha)*Cos(betha)+(rhosyi+rhosys)*Sin(betha)*Sin(betha)
        rhoCALC = rhoCALC/d
        VEd = Sqrt(Qx*Qx+Qy*Qy)
        NEd = Nxx*((Cos(betha))**2)+Nyy*((Sin(betha))**2)+2*Nxy*(Cos(betha))*(Sin(betha))
        Scp = Max(-NEd/ht, 0.0)
        if ((Scp .lt. 0) .or. (compress .eq. 0)) then
            Scp = 0
        end if
        if (Scp .gt. 0) then
            Scp = min(Scp, 0.2*fbeton/gammac)
        end if

!Calcul de alphaCW
        if (Scp .le. 0) then
            alphaCW = 1.d0
        elseif (Scp .le. (0.25*fcd)) then
            alphaCW = 1.d0+Scp/fcd
        elseif (Scp .le. (0.5*fcd)) then
            alphaCW = 1.25
        elseif (Scp .lt. fcd) then
            alphaCW = 2.5*(1-Scp/fcd)
        else
            !A VOIR!
            alphaCW = 1.d0
        end if

        Nu = 0.6*(1-fbeton*unite_pa/250.d0)
        if ((fbeton*unite_pa) .le. 60) Then
            Nu1 = 0.6
        else
            Nu1 = 0.9-fbeton*unite_pa/200.d0
        end if

!Calcul du bras de levier des efforts internes
        tC = ht-0.5*t_inf-0.5*t_sup
        zS = (0.5*t_sup)*(1.d0-cond_trac_sup)*(t_sup*fcd)+(t_sup-enrobs)*(dnsxs+dnsys)*fyd
        denom = (1.d0-cond_trac_sup)*(t_sup*fcd)+(dnsxs+dnsys)*fyd
        if (abs(denom) .gt. epsilon(denom)) then
            zS = zS/denom
        else
            zS = 0
        end if

        zI = (0.5*t_inf)*(1.d0-cond_trac_inf)*(t_inf*fcd)+(t_inf-enrobi)*(dnsxi+dnsyi)*fyd
        denom = (1.d0-cond_trac_inf)*(t_inf*fcd)+(dnsxi+dnsyi)*fyd
        if (abs(denom) .gt. epsilon(denom)) then
            zI = zI/denom
        else
            zI = 0
        end if

        zMOY = tC+zS+zI

!Calcul de la resistance du Beton SEUL!

        CRdc = 0.18/gammac
        kBAR = 1.d0+(200.d0/(unite_m*d))**0.5
        kBAR = Min(kBAR, 2.0)
        rho = rhoCALC
        k1 = 0.15
        vmin = (0.34/gammac)*((fbeton*unite_pa)**0.5)
        vmin = vmin/unite_pa

        if (rho .gt. 0) then
            vCALC = CRdc*kBAR*((100*rho*(fbeton*unite_pa))**(1./3.))
            vCALC = vCALC/unite_pa
        else
            vCALC = 0.d0
        end if

        VRdc = (Max(vCALC, vmin)+k1*Scp)*tC

!Calcul de AsT!

        if (abs(VRdc) .lt. epsilon(VRdc)) then
            ratio = 1.01
        else
            ratio = VEd/VRdc
        end if

        if (ratio .le. 1) then

            AsT = 0.d0
            ThetaB = -1.d0
            VRdmax = -1.d0

        else

            do j = 1, 233
                thetaB_ITER(j) = 21.8+(j-1)*0.1
                thetaB = ThetaB_ITER(j)*r8dgrd()
                VRdmax = alphaCW*zMOY*Nu1*fcd/(Tan(ThetaB)+1/(Tan(ThetaB)))
                EQ_ITER(j) = VRdmax-VEd
                AsT_ITER(j) = (VEd*Tan(ThetaB))/(zMOY*fyd)
            end do

            COUNT_ITER = 0
            AsT = -1.d0
            ThetaB = -1.d0
            VRdmax = -1.d0

            do j = 1, 233
            if (EQ_ITER(j) .ge. 0) then
                COUNT_ITER = COUNT_ITER+1
                if (COUNT_ITER .eq. 1) then
                    ThetaB = ThetaB_ITER(j)*r8dgrd()
                    AsT = AsT_ITER(j)
                    VRdmax = EQ_ITER(j)+VEd
                else
                    if (AsT_ITER(j) .lt. AsT) then
                        AsT = AsT_ITER(j)
                        ThetaB = ThetaB_ITER(j)*r8dgrd()
                        VRdmax = EQ_ITER(j)+VEd
                    end if
                end if
            end if
            end do

        end if

        if (AsT .eq. 0) then

            dnsxt = 0.d0
            dnsyt = 0.d0
            goto 999

        elseif (AsT .eq. (-1)) then

            dnsxt = -1.d0
            dnsyt = -1.d0
            ierrt = 4
            goto 999

        else

            Dnsx_NEW = (Qx**2)/(VEd*Tan(ThetaB))
            Dnsy_NEW = (Qy**2)/(VEd*Tan(ThetaB))
            Dnsxy_NEW = (Qx*Qy)/(VEd*Tan(ThetaB))
            dnsxt = AsT*Cos(betha)*Cos(betha)
            dnsyt = AsT*Sin(betha)*Sin(betha)
            if (ierrl .gt. 0) then
                goto 999
            end if

        end if

        if (epucisa .eq. 0) then
            goto 999
        end if

        if (PREMIERE_ITERATION .eqv. (.true.)) then
            PREMIERE_ITERATION = .false.
        end if

    end do

999 continue

!  -- VERIFICATION DU FERRAILLAGE MINIMUM :
!  ----------------------------------------

    if ((ferrmin .eq. 1) .or. (ferrmin .eq. 2)) then
        if (fbeton .le. (50*unite_pa)) then
            fctm = 0.30*((fbeton/unite_pa)**(2.0/3.0))
        else
            fctm = 2.12*LOG(1.0+((fbeton/unite_pa)+8.0)/10.0)
        end if

        if (ferrmin .eq. 2) then
            rholmin = max(0.26*(fctm/facier), 0.0013)
            rhotmin = 0
        end if

        !ferraillage longitudinal
        if (ierrl .eq. 0) then
            d = ht-enrobi
            if ((dnsxi+dnsyi) .lt. (rholmin*d)) then
                dnsxi = 0.5*rholmin*d
                dnsyi = 0.5*rholmin*d
            end if
            d = ht-enrobs
            if ((dnsxs+dnsys) .lt. (rholmin*d)) then
                dnsxs = 0.5*rholmin*d
                dnsys = 0.5*rholmin*d
            end if
        end if

        if (ierrt .eq. 0) then
            if ((dnsxt+dnsyt) .lt. (rhotmin*ht)) then
                dnsxt = 0.5*rhotmin*ht
                dnsyt = 0.5*rhotmin*ht
            end if
        end if
    end if

!  -- GESTION DES MESSAGES D'ERREURS ET LIEN AVEC TE0146 :
!  -------------------------------------------------------
    if (ierrl .eq. 1) then
!      Equilibre Sandwich non possible
!      Essayer de changer les critères de précision
        ierrl = 2001
    end if

    if (ierrl .eq. 2) then
!      Forte Compression !
!      Alarme dans te0146 + on sort de la boucle + densité = -1 pour l'élément
        ierrl = 2002
    end if

    if (ierrl .eq. 3) then
!      Ferraillage symétrique non possible!
!      Alarme dans te0146 + on sort de la boucle + densité = -1 pour l'élément
        ierrl = 2003
    end if

    if (ierrt .eq. 4) then
!      Béton trop cisaillé !
!      Alarme dans te0146 + on sort de la boucle + dnstra = -1 pour l'élément
        ierrt = 2004
    end if

!   FER2_R =  DNSXI DNSXS DNSYI DNSYS DNSXT DNSYT DNSVOL CONSTRUC
!               1     2     3     4     5     6      7       8

    if (ierrl .eq. 0) then
        dnsits(1) = dnsxi
        dnsits(2) = dnsyi
        dnsits(3) = dnsxs
        dnsits(4) = dnsys
    end if

    if (ierrt .eq. 0) then
        dnsits(5) = dnsxt
        dnsits(6) = dnsyt
    end if

998 continue

end subroutine
