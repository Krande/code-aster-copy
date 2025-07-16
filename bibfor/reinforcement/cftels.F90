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

subroutine cftels(typco, typstru, effrts, effm, effn, efft, effmt, &
                  dnsinf, dnssup, &
                  sigmsi, sigmss, sigmci, sigmcs, alpha, &
                  ht, bw, enrobi, enrobs, facier, fbeton, &
                  scmaxi, scmaxs, ssmax, uc, um, &
                  compress, dnstra, thetab, ak, uk, ierr)
!______________________________________________________________________
!
!     CFTELS
!
!      CALCUL DU FERRAILLAGE TRANSVERSAL A L'ELS
!
!      I TYPCO       CODIFICATION UTILISEE (1 = BAEL91, 2 = EC2)
!      I TYPSTRU     TYPE DE STRUCTURE : 0 = 2D, 1 = 1D
!      I EFFRTS      (DIM 8) TORSEUR DES EFFORTS, MOMENTS, ...
!      I EFFM        MOMENT DE FLEXION
!      I EFFN        EFFORT NORMAL
!      I EFFT        EFFORT TRANCHANT
!      I EFFMT       MOMENT DE TORSION
!      I DNSINF      DENSITE DE L'ACIER LONGITUDINAL INFERIEUR
!      I DNSSUP      DENSITE DE L'ACIER LONGITUDINAL SUPERIEUR
!      I SIGMSI      CONTRAINTE AU NIVEAU DE L'ACIER INFERIEUR
!      I SIGMSS      CONTRAINTE AU NIVEAU DE L'ACIER SUPERIEUR
!      I SIGMCI      CONTRAINTE AU NIVEAU DE LA FIBRE SUPERIEURE DE BETON
!      I SIGMCS      CONTRAINTE AU NIVEAU DE LA FIBRE INFERIEURE DE BETON
!      I ALPHA       COEFFICIENT DE PROFONDEUR DE L'AN
!      I HT          EPAISSEUR DE LA SECTION
!      I BW          LARGEUR DE LA SECTION
!      I ENROBI      ENROBAGE DES ARMATURES INFERIEURES
!      I ENROBS      ENROBAGE DES ARMATURES SUPERIEURES
!      I FACIER      LIMITE D'ELASTICITE DES ACIERS (CONTRAINTE)
!      I FBETON      RESISTANCE EN COMPRESSION DU BETON (CONTRAINTE)
!      I SCMAXI      CONTRAINTE DE COMPRESSION MAXI DU BETON EN FIBRE INF
!      I SCMAXS      CONTRAINTE DE COMPRESSION MAXI DU BETON EN FIBRE SUP
!      I SSMAX       CONTRAINTE MAXI DE L'ACIER DE FLEXION
!      I UC          UNITE DES CONTRAINTES :
!                        UC = 0 CONTRAINTES EN Pa
!                        UC = 1 CONTRAINTES EN MPa
!      I UM          UNITE DES DIMENSIONS :
!                        UM = 0 DIMENSIONS EN m
!                        UM = 1 DIMENSIONS EN mm
!      I COMPRESS    PRISE EN COMPTE DE LA COMPRESSION
!                        COMPRESS = 0 NON
!                        COMPRESS = 1 OUI
!
!      O DNSTRA      DENSITE DE FERRAILLAGE TRANSVERSAL
!      O THETAB      ANGLE D'INCLINAISON DES BIELLES DE COMPRESSION
!      O AK          AIRE INTERIEURE AU FEUILLET MOYEN DE
!                    RESISTANCE EN TORSION
!      O UK          PERIMETRE DE L'AIRE AK
!      O IERR        CODE RETOUR (0 = OK)
!
!______________________________________________________________________
!
!
    implicit none
!
!
#include "asterc/r8pi.h"
#include "asterc/r8dgrd.h"
!
!
    integer(kind=8) ::typco
    integer(kind=8) :: typstru
    real(kind=8) :: effrts(8)
    real(kind=8) :: effm
    real(kind=8) :: effn
    real(kind=8) :: efft
    real(kind=8) :: effmt
    real(kind=8) :: dnsinf
    real(kind=8) :: dnssup
    real(kind=8) :: sigmsi
    real(kind=8) :: sigmss
    real(kind=8) :: sigmci
    real(kind=8) :: sigmcs
    real(kind=8) :: alpha
    real(kind=8) :: ht
    real(kind=8) :: bw
    real(kind=8) :: enrobi
    real(kind=8) :: enrobs
    real(kind=8) :: facier
    real(kind=8) :: fbeton
    real(kind=8) :: scmaxi
    real(kind=8) :: scmaxs
    real(kind=8) :: ssmax
    integer(kind=8) :: uc
    integer(kind=8) :: um
    integer(kind=8) :: compress
    real(kind=8) :: dnstra
    real(kind=8) :: thetab
    real(kind=8) :: ak
    real(kind=8) :: uk
    integer(kind=8) :: ierr

!-----------------------------------------------------------------------
!!!!VARIABLES DE CALCUL
!-----------------------------------------------------------------------

    real(kind=8) :: d, d0, tk, fctm, fctd, fcd, fyd, TRdmax, VRdmax
    real(kind=8) :: sigmat, unite_pa, unite_m, VEd, TEd, VEdT, VRdc, TRdc, Scp
    real(kind=8) :: alphaCW, z1, z2, zMOY
    real(kind=8) :: Nu, Nu1, kBAR, vmin, CRdc, k1, rhoL, vCALC
    real(kind=8), dimension(0:232) :: thetab_ITER, eq_ITER, dnstra_ITER
    real(kind=8) :: gammac, gammas, denom, pi
    integer(kind=8) :: countV, j
!------------------------------------------------------------------------

    pi = r8pi()
    if (effm .ge. 0.) then
        d = ht-enrobi
        d0 = enrobs
    else
        d = ht-enrobs
        d0 = enrobi
    end if

    fcd = 0.5*(scmaxi+scmaxs)
    fyd = ssmax
    gammac = fbeton/fcd
    gammas = facier/fyd
    VEd = Abs(efft)
    TEd = Abs(effmt)

!INITIALISATION DU CODE RETOUR
    ierr = 0

!IMPACT DE LA TORSION SUR CISAILLEMENT
    tk = Max((bw*ht/(2*(bw+ht))), 2*(ht-d), 2*d0)
    ak = (bw-tk)*(ht-tk)
    uk = 2*(bw-tk+ht-tk)
    VEdT = (TEd/(2*ak))*(ht-tk)

!CALCUL POUR CODIFICATION = BAEL91

    if (typco .eq. 1) then

        if (typstru .eq. 0) then
            sigmat = sqrt(effrts(7)*effrts(7)+effrts(8)*effrts(8))/d
        elseif (typstru .eq. 1) then
            sigmat = (VEd+VEdT)/(bw*d)
        end if
        dnstra = sigmat/fyd
        thetab = 45.0*r8dgrd()

!CALCUL POUR CODIFICATION = EC2

    elseif (typco .eq. 2) then

!   CALCULS INTERMEDIAIRES
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

        if ((fbeton*unite_pa) .le. 50) then
            fctm = 0.7*0.3*((fbeton*unite_pa)**(2./3.))
        else
            fctm = 0.7*2.12*log(1+0.1*(fbeton*unite_pa+8))
        end if
        fctd = fctm/(gammac*unite_pa)

!   Ajout Effet de la Torsion

        tk = Max((bw*ht/(2*(bw+ht))), 2*(ht-d), 2*d0)
        ak = (bw-tk)*(ht-tk)
        uk = 2*(bw-tk+ht-tk)

        TRdc = fctd*tk*2*Ak

        Scp = effn/(bw*ht)
        if ((Scp .lt. 0) .or. (compress .eq. 0)) then
            Scp = 0
        end if
        if (Scp .gt. 0) then
            Scp = min(Scp, 0.2*fbeton/gammac)
        end if

!   Impact de la torsion sur Cisaillement

        VEdT = (TEd/(2*ak))*(ht-tk)

!   Calcul de alphaCW
        if (compress .eq. 0) then
            alphaCW = 1
        elseif (compress .eq. 1) then
            if (Scp .le. 0) then
                alphaCW = 1
            elseif (Scp .le. (0.25*fbeton)) then
                alphaCW = 1+Scp/fbeton
            elseif (Scp .le. (0.5*fbeton)) then
                alphaCW = 1.25
            elseif (Scp .lt. fbeton) then
                alphaCW = 2.5*(1-Scp/fbeton)
            else
                alphaCW = 0
            end if
        end if

        Nu = 0.6*(1-fbeton*unite_pa/250.d0)
        if ((fbeton*unite_pa) .le. 60) Then
            Nu1 = 0.6
        else
            Nu1 = 0.9-fbeton*unite_pa/200.d0
        end if

!   Calcul du bras de levier des efforts internes

        if ((alpha .ge. 0) .and. (alpha .le. 1) .and. (effm .gt. 0)) then
            z1 = d-d0
            z2 = (1-alpha/3.0)*d
            zMOY = (dnssup*sigmss)*z1+(0.5*sigmcs*alpha*d*bw)*z2
            denom = dnssup*sigmss+0.5*sigmcs*alpha*d*bw
            !if (denom.ne.0) then
            if (abs(denom) .gt. epsilon(denom)) then
                zMOY = zMOY/denom
            else
                zMOY = 0.9*d
            end if
        elseif ((alpha .ge. 0) .and. (alpha .le. 1) .and. (effm .lt. 0)) then
            z1 = d-d0
            z2 = (1-alpha/3.0)*d
            zMOY = (dnsinf*sigmsi)*z1+(0.5*sigmci*alpha*d*bw)*z2
            denom = dnsinf*sigmsi+0.5*sigmci*alpha*d*bw
            !if (denom.ne.0) then
            if (abs(denom) .gt. epsilon(denom)) then
                zMOY = zMOY/denom
            else
                zMOY = 0.9*d
            end if
        else
            zMOY = 0.9*d
        end if

!   Calcul de la resistance du Beton SEUL

        kBAR = 1.+(200.d0/max(unite_m*d, unite_m*(ht-d0)))**0.5
        kBAR = min(kBAR, 2.0)
        if (typstru .eq. 0) then
            vmin = (0.34/gammac)*((fbeton*unite_pa)**0.5)
        elseif (typstru .eq. 1) then
            vmin = (0.053/gammac)*(kBAR**1.5)*((fbeton*unite_pa)**0.5)
        end if
        vmin = vmin/unite_pa
        CRdc = 0.18/gammac
        k1 = 0.15
        rhoL = 0
        if (sigmss .lt. 0) then
            rhoL = rhoL+dnssup
        end if
        if (sigmsi .lt. 0) then
            rhoL = rhoL+dnsinf
        end if
        rhoL = rhoL/(bw*d)
        if (rhoL .gt. 0) then
            vCALC = CRdc*kBAR*((100*rhoL*(fbeton*unite_pa))**(1./3.))
            vCALC = vCALC/unite_pa
        else
            vCALC = 0
        end if

        VRdc = (Max(vCALC, vmin)+k1*Scp)*bw*d

!   Calcul de DNSTRA

        if ((TEd/TRdc+VEd/VRdc) .le. 1) Then

            dnstra = 0
            thetab = -1

        else

            do j = 0, 232
                thetab_ITER(j) = 21.8+j*0.1
                thetab = Thetab_ITER(j)*r8dgrd()
                VRdmax = alphaCW*bw*zMOY*Nu1*fcd/(tan(thetab)+1.0/(tan(thetab)))
                TRdmax = 2*Nu*alphaCW*fcd*Ak*tk*sin(thetab)*cos(thetab)
                eq_ITER(j) = VEd/VRdmax+TEd/TRdmax
                dnstra_ITER(j) = (VEd+2*VEdT)*tan(thetab)/(zMOY*fyd)
            end do

            countV = 0
            dnstra = -1
            thetab = -1

            do j = 0, 232
                if (eq_ITER(j) .le. 1) then
                    countV = countV+1
                    if (countV .eq. 1) then
                        thetab = thetab_ITER(j)*r8dgrd()
                        dnstra = dnstra_ITER(j)
                    else
                        if (dnstra_ITER(j) .lt. dnstra) then
                            thetab = thetab_ITER(j)*r8dgrd()
                            dnstra = dnstra_ITER(j)
                        end if
                    end if
                end if
            end do

            if (countV .eq. 0) then
                ierr = 1
            end if

        end if

    end if

end subroutine
