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

subroutine glbelsqp(typco, cequi, effrts, ht, bw, &
                    enrobyi, enrobys, enrobzi, enrobzs, &
                    facier, fbeton, sigelsqp, kt, eys, &
                    wmaxyi, wmaxys, wmaxzi, wmaxzs, &
                    phiyi, phiys, phizi, phizs, &
                    precs, ferrsyme, slsyme, ferrcomp, &
                    epucisa, ferrmin, rholmin, rhotmin, compress, uc, um, &
                    dnsits, ierrl, ierrt)
!______________________________________________________________________
!
!      GLBELSQP

!      CALCUL GLOBAL DU FERRAILLAGE DES POUTRES A L'ELS QP
!
!      I TYPCO     CODIFICATION UTILISEE (1 = BAEL91, 2 = EC2)
!      I CEQUI     COEFFICIENT D'EQUIVALENCE ACIER/BETON
!      I EFFM      MOMENT DE FLEXION
!      I EFFN      EFFORT NORMAL
!      I HT        HAUTEUR DE LA SECTION
!      I BW        LARGEUR DE LA SECTION
!      I ENROBYI   ENROBAGE DES ARMATURES INF SUIVANT L'AXE Y
!      I ENROBYS   ENROBAGE DES ARMATURES SUP SUIVANT L'AXE Y
!      I ENROBZI   ENROBAGE DES ARMATURES INF SUIVANT L'AXE Z
!      I ENROBZS   ENROBAGE DES ARMATURES SUP SUIVANT L'AXE Z
!      I FACIER    LIMITE D'ELASTICITE DES ACIERS (CONTRAINTE)
!      I FBETON    RESISTANCE EN COMPRESSION DU BETON (CONTRAINTE)
!      I SIGELSQP  CONTRAINTE ADMISSIBLE DANS LE BETON À L'ELS QP
!      I KT        COEFFICIENT DE DURÉE DE CHARGEMENT
!      I EYS       MODULE D'YOUNG DE L'ACIER
!      I WMAXYI    OUVERTURE MAXIMALE DES FISSURES
!                     EN FACE INFÉRIEURE SUIVANT L'AXE Y (1D)
!      I WMAXYS    OUVERTURE MAXIMALE DES FISSURES
!                     EN FACE SUPÉRIEURE SUIVANT L'AXE Y (1D)
!      I WMAXZI    OUVERTURE MAXIMALE DES FISSURES
!                     EN FACE INFÉRIEURE SUIVANT L'AXE Z (1D)
!      I WMAXZS    OUVERTURE MAXIMALE DES FISSURES
!                     EN FACE SUPÉRIEURE SUIVANT L'AXE Z (1D)
!      I PHIXI     DIAMÈTRE APPROXIMATIF DES ARMATURES INFÉRIEURES SUIVANT X
!      I PHIXS     DIAMÈTRE APPROXIMATIF DES ARMATURES SUPÉRIEURES SUIVANT X
!      I PHIYI     DIAMÈTRE APPROXIMATIF DES ARMATURES INFÉRIEURES SUIVANT Y
!      I PHIYS     DIAMÈTRE APPROXIMATIF DES ARMATURES SUPÉRIEURES SUIVANT Y
!      I PHIZI     DIAMÈTRE APPROXIMATIF DES ARMATURES INFÉRIEURES SUIVANT Z
!      I PHIZS     DIAMÈTRE APPROXIMATIF DES ARMATURES SUPÉRIEURES SUIVANT Z
!      I PRECS     PRECISION SUPPLEMENTAIRE DANS LA RECHERCHE DE L'OPTIMUM
!                   POUR LA METHODE DES 3 PIVOTS (Intervention du 03/2023)
!                     PRECS = 0 (NON)
!                     PRECS = 1 (OUI)
!      I FERRSYME  FERRAILLAGE SYMETRIQUE?
!                     FERRSYME = 0 (NON)
!                     FERRSYME = 1 (OUI)
!      I SLSYME    SECTION SEUIL DE TOLERANCE POUR UN FERRAILLAGE SYMETRIQUE
!      I FERRCOMP  PRISE EN COMPTE DU FERRAILLAGE DE COMPRESSION
!                     FERRCOMP = 0 (NON)
!                     FERRCOMP = 1 (OUI)
!      I UC        UNITE DES CONTRAINTES :
!                     UC = 0 CONTRAINTES EN Pa
!                     UC = 1 CONTRAINTES EN MPa
!      I UM        UNITE DES DIMENSIONS :
!                     UM = 0 DIMENSIONS EN m
!                     UM = 1 DIMENSIONS EN mm
!
!      O DNSITS    DENSITE DES ACIERS CALCULES :
!                     DNSITS(1) = AYI
!                     DNSITS(2) = AYS
!                     DNSITS(3) = AZI
!                     DNSITS(4) = AZS
!                     DNSITS(5) = AST
!                     DNSITS(6) = ATOT = AYI+AYS+AZI+AZS
!      O IERRL     CODE RETOUR LONGI (0 = OK)
!      O IERRT     CODE RETOUR TRNSV (0 = OK)
!
!______________________________________________________________________
!
    implicit none
#include "asterfort/breselsqp.h"
#include "asterfort/cftels.h"
!
    integer(kind=8) :: typco
    real(kind=8) :: cequi
    real(kind=8) :: effrts(6)
    real(kind=8) :: ht
    real(kind=8) :: bw
    real(kind=8) :: enrobyi
    real(kind=8) :: enrobys
    real(kind=8) :: enrobzi
    real(kind=8) :: enrobzs
    real(kind=8) :: facier
    real(kind=8) :: fbeton
    real(kind=8) :: sigelsqp
    real(kind=8) :: kt
    real(kind=8) :: eys
    real(kind=8) :: wmaxyi
    real(kind=8) :: wmaxys
    real(kind=8) :: wmaxzi
    real(kind=8) :: wmaxzs
    real(kind=8) :: phiyi
    real(kind=8) :: phiys
    real(kind=8) :: phizi
    real(kind=8) :: phizs
    integer(kind=8) :: precs
    integer(kind=8) :: ferrsyme
    real(kind=8) :: slsyme
    integer(kind=8) :: ferrcomp
    integer(kind=8) :: epucisa
    integer(kind=8) :: ferrmin
    real(kind=8) :: rholmin
    real(kind=8) :: rhotmin
    integer(kind=8) :: compress
    integer(kind=8) :: uc
    integer(kind=8) :: um
    real(kind=8) :: dnsits(6)
    integer(kind=8) :: ierrl
    integer(kind=8) :: ierrt

!   DEFINITION DES EFFORTS 1D
    real(kind=8) :: effn, effmy, effmz, effty, efftz, effmt
!   DEFINITION DES TERMES DE FERRAILLAGE
    real(kind=8) :: dnsyi, dnsys, dnszi, dnszs, dnstra
!   DEFINITION VARIABLES INTERMEDIAIRES
    real(kind=8) :: sigmsyi, sigmsys, sigmcyi, sigmcys
    real(kind=8) :: sigmszi, sigmszs, sigmczi, sigmczs
    real(kind=8) :: alphay, alphaz
    integer(kind=8) :: pivoty, pivotz, etaty, etatz
    real(kind=8) :: dnstray, thetaby, aky, uky
    real(kind=8) :: dnstraz, thetabz, akz, ukz
    real(kind=8) :: wfinyi, wfinys, wfinzi, wfinzs
    real(kind=8) :: kvarfy, kvarfz
    integer(kind=8) :: ierry, ierrz
    real(kind=8) :: Asl, thetab, ak, uk, unite_pa, unite_m
    real(kind=8) :: Sacier, d, Smoy, fctm
    real(kind=8) :: Calc, ab1, ab2, ab3, ab4
    real(kind=8) :: effrts_fake(8)
    integer(kind=8) :: i, ierr

    do i = 1, 8
        effrts_fake(i) = 0.0
    end do

    !Initialisation des valeurs
    dnsyi = 0
    dnsys = 0
    dnszi = 0
    dnszs = 0
    dnstray = 0
    dnstraz = 0
    alphay = -1
    alphaz = -1
    sigmcyi = 0
    sigmcys = 0
    sigmczi = 0
    sigmczs = 0
    sigmsyi = 0
    sigmsys = 0
    sigmszi = 0
    sigmszs = 0
    pivoty = 0
    pivotz = 0
    etaty = 0
    etatz = 0
    kvarfy = 1
    kvarfz = 1

    effn = effrts(1)
    effmy = effrts(2)
    effmz = effrts(3)
    effty = effrts(4)
    efftz = effrts(5)
    effmt = effrts(6)
    ierrl = 0
    ierrt = 0

    dnsits(1) = -1.d0
    dnsits(2) = -1.d0
    dnsits(3) = -1.d0
    dnsits(4) = -1.d0
    dnsits(5) = -1.d0
    dnsits(6) = -1.d0

    call breselsqp(cequi, effmy, effmz, effn, ht, bw, &
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
!           GESTION DES ALARMES EMISES POUR LES ACIERS DE FLEXION A L'ELS QP
    if (ierr .eq. 1) then
!               Facette en pivot B trop comprimée !
!               Alarme dans te0265 + on sort de la boucle + densité = -1 pour l'élément
        ierrl = 1005
    end if
    if (ierr .eq. 2) then
!               Ferraillage symétrique non possible
!               Alarme dans te0265 + on sort de la boucle + densité = -1 pour l'élément
        ierrl = 10011
    end if
    if (ierr .eq. 3) then
!               Résolution itérative impossible à l'els qp !
!               Alarme dans te0265 + on sort de la boucle + densité = -1 pour l'élément
        ierrl = 1006
    end if
    if (ierr .eq. 4) then
!               Resolution itérative par la methode de Bresler non possible pour FCD
!               Alarme dans te0265 + on sort de la boucle + densité = -1 pour l'élément
        ierrl = 10012
    end if

    !Calcul du ferraillage transversal

    if (ierrl .gt. 0) then
        dnsyi = 0
        dnsys = 0
        dnszi = 0
        dnszs = 0
        alphay = -1
        alphaz = -1
        kvarfy = 1
        kvarfz = 1
        sigmsyi = facier
        sigmsys = facier
        sigmszi = facier
        sigmszs = facier
        sigmcyi = sigelsqp
        sigmcys = sigelsqp
        sigmczi = sigelsqp
        sigmczs = sigelsqp
    end if

    !1e calcul avec MFY et VZ
    Sacier = kvarfz*facier
    call cftels(typco, 1, effrts_fake, effmy, effn, efftz, effmt, &
                dnszi, dnszs, &
                sigmszi, sigmszs, sigmczi, sigmczs, alphaz, &
                ht, bw, enrobzi, enrobzs, facier, fbeton, &
                sigelsqp, sigelsqp, Sacier, uc, um, &
                compress, dnstraz, thetabz, akz, ukz, ierrz)
!               GESTION DES ALARMES EMISES POUR LE FERRAILLAGE TRANSVERSAL A L'ELS
    if (ierrz .eq. 1) then
!                   Béton trop cisaillé !
!                   Alarme dans te0265 + on sort de la boucle + dnstra = -1 pour l'élément
        ierrt = 1007
    end if

    !2e calcul avec MFZ et VY
    Sacier = kvarfy*facier
    call cftels(typco, 1, effrts_fake, effmz, effn, effty, effmt, &
                dnsyi, dnsys, &
                sigmsyi, sigmsys, sigmcyi, sigmcys, alphay, &
                bw, ht, enrobyi, enrobys, facier, fbeton, &
                sigelsqp, sigelsqp, Sacier, uc, um, &
                compress, dnstray, thetaby, aky, uky, ierry)
!               GESTION DES ALARMES EMISES POUR LE FERRAILLAGE TRANSVERSAL A L'ELS
    if (ierry .eq. 1) then
!                   Béton trop cisaillé !
!                   Alarme dans te0265 + on sort de la boucle + dnstra = -1 pour l'élément
        ierrt = 1007
    end if

    !Prise en compte de l'impact de l'effort tranchant sur le ferraillage longitudinal
    !MFY et VZ
    if ((ierrl .eq. 0) .and. (ierrz .eq. 0) .and. &
        & (dnstraz .gt. 0) .and. (epucisa .eq. 1)) then
        Asl = abs(efftz)/(tan(thetabz))
        if (effmy .ge. 0) then
            Calc = abs(sigmszi)
            if (Calc .gt. epsilon(Calc)) then
                Asl = Asl/Calc
            else
                Asl = Asl/(kvarfz*facier)
            end if
            dnszi = dnszi+Asl
        else
            Calc = abs(sigmszs)
            if (Calc .gt. epsilon(Calc)) then
                Asl = Asl/Calc
            else
                Asl = Asl/(kvarfz*facier)
            end if
            dnszs = dnszs+Asl
        end if
    end if

    !MFZ et VY
    if ((ierrl .eq. 0) .and. (ierry .eq. 0) .and. &
        & (dnstray .gt. 0) .and. (epucisa .eq. 1)) then
        Asl = abs(effty)/(tan(thetaby))
        if (effmz .ge. 0) then
            Calc = abs(sigmsyi)
            if (Calc .gt. epsilon(Calc)) then
                Asl = Asl/Calc
            else
                Asl = Asl/(kvarfy*facier)
            end if
            dnsyi = dnsyi+Asl
        else
            Calc = abs(sigmsys)
            if (Calc .gt. epsilon(Calc)) then
                Asl = Asl/Calc
            else
                Asl = Asl/(kvarfy*facier)
            end if
            dnsys = dnsys+Asl
        end if
    end if

    !Choix final du ferraillage transversal
    !Prise en compte de l'impact de la torsion sur le ferraillage longitudinal

    if ((ierry .eq. 0) .and. (ierrz .eq. 0)) then

        if (dnstray .ge. dnstraz) then
            dnstra = dnstray
            thetab = thetaby
            ak = aky
            uk = uky
        else
            dnstra = dnstraz
            thetab = thetabz
            ak = akz
            uk = ukz
        end if

        if ((epucisa .eq. 1) .and. (ierrl .eq. 0)) then
            Asl = (abs(effmt)/(2*ak))*(1/(tan(thetab)))*uk
            if ((abs(effmy) .gt. epsilon(effmy)) .and. (abs(effmz) .gt. epsilon(effmz))) then
                Calc = abs(sigmsyi)
                if (Calc .gt. epsilon(Calc)) then
                    dnsyi = dnsyi+Asl/(4.0*Calc)
                else
                    dnsyi = dnsyi+Asl/(4.0*kvarfy*facier)
                end if
                Calc = abs(sigmsys)
                if (Calc .gt. epsilon(Calc)) then
                    dnsys = dnsys+Asl/(4.0*Calc)
                else
                    dnsys = dnsys+Asl/(4.0*kvarfy*facier)
                end if
                Calc = abs(sigmszi)
                if (Calc .gt. epsilon(Calc)) then
                    dnszi = dnszi+Asl/(4.0*Calc)
                else
                    dnszi = dnszi+Asl/(4.0*kvarfz*facier)
                end if
                Calc = abs(sigmszs)
                if (Calc .gt. epsilon(Calc)) then
                    dnszs = dnszs+Asl/(4.0*Calc)
                else
                    dnszs = dnszs+Asl/(4.0*kvarfz*facier)
                end if
            elseif (abs(effmy) .gt. epsilon(effmy)) then
                Calc = abs(sigmszi)
                if (Calc .gt. epsilon(Calc)) then
                    dnszi = dnszi+Asl/(2.0*Calc)
                else
                    dnszi = dnszi+Asl/(2.0*kvarfz*facier)
                end if
                Calc = abs(sigmszs)
                if (Calc .gt. epsilon(Calc)) then
                    dnszs = dnszs+Asl/(2.0*Calc)
                else
                    dnszs = dnszs+Asl/(2.0*kvarfz*facier)
                end if
            elseif (abs(effmz) .gt. epsilon(effmz)) then
                Calc = abs(sigmsyi)
                if (Calc .gt. epsilon(Calc)) then
                    dnsyi = dnsyi+Asl/(2.0*Calc)
                else
                    dnsyi = dnsyi+Asl/(2.0*kvarfy*facier)
                end if
                Calc = abs(sigmsys)
                if (Calc .gt. epsilon(Calc)) then
                    dnsys = dnsys+Asl/(2.0*Calc)
                else
                    dnsys = dnsys+Asl/(2.0*kvarfy*facier)
                end if
            else
                Smoy = 0.5*(kvarfy+kvarfz)*facier
                dnsyi = dnsyi+Asl/(4.0*Smoy)
                dnsys = dnsys+Asl/(4.0*Smoy)
                dnszi = dnszi+Asl/(4.0*Smoy)
                dnszs = dnszs+Asl/(4.0*Smoy)
            end if
        end if

    else

        dnstra = -1
        thetab = -1

    end if

!  -- VERIFICATION DU FERRAILLAGE MINIMUM :
!  ----------------------------------------

    if ((ferrmin .eq. 1) .or. (ferrmin .eq. 2)) then

        if (uc .eq. 0) then
            unite_pa = 1.e6
            unite_m = 1.
        elseif (uc .eq. 1) then
            unite_pa = 1.
            unite_m = 1.e-3
        end if

        if (fbeton .le. (50*unite_pa)) then
            fctm = 0.30*((fbeton/unite_pa)**(2.0/3.0))
        else
            fctm = 2.12*LOG(1.0+((fbeton/unite_pa)+8.0)/10.0)
        end if

        if (ferrmin .eq. 2) then
            rholmin = max(0.26*(fctm/facier), 0.0013)
            rhotmin = 0.08*(fbeton**(0.5))/facier
        end if

        ierr = ierrl

        !ASYI
        d = bw-enrobyi
        if ((dnsyi .lt. (rholmin*d*ht)) .and. (ierr .eq. 0)) then
            dnsyi = rholmin*d*ht
        end if

        !ASYS
        d = bw-enrobys
        if ((dnsys .lt. (rholmin*d*ht)) .and. (ierr .eq. 0)) then
            dnsys = rholmin*d*ht
        end if

        !ASZI
        d = ht-enrobzi
        if ((dnszi .lt. (rholmin*d*bw)) .and. (ierr .eq. 0)) then
            dnszi = rholmin*d*bw
        end if

        !ASZS
        d = ht-enrobzs
        if ((dnszs .lt. (rholmin*d*bw)) .and. (ierr .eq. 0)) then
            dnszs = rholmin*d*bw
        end if

        ierr = ierrt

        !AST
        if ((dnstra .lt. (rhotmin*max(bw, ht))) .and. (ierr .eq. 0)) then
            dnstra = rhotmin*max(bw, ht)
        end if

    end if

!   RESTITUTION FINALE DES SECTIONS DE FERRAILLAGE
!  -----------------------------------------------
    if (ierrl .eq. 0) then
        dnsits(1) = dnsyi
        dnsits(2) = dnsys
        dnsits(3) = dnszi
        dnsits(4) = dnszs
        dnsits(6) = dnsyi+dnsys+dnszi+dnszs
    end if
    if (ierrt .eq. 0) then
        dnsits(5) = dnstra
    end if

end subroutine
