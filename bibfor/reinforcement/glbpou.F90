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

subroutine glbpou(typcmb, typco, cequi, effrts, ht, bw, &
                  enrobyi, enrobys, enrobzi, enrobzs, &
                  facier, fbeton, sigelsqp, kt, eys, &
                  alphacc, clacier, gammas, gammac, typdiag, &
                  sigcyi, sigcys, sigczi, sigczs, sigs, &
                  wmaxyi, wmaxys, wmaxzi, wmaxzs, &
                  phiyi, phiys, phizi, phizs, &
                  precs, ferrsyme, slsyme, ferrcomp, &
                  epucisa, ferrmin, rholmin, rhotmin, compress, uc, um, &
                  rhoacier, areinf, ashear, astirr, rhocrit, datcrit, lcrit, &
                  dnsits, dnsvol, construc, ierrl, ierrt)
!______________________________________________________________________
!
!      GLBPOU

!      CALCUL GLOBAL DU FERRAILLAGE DES POUTRES
!
!      I TYPCMB    TYPE DE COMBINAISON :
!                       0 = ELU, 1 = ELS, 2 = ELS QP
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
!      I ALPHACC   COEFFICIENT DE SECURITE SUR LA RESISTANCE
!                  DE CALCUL DU BETON EN COMPRESSION
!      I CLACIER   CLASSE DE DUCTILITE DES ACIERS (UTILISE POUR EC2) :
!                     CLACIER = 0 ACIER PEU DUCTILE (CLASSE A)
!                     CLACIER = 1 ACIER MOYENNEMENT DUCTILE (CLASSE B)
!                     CLACIER = 3 ACIER FORTEMENT DUCTILE (CLASSE C)
!      I GAMMAS    COEFFICIENT DE SECURITE SUR LA RESISTANCE
!                  DE CALCUL DES ACIERS
!      I GAMMAC    COEFFICIENT DE SECURITE SUR LA RESISTANCE
!                  DE CALCUL DU BETON
!      I TYPDIAG   TYPE DE DIAGRAMME UTILISÉ POUR L'ACIER
!                     TYPDIAG = 1 ("B1" ==> PALIER INCLINÉ)
!                     TYPDIAG = 2 ("B2" ==> PALIER HORIZONTAL)
!      I SGICYI    CONTRAINTE ULTIME DU BÉTON COMPRIME
!                     EN FIBRE INFERIEURE SUIVANT L'AXE Y À L'ELS
!      I SGICYS    CONTRAINTE ULTIME DU BÉTON COMPRIME
!                     EN FIBRE SUPERIEURE SUIVANT L'AXE Y À L'ELS
!      I SGICZI    CONTRAINTE ULTIME DU BÉTON COMPRIME
!                     EN FIBRE INFERIEURE SUIVANT L'AXE Z À L'ELS
!      I SGICZS    CONTRAINTE ULTIME DU BÉTON COMPRIME
!                     EN FIBRE SUPERIEURE SUIVANT L'AXE Z À L'ELS
!      I SGIS      CONTRAINTE ULTIME DE L'ACIER À L'ELS
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
!      I PRECS     PRECISION ITERATION
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
!      I RHOACIER  MASSE VOLUMIQUE DES ACIERS
!      I AREINF    COEFF DE PONDER DU RATIO DE DENSITÉ D'ACIER PAR MÈTRE CUBE DE BÉTON
!      I ASHEAR    COEFF DE PONDER DU RATIO DE DENSITÉ D'ACIER D'EFFORT TRANCHANT
!      I ASTIRR    COEFF DE PONDER DU RATIO DE LONGUEUR DES ÉPINGLES D'ACIER EFF TRANC
!      I RHOCRIT   DENSITÉ VOLUMIQUE D'ARMATURE CRITIQUE
!      I DATCRIT   FERRAILLAGE D'EFFORT TRANCHANT CRITIQUE
!      I LCRIT     LONGUEUR CRITIQUE DES EPINGLE D'ACIERS D'EFFORT TRANCHANT
!
!      O DNSITS    DENSITE DES ACIERS CALCULES :
!                     DNSITS(1) = AYI
!                     DNSITS(2) = AYS
!                     DNSITS(3) = AZI
!                     DNSITS(4) = AZS
!                     DNSITS(5) = AST
!                     DNSITS(6) = ATOT = AYI+AYS+AZI+AZS
!      O DNSVOL    DENSITE VOLUMIQUE D'ARMATURE (Kg/M3)
!      O CONSTRUC  INDICATEUR DE COMPLEXITE DE CONSTRUCTIBILITE (-)
!      O IERRL         CODE RETOUR LONGI (0 = OK)
!      O IERRT         CODE RETOUR TRNSV (0 = OK)
!
!______________________________________________________________________
!
    implicit none
#include "asterfort/glbelu.h"
#include "asterfort/glbels.h"
#include "asterfort/glbelsqp.h"
!
    integer(kind=8) :: typcmb
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
    real(kind=8) :: alphacc
    integer(kind=8) :: clacier
    real(kind=8) :: gammas
    real(kind=8) :: gammac
    integer(kind=8) :: typdiag
    real(kind=8) :: sigcyi
    real(kind=8) :: sigcys
    real(kind=8) :: sigczi
    real(kind=8) :: sigczs
    real(kind=8) :: sigs
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
    real(kind=8) :: rhoacier
    real(kind=8) :: areinf
    real(kind=8) :: ashear
    real(kind=8) :: astirr
    real(kind=8) :: rhocrit
    real(kind=8) :: datcrit
    real(kind=8) :: lcrit
    real(kind=8) :: dnsits(6)
    real(kind=8) :: dnsvol
    real(kind=8) :: construc
    integer(kind=8) :: ierrl
    integer(kind=8) :: ierrt

!! AUTRES VARIABLES
    real(kind=8) :: shear, reinf, stirrups, Calc

    if (typcmb .eq. 0) then

        call glbelu(typco, alphacc, effrts, ht, bw, &
                    enrobyi, enrobys, enrobzi, enrobzs, &
                    facier, fbeton, gammas, gammac, &
                    clacier, eys, typdiag, precs, &
                    ferrsyme, slsyme, ferrcomp, &
                    epucisa, ferrmin, rholmin, rhotmin, compress, uc, um, &
                    dnsits, ierrl, ierrt)

    elseif (typcmb .eq. 1) then

        call glbels(typco, cequi, effrts, ht, bw, &
                    enrobyi, enrobys, enrobzi, enrobzs, &
                    facier, fbeton, sigcyi, sigcys, sigczi, sigczs, sigs, &
                    precs, ferrsyme, slsyme, ferrcomp, &
                    epucisa, ferrmin, rholmin, rhotmin, compress, uc, um, &
                    dnsits, ierrl, ierrt)

    elseif (typcmb .eq. 2) then

        call glbelsqp(typco, cequi, effrts, ht, bw, &
                      enrobyi, enrobys, enrobzi, enrobzs, &
                      facier, fbeton, sigelsqp, kt, eys, &
                      wmaxyi, wmaxys, wmaxzi, wmaxzs, &
                      phiyi, phiys, phizi, phizs, &
                      precs, ferrsyme, slsyme, ferrcomp, &
                      epucisa, ferrmin, rholmin, rhotmin, compress, uc, um, &
                      dnsits, ierrl, ierrt)
    end if

!   -- CALCUL DE LA DENSITE VOLUMIQUE D'ARMATURE :
!   ----------------------------------------------
!
    if ((rhoacier .gt. 0) .and. (ierrl .eq. 0) .and. (ierrt .eq. 0)) then
        dnsvol = rhoacier*(dnsits(1)+dnsits(2)+dnsits(3)+dnsits(4)+dnsits(5)*max(ht, bw))/(ht*bw)
    else
        dnsvol = -1.d0
    end if
!
!   -- CALCUL DE L'INDICATEUR DE CONSTRUCTIBILITE :
!   -----------------------------------------------
!
    if (rhoacier .gt. 0) then
        reinf = areinf*DNSVOL/rhocrit
        shear = ashear*dnsits(5)/datcrit
        if (shear .lt. 0.d0) shear = 0.d0
        stirrups = astirr*dnsits(5)*(ht-enrobys-enrobyi)*(bw-enrobzs-enrobzi)/(datcrit*lcrit)
        if (stirrups .lt. 0.d0) stirrups = 0.d0
        construc = (reinf+shear+stirrups)/(areinf+ashear+astirr)
    else
        construc = -1.d0
    end if

end subroutine
