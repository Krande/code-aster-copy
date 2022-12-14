! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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

subroutine te0146(option, nomte)
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/clcplq.h"
#include "asterfort/jevech.h"
#include "asterfort/tecach.h"
#include "asterfort/tecael.h"
#include "asterfort/utmess.h"
!
    character(len=16) :: option, nomte
!
!.....................................................................
!  BUT: CALCUL DE L'OPTION FERRAILLAGE POUR LES ELEMENTS DE COQUE
!.....................................................................
!_____________________________________________________________________
!
! CALCUL DES DENSITES DE FERRAILLAGE DANS LE BETON ARME
!              (METHODE DE CAPRA ET MAURY)
!
! VERSION DU 24/09/2021
!_____________________________________________________________________
!
! PARAMETRES D'ECHANGE ENTRE CODE_ASTER ET CLCPLQ 
! (POINT D'ENTREE DU CALCUL DE FERRAILLAGE PAR CAPRA ET MAURY)
!
!   PARAMETRES D'ENTREE (FOURNIS PAR CODE_ASTER)
!
!     TYPCMB     TYPE DE COMBINAISON :
!                   0 = ELU, 1 = ELS, 2 = ELS QP
!     TYPCO      TYPE DE CODIFICATION :
!                   1 = BAEL91, 2 = EUROCODE 2
!     TYPSTRU    TYPE DE STRUCTURE :
!                   0 = 2D, 1 = 1D
!     FERRSYME   FERRAILLAGE SYMETRIQUE?
!                   0 = NON, 1 = OUI
!     SLSYME     SECTION SEUIL DE TOLERANCE POUR UN FERRAILLAGE SYMETRIQUE 
!     FERRCOMP   FERRAILLAGE DE COMPRESSION ADMIS?
!                   0 = NON, 1 = OUI
!     EPUCISA    IMPACT DE L'EFFORT TRANCHANT ET DE LA TORSION SUR LE
!                   FERRAILLAGE LONGITUDINAL?
!                   0 = NON, 1 = OUI
!     FERRMIN    PRISE EN COMPTE D'UN FERRAILLAGE MINIMUN
!                   0 = NON, 1 = OUI, 2 = CODE
!     RHOLMIN    RATIO DE FERRAILLAGE LONGI MINI (A RENSEIGNER SI FERMIN='OUI')
!     RHOTMIN    RATIO DE FERRAILLAGE TRNSV MINI (A RENSEIGNER SI FERMIN='OUI')
!     COMPRESS   VALORISATION DE LA COMPRESSION POUR LES ACIERS TRANSVERSAUX
!                   0 = COMPRESSION NON PRISE EN COMPTE
!                   1 = COMPRESSION PRISE EN COMPTE
!     CEQUI      COEFFICIENT D'EQUIVALENCE ACIER/BETON
!     ENROBI     ENROBAGE DES ARMATURES INFERIEURES (2D)
!     ENROBS     ENROBAGE DES ARMATURES SUPERIEURES (2D)
!     ENROBYI    ENROBAGE DES ARMATURES INFERIEURES SUIVANT L'AXE Y (1D)
!     ENROBYS    ENROBAGE DES ARMATURES SUPERIEURES SUIVANT L'AXE Y (1D)
!     ENROBZI    ENROBAGE DES ARMATURES INFERIEURES SUIVANT L'AXE Z (1D)
!     ENROBZS    ENROBAGE DES ARMATURES SUPERIEURES SUIVANT L'AXE Z (1D)
!     SIGS       CONTRAINTE ULTIME DES ACIERS ?? L'ELS
!     SGICI      CONTRAINTE ULTIME DU B??TON COMPRIME EN FIBRE INFERIEURE ?? L'ELS (2D)
!     SGICS      CONTRAINTE ULTIME DU B??TON COMPRIME EN FIBRE SUPERIEURE ?? L'ELS (2D)
!     SGICYI     CONTRAINTE ULTIME DU B??TON COMPRIME
!                   EN FIBRE INFERIEURE SUIVANT L'AXE Y ?? L'ELS (1D)
!     SGICYS     CONTRAINTE ULTIME DU B??TON COMPRIME
!                   EN FIBRE SUPERIEURE SUIVANT L'AXE Y ?? L'ELS (1D)
!     SGICZI     CONTRAINTE ULTIME DU B??TON COMPRIME
!                   EN FIBRE INFERIEURE SUIVANT L'AXE Z ?? L'ELS (1D)
!     SGICZS     CONTRAINTE ULTIME DU B??TON COMPRIME
!                   EN FIBRE SUPERIEURE SUIVANT L'AXE Z ?? L'ELS (1D)
!     ALPHACC    COEFFICIENT DE SECURITE SUR LA RESISTANCE
!                   DE CALCUL DU BETON EN COMPRESSION
!     GAMMAS     COEFFICIENT DE SECURITE SUR LA RESISTANCE
!                   DE CALCUL DES ACIERS
!     GAMMAC     COEFFICIENT DE SECURITE SUR LA RESISTANCE
!                   DE CALCUL DU BETON
!     FACIER     LIMITE D'ELASTICITE DES ACIERS (CONTRAINTE)
!     EYS        MODULE D'YOUNG DE L'ACIER
!     TYPDIAG    TYPE DE DIAGRAMME UTILIS?? POUR L'ACIER
!                   TYPDIAG = 1 ("B1" ==> PALIER INCLIN??)
!                   TYPDIAG = 2 ("B2" ==> PALIER HORIZONTAL)
!     FBETON     RESISTANCE EN COMPRESSION DU BETON (CONTRAINTE)
!     CLACIER    CLASSE DE DUCTILITE DES ACIERS (POUR L'EC2) :
!                   CLACIER = 0 ACIER PEU DUCTILE (CLASSE A)
!                   CLACIER = 1 ACIER MOYENNEMENT DUCTILE (CLASSE B)
!                   CLACIER = 2 ACIER FORTEMENT DUCTILE (CLASSE C)
!     UC         UNITE DES CONTRAINTES :
!                   0 = CONTRAINTES EN Pa
!                   1 = CONTRAINTES EN MPa
!     UM         UNITE DES DIMENSIONS :
!                   0 = DIMENSIONS EN m
!                   1 = DIMENSIONS EN mm
!     RHOACIER   MASSE VOLUMIQUE DES ACIERS
!     AREINF     COEFF DE PONDER DU RATIO DE DENSIT?? D'ACIER PAR M??TRE CUBE DE B??TON
!     ASHEAR     COEFF DE PONDER DU RATIO DE DENSIT?? D'ACIER D'EFFORT TRANCHANT
!     ASTIRR     COEFF DE PONDER DU RATIO DE LONGUEUR DES ??PINGLES D'ACIER EFF TRANC
!     RHOCRIT    DENSIT?? VOLUMIQUE D'ARMATURE CRITIQUE
!     DATCRIT    FERRAILLAGE D'EFFORT TRANCHANT CRITIQUE
!     LCRIT      LONGUEUR CRITIQUE DES EPINGLE D'ACIERS D'EFFORT TRANCHANT
!     WMAXI      OUVERTURE MAXIMALE DES FISSURES EN FACE INF??RIEURE (2D)
!     WMAXS      OUVERTURE MAXIMALE DES FISSURES EN FACE SUP??RIEURE (2D)
!     WMAXYI     OUVERTURE MAXIMALE DES FISSURES
!                   EN FACE INF??RIEURE SUIVANT L'AXE Y (1D)
!     WMAXYS     OUVERTURE MAXIMALE DES FISSURES
!                   EN FACE SUP??RIEURE SUIVANT L'AXE Y (1D)
!     WMAXZI     OUVERTURE MAXIMALE DES FISSURES
!                   EN FACE INF??RIEURE SUIVANT L'AXE Z (1D)
!     WMAXZS     OUVERTURE MAXIMALE DES FISSURES
!                   EN FACE SUP??RIEURE SUIVANT L'AXE Z (1D)
!     SIGELSQP   CONTRAINTE ADMISSIBLE DANS LE BETON ?? L'ELS QP
!     KT         COEFFICIENT DE DUR??E DE CHARGEMENT
!     PHIXI      DIAM??TRE APPROXIMATIF DES ARMATURES INF??RIEURES SUIVANT X
!     PHIXS      DIAM??TRE APPROXIMATIF DES ARMATURES SUP??RIEURES SUIVANT X
!     PHIYI      DIAM??TRE APPROXIMATIF DES ARMATURES INF??RIEURES SUIVANT Y
!     PHIYS      DIAM??TRE APPROXIMATIF DES ARMATURES SUP??RIEURES SUIVANT Y
!     PHIZI      DIAM??TRE APPROXIMATIF DES ARMATURES INF??RIEURES SUIVANT Z
!     PHIZS      DIAM??TRE APPROXIMATIF DES ARMATURES SUP??RIEURES SUIVANT Z

!     HT         EPAISSEUR DE LA COQUE
!     EFFRTS     TORSEUR DES EFFORTS ET DES MOMENTS (DIM 8)
!
!   PARAMETRES DE SORTIE (RENVOYES A CODE_ASTER)
!
!     DNSITS     DENSITES DE FERRAILLAGE (DIM 5) :
!                   1 A 4 : ACIER LONGITUDINAL (EN M2/M)
!                   5 A 6 : ACIERS TRANSVERSAUX (EN M2/M2)
!     DNSVOL     DENSITE VOLUMIQUE D'ARMATURE (Kg/M3)
!     CONSTRUC   INDICATEUR DE COMPLEXITE DE CONSTRUCTIBILITE (-)
!     IERR       CODE RETOUR (0 = OK)
!---------------------------------------------------------------------
!
    real(kind=8) :: cequi, sigs, sigci, sigcs, sigcyi, sigcys, sigczi, sigczs
    real(kind=8) :: alphacc, effrts(8), dnsits(6)
    real(kind=8) :: ht, enrobi, enrobs, enrobyi, enrobys, enrobzi, enrobzs
    real(kind=8) :: gammac, gammas, rholmin, rhotmin, slsyme
    real(kind=8) :: facier, fbeton, eys, rhoacier, dnsvol, areinf, ashear
    real(kind=8) :: astirr, rhocrit, datcrit, lcrit, construc
    real(kind=8) :: wmaxi, wmaxs, wmaxyi, wmaxys, wmaxzi, wmaxzs, sigelsqp, kt
    real(kind=8) :: phixi, phixs, phiyi, phiys, phizi, phizs
    real(kind=8) :: reinf, shear, stirrups
    integer :: ierr, jepais, jefge, jfer1, jfer2, itab(7), nno
    integer :: typcmb, typco, ferrmin, typdiag, ferrsyme, epucisa, clacier, uc, um
    integer :: ino, icmp, iret, k
    integer :: iadzi, iazk24, compress, ferrcomp, typstru
!
    call tecael(iadzi, iazk24, noms=0)
!
    call jevech('PCACOQU', 'L', jepais)
    call jevech('PFERRA1', 'L', jfer1)
    call jevech('PFERRA2', 'E', jfer2)
    ht = zr(jepais)
!
    call jevech('PEFFORR', 'L', jefge)
    call tecach('OOO', 'PEFFORR', 'L', iret, nval=7, itab=itab)
    ASSERT(iret.eq.0)
    nno = itab(3)
    ASSERT(nno.gt.0.and.nno.le.9)
    ASSERT(itab(2).eq.8*nno)
    
!
!       -- CALCUL DE LA CONTRAINTE MOYENNE :
!       ------------------------------------
!
    do icmp = 1, 8
        effrts(icmp) = 0.d0
        do ino = 1, nno
            effrts(icmp) = effrts(icmp) + zr(jefge-1+(ino-1)*8+icmp)/nno
        end do
    end do
!
!       -- RECUPERATION DES DONNEES DE L'UTILISATEUR :
!       ----------------------------------------------
!     FER1_R = 'TYPCOMB','CODIF','TYPSTRU','FERRSYME','SLSYME','FERRCOMP',
!                  1        2        3          4         5        6
!              'EPUCISA','FERRMIN','RHOLMIN','RHOTMIN',
!                  7         8         9        10
!              'COMPRESS','CEQUI','ENROBI','ENROBS','ENROBYI','ENROBYS',
!                  11        12      13       14       15        16
!              'ENROBZI','ENROBZS','SIGS','SIGCI','SIGCS','SIGCYI','SIGCYS',
!                  17        18      19      20     21       22       23
!              'SIGCZI','SIGCZS','ALPHACC','GAMMAS','GAMMAC','FACIER','EYS',
!                  24      25        26       27       28       29     30
!              'TYPDIAG','FBETON','CLACIER','UC','UM','RHOACIER','AREINF',
!                  31       32       33      34   35      36        37
!              'ASHEAR','ASTIRR','RHOCRIT','DATCRIT','LCRIT','WMAXI','WMAXS',
!                  38      39        40        41      42      43      44
!              'WMAXYI','WMAXYS','WMAXZI','WMAXZS','SIGELSQP','KT',
!                 45       46       47       48        49      50
!              'PHIXI','PHIXS','PHIYI','PHIYS','PHIZI','PHIZS'
!                 51      52      53      54      55      56
!
    typcmb = nint(zr(jfer1-1+1))
    typco = nint(zr(jfer1-1+2))
    typstru = nint(zr(jfer1-1+3))
    ferrsyme = nint(zr(jfer1-1+4))
    slsyme = zr(jfer1-1+5)
    ferrcomp = nint(zr(jfer1-1+6))
    epucisa = nint(zr(jfer1-1+7))
    ferrmin = nint(zr(jfer1-1+8))
    rholmin = zr(jfer1-1+9)
    rhotmin = zr(jfer1-1+10)
    compress = int(zr(jfer1-1+11))
    cequi = zr(jfer1-1+12)
    enrobi = zr(jfer1-1+13)
    enrobs = zr(jfer1-1+14)
    enrobyi = zr(jfer1-1+15)
    enrobys = zr(jfer1-1+16)
    enrobzi = zr(jfer1-1+17)
    enrobzs = zr(jfer1-1+18)
    sigs = zr(jfer1-1+19)
    sigci = zr(jfer1-1+20)
    sigcs = zr(jfer1-1+21)
    sigcyi = zr(jfer1-1+22)
    sigcys = zr(jfer1-1+23)
    sigczi = zr(jfer1-1+24)
    sigczs = zr(jfer1-1+25)
    alphacc = zr(jfer1-1+26)
    gammas = zr(jfer1-1+27)
    gammac = zr(jfer1-1+28)
    facier = zr(jfer1-1+29)
    eys = zr(jfer1-1+30)
    typdiag = int(zr(jfer1-1+31))
    fbeton = zr(jfer1-1+32)
    clacier = int(zr(jfer1-1+33))
    uc = int(zr(jfer1-1+34))
    um = int(zr(jfer1-1+35))
    rhoacier = zr(jfer1-1+36)
    areinf = zr(jfer1-1+37)
    ashear = zr(jfer1-1+38)
    astirr = zr(jfer1-1+39)
    rhocrit = zr(jfer1-1+40)
    datcrit = zr(jfer1-1+41)
    lcrit = zr(jfer1-1+42)
    wmaxi = zr(jfer1-1+43)
    wmaxs = zr(jfer1-1+44)
    wmaxyi = zr(jfer1-1+45)
    wmaxys = zr(jfer1-1+46)
    wmaxzi = zr(jfer1-1+47)
    wmaxzs = zr(jfer1-1+48)
    sigelsqp = zr(jfer1-1+49)
    kt = zr(jfer1-1+50)
    phixi = zr(jfer1-1+51)
    phixs = zr(jfer1-1+52)
    phiyi = zr(jfer1-1+53)
    phiys = zr(jfer1-1+54)
    phizi = zr(jfer1-1+55)
    phizs = zr(jfer1-1+56)
    
   !Only option '2D'
    if (typstru.eq.1.d0) then
    call utmess('A', 'CALCULEL_78')
    goto 998
    endif

!
!   -- CALCUL PROPREMENT DIT :
!   --------------------------
!
!   VERIFICATION DE LA COHERENCE DES PARAMETRES
    if (ht.le.enrobi .or. ht.le.enrobs) then
        call utmess('F', 'CALCULEL_72')
    endif
!
!       -- INITIALISATION DES VARIABLES DE SORTIES :
!       --------------------------------------------
!
    dnsvol = 0.d0
    construc = 0.d0
    do k = 1, 6
        dnsits(k) = 0.d0
    end do
!
    call clcplq(typcmb, typco, ferrsyme, slsyme, ferrcomp, epucisa,&
                  ferrmin, rholmin, rhotmin, compress, cequi,&
                  enrobi, enrobs, sigs, sigci, sigcs,&
                  alphacc, gammas, gammac, facier, eys, typdiag,&
                  fbeton, clacier, uc, um,&
                  wmaxi, wmaxs, sigelsqp, kt, phixi, phixs, phiyi, phiys,&
                  ht, effrts, dnsits, ierr)

!
!       -- CALCUL DE LA DENSITE VOLUMIQUE D'ARMATURE :
!       ----------------------------------------------
!
    if (rhoacier.gt.0) then
        dnsvol = rhoacier*(dnsits(1)+dnsits(2)+dnsits(3)+dnsits(4)+dnsits(5)*ht+dnsits(6)*ht)/ht
        if (dnsits(5).eq.-1.d0 .or. dnsits(6).eq.-1.d0) then
!           Vrai uniquement pour le calcul du ferraillage transversal au BAEL
!           (pour lequel les aciers d'effort tranchant ne sont pas calcul??s)
            dnsvol = rhoacier*((dnsits(1)+dnsits(2)+dnsits(3)+dnsits(4)))/ht
        endif
    else
        dnsvol = -1.d0
    endif
!
!       -- CALCUL DE L'INDICATEUR DE CONSTRUCTIBILITE :
!       -----------------------------------------------
!
    if (rhoacier.gt.0) then
        reinf = areinf*dnsvol/rhocrit
        shear = ashear*dnsits(5)/datcrit
        if (shear.lt.0.d0) shear = 0.d0
        stirrups = astirr*((dnsits(5)*dnsits(5)+dnsits(6)*dnsits(6))**0.5)
        stirrups = stirrups*(ht-enrobs-enrobi)/(datcrit*lcrit)
        if (stirrups.lt.0.d0) stirrups = 0.d0
        construc = (reinf+shear+stirrups)/(areinf+ashear+astirr)
    else
        construc = -1.d0
    endif
!

!       -- GESTION DES ALARMES EMISES :
!       -------------------------------
!
    if (ierr.eq.1001) then
!       ELU : section trop comprim??e
!       on fixe toutes les densit??s de ferraillage de l'??l??ment ?? -1
        call utmess('A', 'CALCULEL_83')
        do k = 1, 6
            dnsits(k) = -1.d0
        end do
        dnsvol = -1.d0
        construc = -1.d0
    endif
    
    if (ierr.eq.10011) then
!       ELU : ferraillage sym??trique non possible!
!       on fixe toutes les densit??s de ferraillage de l'??l??ment ?? -1
        call utmess('A', 'CALCULEL7_28')
        do k = 1, 6
            dnsits(k) = -1.d0
        end do
        dnsvol = -1.d0
        construc = -1.d0
    endif
!
    if (ierr.eq.1003) then
!       ELS : section trop comprim??e
!       on fixe toutes les densit??s de ferraillage de l'??l??ment ?? -1
        call utmess('A', 'CALCULEL_84')
        do k = 1, 6
            dnsits(k) = -1.d0
        end do
        dnsvol = -1.d0
        construc = -1.d0
    endif
    
    if (ierr.eq.1005) then
!       ELS_QP : section trop comprim??e
!       on fixe toutes les densit??s de ferraillage de l'??l??ment ?? -1
        call utmess('A', 'CALCULEL_85')
        do k = 1, 6
            dnsits(k) = -1.d0
        end do
        dnsvol = -1.d0
        construc = -1.d0
    endif
!
    if (ierr.eq.1002) then
!       ELU BETON TROP CISAILLE : densit?? transversale fix??e ?? -1 pour l'??l??ment
        call utmess('A', 'CALCULEL_81')
        dnsits(5) = -1.d0
        dnsits(6) = -1.d0
        dnsvol = -1.d0
        construc = -1.d0
    endif

    if (ierr.eq.1004) then
!       ELS BETON TROP CISAILLE : densit?? transversale fix??e ?? -1 pour l'??l??ment
        call utmess('A', 'CALCULEL_86')
        dnsits(5) = -1.d0
        dnsits(6) = -1.d0
        dnsvol = -1.d0
        construc = -1.d0
    endif

    if (ierr.eq.1007) then
!       ELS_QP BETON TROP CISAILLE : densit?? transversale fix??e ?? -1 pour l'??l??ment
        call utmess('A', 'CALCULEL_87')
        dnsits(5) = -1.d0
        dnsits(6) = -1.d0
        dnsvol = -1.d0
        construc = -1.d0
    endif

    if (ierr.eq.1006) then
!       ELS QP SOLLICITATION TROP IMPORTANTE : R??solution it??rative impossible ?? l'els qp !
        call utmess('A', 'CALCULEL_77')
        do k = 1, 6
            dnsits(k) = -1.d0
        end do
        dnsvol = -1.d0
        construc = -1.d0
    endif
!
!       -- STOCKAGE DES RESULTATS DANS FER2 :
!       -------------------------------------
!   FER2_R =  DNSXI DNSXS DNSYI DNSYS DNSXT DNSYT DNSVOL CONSTRUC
!               1     2     3     4     5     6      7       8
!
    zr(jfer2-1+1)= dnsits(1)
    zr(jfer2-1+2)= dnsits(3)
    zr(jfer2-1+3)= dnsits(2)
    zr(jfer2-1+4)= dnsits(4)
    zr(jfer2-1+5)= dnsits(5)
    zr(jfer2-1+6)= dnsits(6)
    zr(jfer2-1+7)= dnsvol
    zr(jfer2-1+8)= construc

998 continue

end subroutine
