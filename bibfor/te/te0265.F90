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

subroutine te0265(nomopt, nomte)
    implicit none
#include "jeveux.h"
#include "asterfort/jevech.h"
#include "asterfort/tecach.h"
#include "asterfort/tecael.h"
#include "asterfort/utmess.h"
#include "asterfort/glbpou.h"
!.....................................................................
!BUT: CALCUL DE L'OPTION FERRAILLAGE POUR LES ELEMENTS POUTRES/POTEAUX
!.....................................................................
!_____________________________________________________________________
!
! CALCUL DES DENSITES DE FERRAILLAGE DANS LE BETON ARME
!
! VERSION DU 24/09/2021
!_____________________________________________________________________

! PARAMETRES D'ECHANGE ENTRE CODE_ASTER ET CAF(ELU/ELS/ES_QP)
!
!   PARAMETRES D'ENTREE (FOURNIS PAR CODE_ASTER)
!
!     TYPCMB     TYPE DE COMBINAISON :
!                   0 = ELU, 1 = ELS, 2 = ELS QP
!     TYPCO      TYPE DE CODIFICATION :
!                   1 = BAEL91, 2 = EUROCODE 2
!     METH2D     CHOIX DE METHODE POUR LE CALCUL DU FERRAILLAGE DES
!                   ELEMENTS 2D
!                   1 = CapraMaury, 2 = Sandwich/Multicouches
!     THITER     ANGLE D'ITERATION DANS LE CAS DU CALCUL DU FERRAILLAGE DES PLAQUES
!                   (POUR LA METHODE SANDWICH, IL S'AGIT DE L'INCLINAISON DES BIELLES)
!                   (POUR LA METHODE CAPRA-MAURY, IL S'AGIT DE L'ORIENTATION DES FACETTES)
!     COND109    PRISE EN COMPTE OU PAS DE LA CLAUSE §109 POUR LE CALCUL
!                   DE RESISTANCE DANS LA METHODE SANDWICH
!                   1 = OUI, 2 = NON
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
!     SIGS       CONTRAINTE ULTIME DES ACIERS À L'ELS
!     SGICI      CONTRAINTE ULTIME DU BÉTON COMPRIME EN FIBRE INFERIEURE À L'ELS (2D)
!     SGICS      CONTRAINTE ULTIME DU BÉTON COMPRIME EN FIBRE SUPERIEURE À L'ELS (2D)
!     SGICYI     CONTRAINTE ULTIME DU BÉTON COMPRIME
!                   EN FIBRE INFERIEURE SUIVANT L'AXE Y À L'ELS (1D)
!     SGICYS     CONTRAINTE ULTIME DU BÉTON COMPRIME
!                   EN FIBRE SUPERIEURE SUIVANT L'AXE Y À L'ELS (1D)
!     SGICZI     CONTRAINTE ULTIME DU BÉTON COMPRIME
!                   EN FIBRE INFERIEURE SUIVANT L'AXE Z À L'ELS (1D)
!     SGICZS     CONTRAINTE ULTIME DU BÉTON COMPRIME
!                   EN FIBRE SUPERIEURE SUIVANT L'AXE Z À L'ELS (1D)
!     ALPHACC    COEFFICIENT DE SECURITE SUR LA RESISTANCE
!                   DE CALCUL DU BETON EN COMPRESSION
!     GAMMAS     COEFFICIENT DE SECURITE SUR LA RESISTANCE
!                   DE CALCUL DES ACIERS
!     GAMMAC     COEFFICIENT DE SECURITE SUR LA RESISTANCE
!                   DE CALCUL DU BETON
!     FACIER     LIMITE D'ELASTICITE DES ACIERS (CONTRAINTE)
!     EYS        MODULE D'YOUNG DE L'ACIER
!     TYPDIAG    TYPE DE DIAGRAMME UTILISÉ POUR L'ACIER
!                   TYPDIAG = 1 ("B1" ==> PALIER INCLINÉ)
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
!     AREINF     COEFF DE PONDER DU RATIO DE DENSITÉ D'ACIER PAR MÈTRE CUBE DE BÉTON
!     ASHEAR     COEFF DE PONDER DU RATIO DE DENSITÉ D'ACIER D'EFFORT TRANCHANT
!     ASTIRR     COEFF DE PONDER DU RATIO DE LONGUEUR DES ÉPINGLES D'ACIER EFF TRANC
!     RHOCRIT    DENSITÉ VOLUMIQUE D'ARMATURE CRITIQUE
!     DATCRIT    FERRAILLAGE D'EFFORT TRANCHANT CRITIQUE
!     LCRIT      LONGUEUR CRITIQUE DES EPINGLE D'ACIERS D'EFFORT TRANCHANT
!     WMAXI      OUVERTURE MAXIMALE DES FISSURES EN FACE INFÉRIEURE (2D)
!     WMAXS      OUVERTURE MAXIMALE DES FISSURES EN FACE SUPÉRIEURE (2D)
!     WMAXYI     OUVERTURE MAXIMALE DES FISSURES
!                   EN FACE INFÉRIEURE SUIVANT L'AXE Y (1D)
!     WMAXYS     OUVERTURE MAXIMALE DES FISSURES
!                   EN FACE SUPÉRIEURE SUIVANT L'AXE Y (1D)
!     WMAXZI     OUVERTURE MAXIMALE DES FISSURES
!                   EN FACE INFÉRIEURE SUIVANT L'AXE Z (1D)
!     WMAXZS     OUVERTURE MAXIMALE DES FISSURES
!                   EN FACE SUPÉRIEURE SUIVANT L'AXE Z (1D)
!     SIGELSQP   CONTRAINTE ADMISSIBLE DANS LE BETON À L'ELS QP
!     KT         COEFFICIENT DE DURÉE DE CHARGEMENT
!     PHIXI      DIAMÈTRE APPROXIMATIF DES ARMATURES INFÉRIEURES SUIVANT X
!     PHIXS      DIAMÈTRE APPROXIMATIF DES ARMATURES SUPÉRIEURES SUIVANT X
!     PHIYI      DIAMÈTRE APPROXIMATIF DES ARMATURES INFÉRIEURES SUIVANT Y
!     PHIYS      DIAMÈTRE APPROXIMATIF DES ARMATURES SUPÉRIEURES SUIVANT Y
!     PHIZI      DIAMÈTRE APPROXIMATIF DES ARMATURES INFÉRIEURES SUIVANT Z
!     PHIZS      DIAMÈTRE APPROXIMATIF DES ARMATURES SUPÉRIEURES SUIVANT Z
!     HT         HAUTEUR DE LA POUTRE
!     BW         LARGEUR DE LA POUTRE
!     EFFRTS     TORSEUR DES EFFORTS ET DES MOMENTS (DIM 6)
!                     EFFRTS(1) = EFFN / EFFORT NORMAL
!                     EFFRTS(2) = EFFMY / MOMENT FLÉCHISSANT SUIVANT Y
!                     EFFRTS(3) = EFFMZ / MOMENT FLÉCHISSANT SUIVANT Z
!                     EFFRTS(4) = EFFTY / EFFORT TRANCHANT SUIVANT Y
!                     EFFRTS(5) = EFFTZ / EFFORT TRANCHANT SUIVANT Z
!                     EFFRTS(6) = EFFMT / MOMENT DE TORSION
!
!   PARAMETRES DE SORTIE (RENVOYES A CODE_ASTER)
!
!     DNSITS     DENSITES DE FERRAILLAGE (DIM 6) :
!                     DNSITS(1) = AYI
!                     DNSITS(2) = AYS
!                     DNSITS(3) = AZI
!                     DNSITS(4) = AZS
!                     DNSITS(5) = AST
!                     DNSITS(6) = ATOT = AYI+AYS+AZI+AZS
!     DNSVOL     DENSITE VOLUMIQUE D'ARMATURE (Kg/M3)
!     CONSTRUC   INDICATEUR DE COMPLEXITE DE CONSTRUCTIBILITE (-)
!     IERR       CODE RETOUR (0 = OK)
!---------------------------------------------------------------------

    character(len=16) :: nomopt, nomte
    ! Return value for differents ops
    integer(kind=8) :: iret
    ! Max is 7
    integer(kind=8) :: itabin(3)
    integer(kind=8) :: EFGE_ELNO_pointer
    integer(kind=8) :: EFGE_ELNO_length
    integer(kind=8) :: EFGE_ELNO_nbOfPoints
    integer(kind=8) :: EFGE_ELNO_nbOfComponents
    integer(kind=8) :: jcagepo, jfer1, jfer2, jefge

!   INFO ON MAILLE AND ELEMENT
    ! For tecael
    integer(kind=8) :: iadzi, iazk24
    character(len=24) :: meshName, elementName
!   For tecael
    integer(kind=8) :: node_nb, node_i_id, node_j_id

!   OUTPUT FOR 1D ELEMENTS
    real(kind=8) :: dnsvol, construc

!   GEOMETRICAL DATA
    real(kind=8) :: HY1, HZ1, TSEC
!   EFGE_ELNO DATA
    real(kind=8) :: Ni, VYi, VZi, MTi, MFYi, MFZi
    real(kind=8) :: Nj, VYj, VZj, MTj, MFYj, MFZj
    real(kind=8) :: Vx, VY, VZ, MT, MFY, MFZ

!   COMMAND DATA
    integer(kind=8) :: typcmb, typco, i, uc, um, typstru
    real(kind=8) :: cequi, sigs, sigci, sigcs, sigcyi, sigcys, sigczi, sigczs
    real(kind=8) :: alphacc, effrts(6), dnsits(6)
    real(kind=8) :: ht, bw, enrobi, enrobs, enrobyi, enrobys, enrobzi, enrobzs
    real(kind=8) :: gammac, gammas, rholmin, rhotmin, slsyme
    real(kind=8) :: facier, fbeton, eys, rhoacier, areinf, ashear
    real(kind=8) :: astirr, rhocrit, datcrit, lcrit, thiter, epiter, aphiter
    real(kind=8) :: wmaxi, wmaxs, wmaxyi, wmaxys, wmaxzi, wmaxzs, sigelsqp, kt
    real(kind=8) :: phixi, phixs, phiyi, phiys, phizi, phizs
    integer(kind=8) :: clacier, compress, epucisa, ferrcomp, ferrmin, ferrsyme, typdiag
    integer(kind=8) :: ierrl, ierrt, meth2D, cond109, precs
    character(len=24) :: valk(2)
    integer(kind=8) :: vali(2)

!   DEFAULT VALUES
    do i = 1, 6
        dnsits(i) = -1
    end do
    dnsvol = -1
    construc = -1

    call tecael(iadzi, iazk24, noms=1)
!
    meshName = zk24(iazk24-1+1)
    elementName = zk24(iazk24-1+3)
    node_nb = zi(iadzi-1+2)
    node_i_id = zi(iadzi-1+3)
    node_j_id = zi(iadzi-1+4)

!       -- RECUPERATION DES DONNEES DE L'UTILISATEUR :
!       ----------------------------------------------
!     FER1_R = 'TYPCOMB','CODIF','METH2D','THITER','EPITER','APHITER',
!                  1        2        3        4        5        6
!              'COND109','TYPSTRU','FERRSYME','SLSYME',
!                  7         8          9        10
!              'FERRCOMP','EPUCISA','FERRMIN','RHOLMIN','RHOTMIN',
!                  11        12         13        14       15
!              'COMPRESS','CEQUI','ENROBI','ENROBS','ENROBYI','ENROBYS',
!                  16        17      18       19       20        21
!              'ENROBZI','ENROBZS','SIGS','SIGCI','SIGCS','SIGCYI','SIGCYS',
!                  22        23      24      25      26      27       28
!              'SIGCZI','SIGCZS','ALPHACC','GAMMAS','GAMMAC','FACIER','EYS',
!                 29       30        31       32       33       34     35
!              'TYPDIAG','FBETON','CLACIER','UC','UM','RHOACIER','AREINF',
!                 36        37       38      39   40      41        42
!              'ASHEAR','ASTIRR','RHOCRIT','DATCRIT','LCRIT','WMAXI','WMAXS',
!                 43       44        45       46       47      48      49
!              'WMAXYI','WMAXYS','WMAXZI','WMAXZS','SIGELSQP','KT',
!                 50       51       52       53        54      55
!              'PHIXI','PHIXS','PHIYI','PHIYS','PHIZI','PHIZS'
!                 56      57      58      59      60      61
!
!                                  PCAGEPO
!  'HY1', 'HZ1', 'EPY1', 'EPZ1', 'HY2','HZ2', 'EPY2', 'EPZ2', 'R1', 'EP1',
!  'R2', 'EP2', 'TSEC',
!
!   TSEC = 0 GENERAL 1 RECTANGLE 2 CERCLE
!
!   Retriving instantied pointer for POUTRE GEOMETRICAL VALUES
    call jevech('PCAGEPO', 'L', jcagepo)
    HY1 = zr(jcagepo-1+1)
    HZ1 = zr(jcagepo-1+2)
    TSEC = zr(jcagepo-1+13)
    if (TSEC .ne. 1) then
        valk(1) = meshName
        valk(2) = elementName
        vali(1) = node_i_id
        vali(2) = node_j_id
        call utmess('F', 'CALCULEL7_1', nk=2, valk=valk, ni=2, vali=vali)
        goto 998
    end if

!   Retriving instantied pointer for INPUT PARAMETER
    call jevech('PFERRA1', 'L', jfer1)

    ht = HZ1
    bw = HY1
    typcmb = nint(zr(jfer1-1+1))
    typco = nint(zr(jfer1-1+2))
    meth2D = nint(zr(jfer1-1+3))
    thiter = zr(jfer1-1+4)
    epiter = zr(jfer1-1+5)
    aphiter = zr(jfer1-1+6)
    cond109 = nint(zr(jfer1-1+7))
    typstru = nint(zr(jfer1-1+8))
    ferrsyme = nint(zr(jfer1-1+9))
    slsyme = zr(jfer1-1+10)
    ferrcomp = nint(zr(jfer1-1+11))
    epucisa = nint(zr(jfer1-1+12))
    ferrmin = nint(zr(jfer1-1+13))
    rholmin = zr(jfer1-1+14)
    rhotmin = zr(jfer1-1+15)
    compress = int(zr(jfer1-1+16))
    cequi = zr(jfer1-1+17)
    enrobi = zr(jfer1-1+18)
    enrobs = zr(jfer1-1+19)
    enrobyi = zr(jfer1-1+20)
    enrobys = zr(jfer1-1+21)
    enrobzi = zr(jfer1-1+22)
    enrobzs = zr(jfer1-1+23)
    sigs = zr(jfer1-1+24)
    sigci = zr(jfer1-1+25)
    sigcs = zr(jfer1-1+26)
    sigcyi = zr(jfer1-1+27)
    sigcys = zr(jfer1-1+28)
    sigczi = zr(jfer1-1+29)
    sigczs = zr(jfer1-1+30)
    alphacc = zr(jfer1-1+31)
    gammas = zr(jfer1-1+32)
    gammac = zr(jfer1-1+33)
    facier = zr(jfer1-1+34)
    eys = zr(jfer1-1+35)
    typdiag = int(zr(jfer1-1+36))
    fbeton = zr(jfer1-1+37)
    clacier = int(zr(jfer1-1+38))
    uc = int(zr(jfer1-1+39))
    um = int(zr(jfer1-1+40))
    rhoacier = zr(jfer1-1+41)
    areinf = zr(jfer1-1+42)
    ashear = zr(jfer1-1+43)
    astirr = zr(jfer1-1+44)
    rhocrit = zr(jfer1-1+45)
    datcrit = zr(jfer1-1+46)
    lcrit = zr(jfer1-1+47)
    wmaxi = zr(jfer1-1+48)
    wmaxs = zr(jfer1-1+49)
    wmaxyi = zr(jfer1-1+50)
    wmaxys = zr(jfer1-1+51)
    wmaxzi = zr(jfer1-1+52)
    wmaxzs = zr(jfer1-1+53)
    sigelsqp = zr(jfer1-1+54)
    kt = zr(jfer1-1+55)
    phixi = zr(jfer1-1+56)
    phixs = zr(jfer1-1+57)
    phiyi = zr(jfer1-1+58)
    phiys = zr(jfer1-1+59)
    phizi = zr(jfer1-1+60)
    phizs = zr(jfer1-1+61)

    !Only option '1D'
    if (typstru .eq. 0.d0) then
        call utmess('A', 'CALCULEL_79')
        goto 998
    end if

    !Retriving instantied pointer for OUT PARAMETER
    call jevech('PFERRA2', 'E', jfer2)

    !Retriving instantied pointer for EFGE_ELNO
    call jevech('PEFFORR', 'L', jefge)
    !Retriving instantied pointer for EFGE_ELNO
    call tecach('OOO', 'PEFFORR', 'L', iret, nval=3, itab=itabin)
    !itabin(1): pointer of local field (in zr, zc,)
    !itabin(2): total length of local field
    !itabin(3): nb points (gauss or nodes)
    EFGE_ELNO_pointer = itabin(1)
    EFGE_ELNO_length = itabin(2)
    EFGE_ELNO_nbOfPoints = itabin(3)
    EFGE_ELNO_nbOfComponents = EFGE_ELNO_length/EFGE_ELNO_nbOfPoints

    Ni = zr(EFGE_ELNO_pointer-1+(1-1)*EFGE_ELNO_nbOfComponents+1)
    VYi = zr(EFGE_ELNO_pointer-1+(1-1)*EFGE_ELNO_nbOfComponents+2)
    VZi = zr(EFGE_ELNO_pointer-1+(1-1)*EFGE_ELNO_nbOfComponents+3)
    MTi = zr(EFGE_ELNO_pointer-1+(1-1)*EFGE_ELNO_nbOfComponents+4)
    MFYi = zr(EFGE_ELNO_pointer-1+(1-1)*EFGE_ELNO_nbOfComponents+5)
    MFZi = zr(EFGE_ELNO_pointer-1+(1-1)*EFGE_ELNO_nbOfComponents+6)

    call tecael(iadzi, iazk24, noms=1)

    Nj = zr(EFGE_ELNO_pointer-1+(2-1)*EFGE_ELNO_nbOfComponents+1)
    VYj = zr(EFGE_ELNO_pointer-1+(2-1)*EFGE_ELNO_nbOfComponents+2)
    VZj = zr(EFGE_ELNO_pointer-1+(2-1)*EFGE_ELNO_nbOfComponents+3)
    MTj = zr(EFGE_ELNO_pointer-1+(2-1)*EFGE_ELNO_nbOfComponents+4)
    MFYj = zr(EFGE_ELNO_pointer-1+(2-1)*EFGE_ELNO_nbOfComponents+5)
    MFZj = zr(EFGE_ELNO_pointer-1+(2-1)*EFGE_ELNO_nbOfComponents+6)

    !Mean value on nodes
    VX = (Ni+Nj)/2
    MFY = (MFYi+MFYj)/2
    MFZ = (MFZi+MFZj)/2
    VY = (VYi+VYj)/2
    VZ = (VZi+VZj)/2
    MT = (MTi+MTj)/2

    !Code Aster uses negative N for compression
    effrts(1) = -VX
    effrts(2) = -MFY
    effrts(3) = -MFZ
    effrts(4) = VY
    effrts(5) = VZ
    effrts(6) = MT

    !Calcul du ferraillage global de la poutre

    precs = ceiling(1/epiter)
    call glbpou(typcmb, typco, cequi, effrts, ht, bw, &
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

!       -- GESTION DES ALARMES EMISES :
!       -------------------------------
!
    if (ierrl .eq. 1001) then
!       ELU : section trop comprimée
        call utmess('A', 'CALCULEL7_12')
        dnsvol = -1.d0
        construc = -1.d0
    end if

    if (ierrl .eq. 10011) then
!       ELU : ferraillage symétrique non possible
        call utmess('A', 'CALCULEL7_28')
        dnsvol = -1.d0
        construc = -1.d0
    end if

    if (ierrl .eq. 10012) then
!       Resolution iterative Bresler non possible FCD
        call utmess('A', 'CALCULEL7_29')
        dnsvol = -1.d0
        construc = -1.d0
    end if
!
    if (ierrl .eq. 1003) then
!       ELS : section trop comprimée
        call utmess('A', 'CALCULEL7_13')
        dnsvol = -1.d0
        construc = -1.d0
    end if

    if (ierrl .eq. 1005) then
!       ELS_QP : section trop comprimée
        call utmess('A', 'CALCULEL7_14')
        dnsvol = -1.d0
        construc = -1.d0
    end if
!
    if (ierrt .eq. 1002) then
!       ELU BETON TROP CISAILLE : densité transversale fixée à -1 pour l'élément
        call utmess('A', 'CALCULEL7_15')
        dnsvol = -1.d0
        construc = -1.d0
    end if

    if (ierrt .eq. 1004) then
!       ELS BETON TROP CISAILLE : densité transversale fixée à -1 pour l'élément
        call utmess('A', 'CALCULEL7_16')
        dnsvol = -1.d0
        construc = -1.d0
    end if

    if (ierrt .eq. 1007) then
!       ELS_QP BETON TROP CISAILLE : densité transversale fixée à -1 pour l'élément
        call utmess('A', 'CALCULEL7_17')
        dnsvol = -1.d0
        construc = -1.d0
    end if

    if (ierrl .eq. 1006) then
!       ELS QP SOLLICITATION TROP IMPORTANTE : Résolution itérative impossible à l'els qp !
        call utmess('A', 'CALCULEL_77')
        dnsvol = -1.d0
        construc = -1.d0
    end if

998 continue

    zr(jfer2-1+1) = dnsits(1)
    zr(jfer2-1+2) = dnsits(2)
    zr(jfer2-1+3) = dnsits(3)
    zr(jfer2-1+4) = dnsits(4)
    zr(jfer2-1+5) = dnsits(5)
    zr(jfer2-1+6) = dnsits(6)
    zr(jfer2-1+7) = dnsvol
    zr(jfer2-1+8) = construc

end subroutine
