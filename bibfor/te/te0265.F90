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

subroutine te0265(nomopt, nomte)
    implicit none
#include "asterfort/assert.h"
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
!     HT         HAUTEUR DE LA POUTRE
!     BW         LARGEUR DE LA POUTRE
!     EFFRTS     TORSEUR DES EFFORTS ET DES MOMENTS (DIM 6)
!                     EFFRTS(1) = EFFN / EFFORT NORMAL
!                     EFFRTS(2) = EFFMY / MOMENT FL??CHISSANT SUIVANT Y
!                     EFFRTS(3) = EFFMZ / MOMENT FL??CHISSANT SUIVANT Z
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
    integer :: iret 
    ! Max is 7
    integer :: itabin(3) 
    integer :: EFGE_ELNO_pointer
    integer :: EFGE_ELNO_length
    integer :: EFGE_ELNO_nbOfPoints
    integer :: EFGE_ELNO_nbOfComponents
    integer :: jcagepo, jfer1, jfer2, jefge
    
!   INFO ON MAILLE AND ELEMENT
    ! For tecael
    integer :: iadzi, iazk24 
    character(len=24) :: meshName, elementName
!   For tecael
    integer :: node_nb, node_i_id, node_j_id 
    
!   OUTPUT FOR 1D ELEMENTS
    real(kind=8) :: dnsvol, construc
    
!   GEOMETRICAL DATA
    real(kind=8) :: HY1, HZ1, TSEC
!   EFGE_ELNO DATA
    real(kind=8) :: Ni, VYi, VZi, MTi, MFYi, MFZi
    real(kind=8) :: Nj, VYj, VZj, MTj, MFYj, MFZj
    real(kind=8) :: Vx, VY, VZ, MT, MFY, MFZ

!   COMMAND DATA
    integer :: typcmb, typco, i, k, uc, um, typstru
    real(kind=8) :: cequi, sigs, sigci, sigcs, sigcyi, sigcys, sigczi, sigczs
    real(kind=8) :: alphacc, effrts(6), dnsits(6)
    real(kind=8) :: ht, bw, enrobi, enrobs, enrobyi, enrobys, enrobzi, enrobzs
    real(kind=8) :: gammac, gammas, rholmin, rhotmin, slsyme
    real(kind=8) :: facier, fbeton, eys, rhoacier, areinf, ashear
    real(kind=8) :: astirr, rhocrit, datcrit, lcrit
    real(kind=8) :: wmaxi, wmaxs, wmaxyi, wmaxys, wmaxzi, wmaxzs, sigelsqp, kt
    real(kind=8) :: phixi, phixs, phiyi, phiys, phizi, phizs
    integer :: clacier, compress, epucisa, ferrcomp, ferrmin, ferrsyme, typdiag
    integer :: ierr
    character(len=24) :: valk(2)
    integer :: vali(2)

    
!   DEFAULT VALUES
    do i=1,6
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
    if (TSEC.ne.1) then
    valk(1) = meshName
    valk(2) = elementName
    vali(1) = node_i_id
    vali(2) = node_j_id
    call utmess('F', 'CALCULEL7_1', nk=2, valk=valk, ni=2, vali=vali)
    goto 998
    endif

!   Retriving instantied pointer for INPUT PARAMETER
    call jevech('PFERRA1', 'L', jfer1)
    
    ht = HZ1
    bw = HY1
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

   !Only option '1D'
    if (typstru.eq.0.d0) then
    call utmess('A', 'CALCULEL_79')
    goto 998
    endif
    
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

    Ni = zr(EFGE_ELNO_pointer-1+(1-1)*EFGE_ELNO_nbOfComponents + 1)
    VYi = zr(EFGE_ELNO_pointer-1+(1-1)*EFGE_ELNO_nbOfComponents + 2)
    VZi = zr(EFGE_ELNO_pointer-1+(1-1)*EFGE_ELNO_nbOfComponents + 3)
    MTi = zr(EFGE_ELNO_pointer-1+(1-1)*EFGE_ELNO_nbOfComponents + 4)
    MFYi = zr(EFGE_ELNO_pointer-1+(1-1)*EFGE_ELNO_nbOfComponents + 5)
    MFZi = zr(EFGE_ELNO_pointer-1+(1-1)*EFGE_ELNO_nbOfComponents + 6)

    call tecael(iadzi, iazk24, noms=1)

    Nj   = zr(EFGE_ELNO_pointer-1+(2-1)*EFGE_ELNO_nbOfComponents + 1)
    VYj  = zr(EFGE_ELNO_pointer-1+(2-1)*EFGE_ELNO_nbOfComponents + 2)
    VZj  = zr(EFGE_ELNO_pointer-1+(2-1)*EFGE_ELNO_nbOfComponents + 3)
    MTj  = zr(EFGE_ELNO_pointer-1+(2-1)*EFGE_ELNO_nbOfComponents + 4)
    MFYj = zr(EFGE_ELNO_pointer-1+(2-1)*EFGE_ELNO_nbOfComponents + 5)
    MFZj = zr(EFGE_ELNO_pointer-1+(2-1)*EFGE_ELNO_nbOfComponents + 6)

   !Mean value on nodes
    VX = (Ni+Nj)/2
    MFY = (MFYi+MFYj)/2
    MFZ = (MFZi+MFZj)/2
    VY = (VYi+VYj)/2
    VZ = (VZi+VZj)/2
    MT = (MTi+MTj)/2
    
   !Code Aster uses negative N for compression
    effrts(1) = -VX
    effrts(2) = MFY
    effrts(3) = MFZ
    effrts(4) = VY
    effrts(5) = VZ
    effrts(6) = MT

   !Calcul du ferraillage global de la poutre
   
    call glbpou(typcmb, typco, cequi, effrts, ht, bw,&
                  enrobyi, enrobys, enrobzi, enrobzs,&
                  facier, fbeton, sigelsqp, kt, eys,&
                  alphacc, clacier, gammas, gammac, typdiag,&
                  sigcyi, sigcys, sigczi, sigczs, sigs,&
                  wmaxyi, wmaxys, wmaxzi, wmaxzs,&
                  phiyi, phiys, phizi, phizs,&
                  ferrsyme, slsyme, ferrcomp,&
                  epucisa, ferrmin, rholmin, rhotmin, compress, uc, um, &
                  rhoacier, areinf, ashear, astirr, rhocrit, datcrit, lcrit,&
                  dnsits, dnsvol, construc, ierr)
    
    
!       -- GESTION DES ALARMES EMISES :
!       -------------------------------
!
    if (ierr.eq.1001) then
!       ELU : section trop comprim??e
!       on fixe toutes les densit??s de ferraillage de l'??l??ment ?? -1
        call utmess('A', 'CALCULEL7_12')
        do k = 1, 6
            dnsits(k) = -1.d0
        end do
        dnsvol = -1.d0
        construc = -1.d0
    endif
    
    if (ierr.eq.10011) then
!       ELU : ferraillage sym??trique non possible
!       on fixe toutes les densit??s de ferraillage de l'??l??ment ?? -1
        call utmess('A', 'CALCULEL7_28')
        do k = 1, 6
            dnsits(k) = -1.d0
        end do
        dnsvol = -1.d0
        construc = -1.d0
    endif
    
    if (ierr.eq.10012) then
!       Resolution iterative Bresler non possible FCD
!       on fixe toutes les densit??s de ferraillage de l'??l??ment ?? -1
        call utmess('A', 'CALCULEL7_29')
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
        call utmess('A', 'CALCULEL7_13')
        do k = 1, 6
            dnsits(k) = -1.d0
        end do
        dnsvol = -1.d0
        construc = -1.d0
    endif
    
    if (ierr.eq.1005) then
!       ELS_QP : section trop comprim??e
!       on fixe toutes les densit??s de ferraillage de l'??l??ment ?? -1
        call utmess('A', 'CALCULEL7_14')
        do k = 1, 6
            dnsits(k) = -1.d0
        end do
        dnsvol = -1.d0
        construc = -1.d0
    endif
!
    if (ierr.eq.1002) then
!       ELU BETON TROP CISAILLE : densit?? transversale fix??e ?? -1 pour l'??l??ment
        call utmess('A', 'CALCULEL7_15')
        dnsits(5) = -1.d0
        dnsvol = -1.d0
        construc = -1.d0
    endif

    if (ierr.eq.1004) then
!       ELS BETON TROP CISAILLE : densit?? transversale fix??e ?? -1 pour l'??l??ment
        call utmess('A', 'CALCULEL7_16')
        dnsits(5) = -1.d0
        dnsvol = -1.d0
        construc = -1.d0
    endif

    if (ierr.eq.1007) then
!       ELS_QP BETON TROP CISAILLE : densit?? transversale fix??e ?? -1 pour l'??l??ment
        call utmess('A', 'CALCULEL7_17')
        dnsits(5) = -1.d0
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
