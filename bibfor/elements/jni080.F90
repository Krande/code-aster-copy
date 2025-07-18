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
! aslint: disable=W1501
subroutine jni080(elrefe, nmaxob, liobj, nbobj)
!
    implicit none
!
#include "jeveux.h"
#include "asterc/r8nnem.h"
#include "asterfort/assert.h"
#include "asterfort/elraga.h"
#include "asterfort/elrfdf.h"
#include "asterfort/elrfvf.h"
#include "asterfort/fcepai.h"
#include "asterfort/fcesnd.h"
#include "asterfort/jeexin.h"
#include "asterfort/mamagi.h"
#include "asterfort/wkvect.h"
!
    character(len=8) :: elrefe
    integer(kind=8) :: nmaxob, nbobj
!
!.......................................................................
!
! BUT :  ROUTINE D'INITIALISATION DES ELEMENTS COQUE (MEC3QU9H,MEC3TR7H)
!.......................................................................
!
!    ON CONSTRUIT 3 TABLEAUX POUR CHAQUE ELEMENT (MEC3QU9H,MEC3TR7H)
!
!    DESI='&INEL.'//ELREFE//'.DESI' CONTIENT :
!                MEC3QU9H    MEC3TR7H
!    NBN1        8           6
!    NBN2        9           7
!    NPGSR       4           3
!    NPGSN       9           7
!
!   DESR='&INEL.'//ELREFE//'.DESR' (LONGUEUR 2100) CONTIENT :
! ======================================================================
! OBJET  DESR                     NBVAL    POINTEUR   NBVAL    POINTEUR
!                               MEC3QU9H   Q9        MEC3TR7H     T7
! ======================================================================
! INTEGRATION REDUITE
! COORDONNEE X POINTS DE GAUSS      4                 3
! COORDONNEE Y POINTS DE GAUSS      4                 3
! POIDS DE GAUSS                    4                 3
! ----------------------------------------------------------------------
! FONCTIONS DE FORME SERENDIP (NBN1) PT GAUSS REDUITS (NPSGR)
! POINTEURS CI-DESSOUS UTILISEES DANS : BTDMSR, JM1DN1, JM1DN2 (VALFOR),
! VECTGT,
!                                          12                  12
! ----------------------------------------------------------------------
! NI (PG1) I=1,NBN1                 8                 6
! ...                               ...               ...
! NI (PGN) I=1,NBN1                 8                 6
!                                         44                    44
! DNI/DKSI I=1,NBN1                 8                 6
! ...                               ...               ...
! DNI/DKSI I=1,NBN1                 8                 6
!                                         76                    76
! DNI/DETA I=1,NBN1                 8                 6
! ...                               ...               ...
! DNI/DETA I=1,NBN1                 8                 6
!      NOMBRE TOTAL 3*NPGSR*NBN1   96     108        54         63
! ----------------------------------------------------------------------
! DEFINITION POINTS DE GAUSS COMPLETS (NPGSN)
! POINTEURS CI-DESSOUS UTILISEES DANS : FORNGR, FORCEN, FORPES, FORSRG
! HSJ1F TE0402 TE0406 TE0417 TE0486 VDGNLR, VDPNLR
!                                         108                  108
! ----------------------------------------------------------------------
! INTEGRATION COMPLETE
! COORDONNEE X POINTS DE GAUSS      9                 7
! COORDONNEE Y POINTS DE GAUSS      9                 7
!                                         126                  126
! POIDS DE GAUSS                    9                 7
!                                  27     135        21        130
! ----------------------------------------------------------------------
! FONCTIONS DE FORME SERENDIP (NBN1) PT GAUSS COMPLETS (NPSGN)
! POINTEURS CI-DESSOUS UTILISEES DANS : FORCEN, FORPES, FORSRG, MATRN
!  TE0417, JM1DN1, JM1DN2 (VALFOR), VECTGT, VECTCI
!                                         135                  135
! ----------------------------------------------------------------------
! NI (PG1) I=1,NBN1                 8                 6
! ...                               ...               ...
! NI (PGN) I=1,NBN1                 8                 6
!                                         207                  207
! DNI/DKSI I=1,NBN1                 8                 6
! ...                               ...               ...
! DNI/DKSI I=1,NBN1                 8                 6
!                                         279                  279
! DNI/DETA I=1,NBN1                 8                 6
! ...                               ...               ...
! DNI/DETA I=1,NBN1                 8                 6
!      NOMBRE TOTAL 3*NPGSN*NBN1    216   351        126       261
! ----------------------------------------------------------------------
! FONCTIONS DE FORME LAGRANGE (NBN2) PT GAUSS REDUITS (NPSGR)
! POINTEURS CI-DESSOUS UTILISEES DANS : JM1DN1, JM1DN2 (VALFOR)
!                                         351                  351
! ----------------------------------------------------------------------
! NI (PG1) I=1,NBN2                 9                 7
! ...                               ...               ...
! NI (PGN) I=1,NBN2                 9                 7
!                                         387                  387
! DNI/DKSI I=1,NBN2                 9                 7
! ...                               ...               ...
! DNI/DKSI I=1,NBN2                 9                 7
!                                         423                  423
! DNI/DETA I=1,NBN2                 9                 7
! ...                               ...               ...
! DNI/DETA I=1,NBN2                 9                 7
!  NOMBRE TOTAL 3*NPGSR*NBN2      108    459         63        414
! ----------------------------------------------------------------------
! FONCTIONS DE FORME LAGRANGE (NBN2) PT GAUSS COMPLETS (NPSGN)
! POINTEURS CI-DESSOUS UTILISEES DANS : FORSRG, TE0486,
! JM1DN1, JM1DN2 (VALFOR), VDXNLR, MATRN, JM1DN3, VDGNLLR, VDPNLR
!                                        459                  459
! ----------------------------------------------------------------------
! NI (PG1) I=1,NBN2                 9                 7
! ...                               ...               ...
! NI (PGN) I=1,NBN2                 9  NPGSN*NBN2
!                                         540                  540
! DNI/DKSI I=1,NBN2                 9
! ...                               ...               ...
! DNI/DKSI I=1,NBN2                 9  NPGSN*NBN2
!                                         621                  621
! DNI/DETA I=1,NBN2                 9
! ...                               ...               ...
! DNI/DETA I=1,NBN2                 9  NPGSN*NBN2
!   NOMBRE TOTAL 3*NPGSN*NBN2      243     702        147      606
! ----------------------------------------------------------------------
! FONCTIONS DE FORME REDUITES (LINEAIRES SUR ELEMENT DONT LES SOMMETS
! SONT LES POINTS DE GAUSS REDUITS) EVALUEES AUX PT GAUSS COMPLETS
! POINTEURS CI-DESSOUS UTILISEES DANS :  MATBSU, BTDMSN, VDGNLR, VDPNLR
!                                          702                  702
! ----------------------------------------------------------------------
! NI (PG1) I=1,NPGSR                4
! ...                               ...               ...
! NI (PGN) I=1,NPGSR                4
! DNI/DKSI I=1,NPGSR                 4
! ...                               ...               ...
! DNI/DKSI I=1,NPGSR                 4
! DNI/DETA I=1,NPGSR                 4
! ...                               ...               ...
! DNI/DETA I=1,NPGSR                 4
!   NOMBRE TOTAL 3*NPGSN*NBN2      108     810
! ----------------------------------------------------------------------
! COORDONNEES DES NBN2 NOEUDS
! POINTEURS CI-DESSOUS UTILISEES DANS :
!                                          810                 810
! ----------------------------------------------------------------------
! COORDONNES DES NOEUDS X           9                  7
! COORDONNES DES NOEUDS Y           9                  7
!   NOMBRE TOTAL 3*NPGSN*NBN2      18     828          14      824
! ----------------------------------------------------------------------
! DERIVEES DES NBN1 FONCTIONS DE FORMES AUX NBN2 NOEUDS
! POINTEURS CI-DESSOUS UTILISEES DANS : VECTAN
!                                         828                 828
! ----------------------------------------------------------------------
! DNI/DKSI I=1,NBN2                 9                 7
! ...                               ...               ...
! DNI/DKSI I=1,NBN2                 9                 7
!                                         900                  900
! DNI/DETA I=1,NBN2                 9                 7
! ...                               ...               ...
! DNI/DETA I=1,NBN2                 9                 7
!   NOMBRE TOTAL 2*NBN1 *NBN2     144     972        84        912
! ----------------------------------------------------------------------
! ZONE DE TRAVAIL LONGUEUR 10*NB1 POINTEUR             1000
! UTILISEE DANS FPRES, FSURF
!                                                      1080
! ----------------------------------------------------------------------
! ZONE DE TRAVAIL LONGUEUR 9*NB2 POINTEUR             1090
! UTILISEE DANS VDREPE, VECTAN
! ----------------------------------------------------------------------
! ZONE DE TRAVAIL LONGUEUR 9     POINTEUR             1180
! UTILISEE DANS HSJ1F
! ----------------------------------------------------------------------
! ESPACE VIDE DANS DESR JUSQU'A :        1250                 1250
! ----------------------------------------------------------------------
! COORDONNES DES 3 POINTS DANS L'EPAISSEUR  1253              1253
! ----------------------------------------------------------------------
! POUR CHAQUE POINT DANS L'EPAISSEUR
! VALEURS DES NPGSN-1 FONCTIONS D'INTERPOLATION AUX NBN1 NOEUDS
! CES FONCTIONS SONT DEFINIES SUR UN ELEMENT DE REFERENCE
! CONSTRUIT SUR LES NPGSN -1 POINTS DE GAUSS
! POINTEURS CI-DESSOUS UTILISEES DANS : VDESNG
!                                         1260                 1260
! ----------------------------------------------------------------------
!  NOMBRE TOTAL 3*NBN1*(NPGSN-1)  192     1452     108        1368
! ----------------------------------------------------------------------
! VALEURS DES NPGSR FONCTIONS D'INTERPOLATION AUX NBN1 NOEUDS
! CES FONCTIONS SONT DEFINIES SUR UN ELEMENT DE REFERENCE
! CONSTRUIT SUR LES NPGSR POINTS DE GAUSS      NBN1*(NPGSR)
! POINTEURS CI-DESSOUS UTILISEES DANS : VDEFGE
!                                        1452                 1452
! ----------------------------------------------------------------------
!  NOMBRE TOTAL 3*NBN1*(NPGSR)    32     1484     18          1470
! ----------------------------------------------------------------------
! DANS FCEPAIS, POIDS ET POINTS DE GAUSS DANS L'EPAISSEUR
! POINTEURS CI-DESSOUS UTILISEES DANS :
!                                        1500                 1500
! ----------------------------------------------------------------------
! COORDONNES ET POIDS                     6                    6
! 3 VALEURS POUR CHACUN DES 3 POINTS      9                    9
! ESPACE VIDE DANS DESR A PARTIR DE      1515                 1515
! ----------------------------------------------------------------------
! ZONE DE TRAVAIL LONGUEUR 1 POINTEUR             1550
! UTILISEE DANS CAURTG,FORNRG, PK2GAU,VDGNLR, VDPNLR, VDXRIG, VDXNLR
! ----------------------------------------------------------------------
! ZONE DE TRAVAIL LONGUEUR 9*NPGSR  POINTEUR        2000
! UTILISEE DANS TE0415, VDXSIG
! ----------------------------------------------------------------------
! ESPACE VIDE DANS DESR JUSQU'A :        2100                 2100
! ----------------------------------------------------------------------
!
! ======================================================================
! LISTE DES ROUTINES COQUE3D UTILISANT LES .DESR
!=======================================================================
!
! BSTHCO UTILISE ZR(459)
! BTDFN  IND=1 L=459 + IND=0 L=351 FF LAGRANGE PG REDUIT
! BTDMSN INTERGRATION REDUITE L=702 FF REDUITES SUR PT GAUSS NORMAUX
    character(len=24) :: liobj(nmaxob)
! BTDMSR DERIVEES FF SERENDIP PG REDUITS L=44,76 FF ET
! BTDMSR DERIVEES LAGRANGE PG REDUITS L=351, 387, 423
! BTLDTH INDIC=0 L=351 FF LAGRANGE PG REDUITS
! BTLDTH INDIC=1 L=459 FF LAGRANGE PG COMPLETS
! CAURTG ZONE TRAVAIL ZR(1550) L=1
! FORCEN L=127 POIDS PG COMPLETS  L=135 FF SERENDIP PG COMPLETS
! FORNGR ZONE TRAVAIL ZR(1550) BTSIG(ZR(127...=POIDS)
! FORPES L=127 POIDS PG COMPLETS  L=135 FF SERENDIP PG COMPLETS
! FORSRG L=127 POIDS PG COMPLETS  L=135 FF SERENDIP PG COMPLETS
! FORSRG L=459 FF LAGRANGE PG COMPLET
! FPRES  ZONE TRAVAIL DE 1000 A 1000+10*NB1 = 1080
! FSURF  ZONE TRAVAIL DE 1000 A 1000+10*NB1 1080
! HSJ1F  L=127 POIDS PG COMPLETS ZONE TRAVAIL XR(1180) LONGUEUR 9
! JM1DN1 UTILISE VALFOR
! JM1DN2   UTILISE VALFOR
! JM1DN3  L=459, 540, 621 FF ET DERIVEES LAGRANGE PG COMPLETS
! MATBSU L=702 L=702 FF REDUITES SUR PT GAUSS NORMAUX
! MATRN  L=135 FF SERENDIP PG COMPLET + L=459 FF LAGRANGE PG COMPLET
! PK2CAU ZONE TRAVAIL ZR(1550)
! TE0402 RIGI_GEOM UTILISE ZR(127...=POIDS) BTSIG(ZR(127...=POIDS)
! TE0406 INTEGRATION NUMERIQUE ZR(127)
! TE0415 SIEF_ELNO ZONE TRAVAIL ZR(2000) LONGUEUR 9*NBGSR
! TE0417 UTILISE ZR(127...) ZR(135..)
! TE0419 CHAR_MECA_TEMP_R
! TE0486 UTILISE ZR(459) B1TDB2(ZR(127)
! VALFOR INDN=0 => LT1=44  LT2=76 DERIVEES DES FF SERENDIP PG RESDUITS
! VALFOR INDN=1 => LT1=207 LT2=279 DERIVEES DES FF SERENDIP PG COMPLETS
! VALFOR L1=351 L2=387 L3=423, FF ET DERIVEES LAGRANGE PG REDUITS
! VALFOR L1=459 L2=540 L3=621, FF ET DERIVEES LAGRANGE PG COMPLETS
! VDEFGE L=1452 FF BASEES SUR NPGSR EVALUEES AUX NPG1 NOEUDS
! VDESND Q9 L=1260 FF BASEES SUR NPGSN EVALUEES AUX NPG1 NOEUDS
! VDESND T7 L=1260 FF BASEES SUR NPGSN EVALUEES AUX NPG1 NOEUDS
! VDGNLR ZONE TRAVAIL ZR(1550)  ZR(459) BTSIG(ZR(127)) ZR(702)
! VDPNLR ZONE TRAVAIL ZR(1550)  ZR(459) BTSIG(ZR(127)) ZR(702)
! VDREPE ZONE TRAVAIL ZR(1090) MATRICE PASSAGE  LONGUEUR 9*NB2
! VDXNLR UTILISE ZR(459) ZR(LZR-1+1550) = COEF
! VDXRIG ZR(LZR-1+1550) = COEF
! VDXSIG ZONE TRAVAIL  ZR(2000) LONGUEUR 9*NPGSR
! VECTAN L=1090 ZONE TRAVAIL  MATRICES DE PASSAGE LONGUEUR 9*NB2
! VECTAN L=828, 900 DERIVEES DES NBN1 FF EVALUEES AUX NBN2 NOEUDS
! VECTCI L=207, 279 DERIVEESS FF SERENDIP PG COMPLETS
! VECTGT IND=0 L=12, 44, 76  FF SERENDIP ET DERIVESS PG REDUITS
! VECTGT IND=1 L=135, 207, 279  FF SERENDIP ET DERIVESS PG COMPLETS
! ======================================================================
!
! CHANTIER ELREFE : QUE FAIRE ?
!  - IL FAUDRAIT VIRER FCEPAI, QUI APPAREMMENT N'EST PAS UTILISE
!  - SI ON VEUT DIMINUER LE NOMBRE DE FF, SUPPRIMER LES FF 1260 ET 1452
!    REMPLACER DANS VDEFGE ET VDESNR L'EXTRAPOLATION AUX NOEUDS ACTUELLE
!  PAR UN APPEL A PPGANO. LE PB C'EST QUE LA MATRICE DE PASSAGE ACTUELLE
!  EST CONSTRUITE SUR LES 3 POINTS D'INTéGRATION PAR COUCHE UTILISéS EN
!  NON LINéAIRE, ALORS QU'EN LINéAIRE, IL N'Y A QUE 2 POINTS.
!  IL FAUDRAIT DONC SOIT METTRE 3 POINTS EN LINéAIRE, SOIT CHANGER LE
!  PASSAGE GAUSS NOEUDS EN NON LINéAIRE ET EN LINéAIRE.
!  - IL FAUDRAIT REMPLACER LES ZONES DE TRAVAIL PAR DES PASSAGES
!  D'ARGUMENTS
!  - ENSUITE ON POURRAIT PEUT ETE UTILISER DES ELREFE...
! ======================================================================
    character(len=16) :: elrefl
    character(len=24) :: desi, desr
    integer(kind=8) :: iret, ndim, nbpg, npgcou, npgsn, npgsr, nso, ino
    integer(kind=8) :: lzi, lzr, i1, i2, i3, i4, i5, k, l, ll, m, nbn1, nbn2, kompt
    integer(kind=8) :: jmas, ldesi, ldesr, ljmas
    real(kind=8) :: a, b, aa
    real(kind=8) :: vfesnd(45)
    real(kind=8) :: xpg(81), poipg(27), x(2), ff(9), dff(3, 9), xi1, xi2, xi3
!
! DEB ------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: l1
!-----------------------------------------------------------------------
    nbobj = 3
    ASSERT(nmaxob .gt. nbobj)
    liobj(1) = '&INEL.'//elrefe//'.DESI'
    liobj(2) = '&INEL.'//elrefe//'.DESR'
    elrefl = elrefe
    liobj(3) = '&INEL.'//elrefl//'.B'
!
    desi = liobj(1)
    desr = liobj(2)
!
    npgcou = 3
    if (elrefe .eq. 'MEC3QU9H') then
        nbn1 = 8
        nbn2 = 9
        npgsr = 4
        npgsn = 9
        nso = 4
    else if (elrefe .eq. 'MEC3TR7H') then
        nbn1 = 6
        nbn2 = 7
        npgsr = 3
        npgsn = 7
        nso = 3
    else
        ASSERT(ASTER_FALSE)
    end if
!
    call jeexin(desi, iret)
    if (iret .ne. 0) goto 999
!
    ldesi = 4
    call wkvect(desi, 'V V I', ldesi, lzi)
    zi(lzi-1+1) = nbn1
    zi(lzi-1+2) = nbn2
    zi(lzi-1+3) = npgsr
    zi(lzi-1+4) = npgsn
!
    ldesr = 2100
    call wkvect(desr, 'V V R', ldesr, lzr)
    zr(lzr-1+1550) = r8nnem()
!
    if (elrefe .eq. 'MEC3QU9H') then
!
!     DEFINITION DES 4=2*2 PTS DE GAUSS REDUIT ET DES POIDS
!     CORRESPONDANTS
        call elraga('QU9', 'FPG4    ', ndim, nbpg, xpg, poipg)
!
        zr(lzr-1+1) = xpg(1)
        zr(lzr-1+2) = xpg(3)
        zr(lzr-1+3) = xpg(5)
        zr(lzr-1+4) = xpg(7)
!
        zr(lzr-1+5) = xpg(2)
        zr(lzr-1+6) = xpg(4)
        zr(lzr-1+7) = xpg(6)
        zr(lzr-1+8) = xpg(8)
!
        zr(lzr-1+9) = poipg(1)
        zr(lzr-1+10) = poipg(2)
        zr(lzr-1+11) = poipg(3)
        zr(lzr-1+12) = poipg(4)
!
!     FONCTIONS SERENDIP - POUR LES TERMES DE
!     TRANSLATION
!
!     VALEURS DES FONCTIONS DE FORME ET DE
!     LEURS DERIVEES AUX 4 PTS DE GAUSS REDUITS
!
        do l = 1, 4
            i1 = l
            i2 = 4+l
            x(1) = zr(lzr-1+i1)
            x(2) = zr(lzr-1+i2)
            call elrfvf('QU8', x, ff)
            call elrfdf('QU8', x, dff)
            ll = 8*(l-1)
            do ino = 1, 8
                zr(lzr-1+12+ll+ino) = ff(ino)
                zr(lzr-1+44+ll+ino) = dff(1, ino)
                zr(lzr-1+76+ll+ino) = dff(2, ino)
            end do
        end do
!
!     DEFINITION DES 9=3*3 PTS DE GAUSS NORMAL
!     POSITION DES POINTS DE GAUSS DANS LE PLAN
!     ET VALEUR DES POIDS ASSOCIES
        call elraga('QU8', 'FPG9COQ ', ndim, nbpg, xpg, &
                    poipg)
!
        zr(lzr-1+109) = xpg(1)
        zr(lzr-1+110) = xpg(3)
        zr(lzr-1+111) = xpg(5)
        zr(lzr-1+112) = xpg(7)
        zr(lzr-1+113) = xpg(9)
        zr(lzr-1+114) = xpg(11)
        zr(lzr-1+115) = xpg(13)
        zr(lzr-1+116) = xpg(15)
        zr(lzr-1+117) = xpg(17)

        zr(lzr-1+118) = xpg(2)
        zr(lzr-1+119) = xpg(4)
        zr(lzr-1+120) = xpg(6)
        zr(lzr-1+121) = xpg(8)
        zr(lzr-1+122) = xpg(10)
        zr(lzr-1+123) = xpg(12)
        zr(lzr-1+124) = xpg(14)
        zr(lzr-1+125) = xpg(16)
        zr(lzr-1+126) = xpg(18)
!
!     VALEUR DES POIDS ASSOCIES
!
        zr(lzr-1+127) = poipg(1)
        zr(lzr-1+128) = poipg(2)
        zr(lzr-1+129) = poipg(3)
        zr(lzr-1+130) = poipg(4)
        zr(lzr-1+131) = poipg(5)
        zr(lzr-1+132) = poipg(6)
        zr(lzr-1+133) = poipg(7)
        zr(lzr-1+134) = poipg(8)
        zr(lzr-1+135) = poipg(9)
!
!     VALEURS DES FONCTIONS DE FORME ET
!     DERVIVEES AUX 9 PTS DE GAUSS NORMAL
!
        do l = 1, 9
            i1 = 108+l
            i2 = 108+9+l
            x(1) = zr(lzr-1+i1)
            x(2) = zr(lzr-1+i2)
            call elrfvf('QU8', x, ff)
            call elrfdf('QU8', x, dff)
            ll = 8*(l-1)
            do ino = 1, 8
                zr(lzr-1+135+ll+ino) = ff(ino)
                zr(lzr-1+207+ll+ino) = dff(1, ino)
                zr(lzr-1+279+ll+ino) = dff(2, ino)
            end do
        end do
!
!     FONCTIONS DE LAGRANGE - POUR LES TERMES
!     DE ROTATION
!
!     VALEURS AUX 4 PTS DE GAUSS REDUITS
!     DES FONCTIONS DE FORME ET DE LEURS
!     DERIVEES
!
        do l = 1, 4
            i1 = l
            i2 = 4+l
            x(1) = zr(lzr-1+i1)
            x(2) = zr(lzr-1+i2)
            call elrfvf('QU9', x, ff)
            call elrfdf('QU9', x, dff)
            ll = 9*(l-1)
            do ino = 1, 9
                zr(lzr-1+351+ll+ino) = ff(ino)
                zr(lzr-1+387+ll+ino) = dff(1, ino)
                zr(lzr-1+423+ll+ino) = dff(2, ino)
            end do
        end do
!
!    VALEURS AUX 9 PTS DE GAUSS NORMAL
!    DES FONCTIONS DE FORME ET DE LEURS
!    DERIVEES
!
        do l = 1, 9
            i1 = 108+l
            i2 = 108+9+l
            x(1) = zr(lzr-1+i1)
            x(2) = zr(lzr-1+i2)
            call elrfvf('QU9', x, ff)
            call elrfdf('QU9', x, dff)
            ll = 9*(l-1)
            do ino = 1, 9
                zr(lzr-1+459+ll+ino) = ff(ino)
                zr(lzr-1+540+ll+ino) = dff(1, ino)
                zr(lzr-1+621+ll+ino) = dff(2, ino)
            end do
        end do
!
!    FONCTIONS ASSOCIEES AUX 4 PTS DE GAUSS REDUITS
!
!    VALEURS AUX 9 PTS DE GAUSS NORMAL
!
        do l = 1, 9
            i1 = 108+l
            i2 = 108+9+l
            x(1) = zr(lzr-1+i1)
            x(2) = zr(lzr-1+i2)
            aa = 0.577350269189626d0
            ff(1) = 0.75d0*(aa-x(1))*(aa-x(2))
            ff(2) = 0.75d0*(aa+x(1))*(aa-x(2))
            ff(3) = 0.75d0*(aa+x(1))*(aa+x(2))
            ff(4) = 0.75d0*(aa-x(1))*(aa+x(2))
            dff(1, 1) = -0.75d0*(aa-x(2))
            dff(1, 2) = 0.75d0*(aa-x(2))
            dff(1, 3) = 0.75d0*(aa+x(2))
            dff(1, 4) = -0.75d0*(aa+x(2))
            dff(2, 1) = -0.75d0*(aa-x(1))
            dff(2, 2) = -0.75d0*(aa+x(1))
            dff(2, 3) = 0.75d0*(aa+x(1))
            dff(2, 4) = 0.75d0*(aa-x(1))
            ll = 4*(l-1)
            do l1 = 1, 4
                i3 = 702+ll+l1
                i4 = 738+ll+l1
                i5 = 774+ll+l1
                zr(lzr-1+i3) = ff(l1)
                zr(lzr-1+i4) = dff(1, l1)
                zr(lzr-1+i5) = dff(2, l1)
            end do
        end do
!
!     VALEURS DES FONCTIONS DE SERENDIP AUX 9 NOEUDS
!     ON DONNE LA POSITION DES NOEUDS
!
        zr(lzr-1+811) = -1.d0
        zr(lzr-1+812) = 1.d0
        zr(lzr-1+813) = 1.d0
        zr(lzr-1+814) = -1.d0
        zr(lzr-1+815) = 0.d0
        zr(lzr-1+816) = 1.d0
        zr(lzr-1+817) = 0.d0
        zr(lzr-1+818) = -1.d0
        zr(lzr-1+819) = 0.d0

        zr(lzr-1+820) = -1.d0
        zr(lzr-1+821) = -1.d0
        zr(lzr-1+822) = 1.d0
        zr(lzr-1+823) = 1.d0
        zr(lzr-1+824) = -1.d0
        zr(lzr-1+825) = 0.d0
        zr(lzr-1+826) = 1.d0
        zr(lzr-1+827) = 0.d0
        zr(lzr-1+828) = 0.d0
!
!     VALEURS DES FONCTIONS DE SERENDIP
!     ET DE LEURS DERIVEES AUX 9 NOEUDS
!
        do l = 1, 9
            i1 = 810+l
            i2 = 810+9+l
            x(1) = zr(lzr-1+i1)
            x(2) = zr(lzr-1+i2)
            call elrfdf('QU8', x, dff)
            ll = 8*(l-1)
            do ino = 1, 8
                zr(lzr-1+828+ll+ino) = dff(1, ino)
                zr(lzr-1+900+ll+ino) = dff(2, ino)
            end do
        end do
!
!     DEFINITION DES 8 FONCTIONS D'INTERPOLATION QUI PERMETTENT
!     D'EXTRAPOLER LES DEFORMATIONS OU CONTRAINTES AUX NOEUDS
!     PARTIR DE LEURS VALEURS AUX POINTS DE GAUSS
!
        zr(lzr-1+1251) = -1.d0
        zr(lzr-1+1252) = 0.d0
        zr(lzr-1+1253) = 1.d0
!
        kompt = 0
!
        do m = 1, 3
            i3 = 1250+m
            xi3 = zr(lzr-1+i3)
            do l = 1, 8
                i1 = 810+l
                i2 = 810+9+l
                xi1 = zr(lzr-1+i1)
                xi2 = zr(lzr-1+i2)
                kompt = kompt+1
                ll = 8*(kompt-1)
                call fcesnd(elrefe, 1, xi1, xi2, xi3, 'LI', vfesnd)
                do k = 1, 8
                    i4 = 1260+ll+k
                    zr(lzr-1+i4) = vfesnd(k)
                end do
            end do
        end do
!
        xi3 = 0.d0
        do l = 1, 8
            i1 = 810+l
            i2 = 810+9+l
            xi1 = zr(lzr-1+i1)
            xi2 = zr(lzr-1+i2)
            call fcesnd(elrefe, 0, xi1, xi2, xi3, 'LI', vfesnd)
            ll = 4*(l-1)
            do k = 1, 4
                i3 = 1452+ll+k
                zr(lzr-1+i3) = vfesnd(k)
            end do
        end do
!
!     EN NON LINEAIRE, CREATION DE LA MATRICE MAGIQUE DE PASSAGE
!     PTS DE GAUSS AUX NOEUDS PAR MOINDRES CARRES
!
        ljmas = npgcou*nso*npgcou*npgsn
        call wkvect('&INEL.'//elrefl//'.B', 'V V R', ljmas, jmas)
!
        call mamagi(elrefe, zr(lzr), zr(jmas))
!
!     DEFINITION DES PTS DE GAUSS DANS L'EPAISSEUR, DE
!     LEURS POIDS
!
        call fcepai(zr(lzr))
!
    else if (elrefe .eq. 'MEC3TR7H') then
!
!     POUR LE TRIANGLE
!
!     DEFINITION DES 3 PTS DE HAMMER REDUIT ET DES POIDS CORRESPONDANT
!
        call elraga('TR7', 'FPG3    ', ndim, nbpg, xpg, &
                    poipg)
        zr(lzr-1+1) = xpg(1)
        zr(lzr-1+2) = xpg(3)
        zr(lzr-1+3) = xpg(5)
!
        zr(lzr-1+5) = xpg(2)
        zr(lzr-1+6) = xpg(4)
        zr(lzr-1+7) = xpg(6)
!
        zr(lzr-1+9) = poipg(1)
        zr(lzr-1+10) = poipg(2)
        zr(lzr-1+11) = poipg(3)
!
!     FONCTIONS DE LAGRANGE (6 FONCTIONS)
!
!     VALEURS AUX 3 PTS DE HAMMER REDUITS
!     DES FONCTIONS DE LAGRANGE ET DE LEURS
!     DERIVEES
!
        do l = 1, 3
            i1 = l
            i2 = 4+l
            xi1 = zr(lzr-1+i1)
            xi2 = zr(lzr-1+i2)
            ll = 8*(l-1)
            ff(1) = (1.d0-xi1-xi2)*(2.d0*(1.d0-xi1-xi2)-1.d0)
            ff(2) = xi1*(2.d0*xi1-1.d0)
            ff(3) = xi2*(2.d0*xi2-1.d0)
            ff(4) = 4.d0*xi1*(1.d0-xi1-xi2)
            ff(5) = 4.d0*xi1*xi2
            ff(6) = 4.d0*xi2*(1.d0-xi1-xi2)
            dff(1, 1) = -4.d0*(1.d0-xi1-xi2)+1.d0
            dff(1, 2) = 4.d0*xi1-1.d0
            dff(1, 3) = 0.d0
            dff(1, 4) = 4.d0*(1.d0-xi1-xi2)-4.d0*xi1
            dff(1, 5) = 4.d0*xi2
            dff(1, 6) = -4.d0*xi2
            dff(2, 1) = -4.d0*(1.d0-xi1-xi2)+1.d0
            dff(2, 2) = 0.d0
            dff(2, 3) = 4.d0*xi2-1.d0
            dff(2, 4) = -4.d0*xi1
            dff(2, 5) = 4.d0*xi1
            dff(2, 6) = 4.d0*(1.d0-xi1-xi2)-4.d0*xi2
            do l1 = 1, 6
                i3 = 12+ll+l1
                i4 = 44+ll+l1
                i5 = 76+ll+l1
                zr(lzr-1+i3) = ff(l1)
                zr(lzr-1+i4) = dff(1, l1)
                zr(lzr-1+i5) = dff(2, l1)
            end do
        end do
!
!     DEFINITION DES 7 PTS DE HAMMER NORMAL ET DES POIDS CORRESPONDANTS
!
        call elraga('TR7', 'FPG7    ', ndim, nbpg, xpg, &
                    poipg)
!
        zr(lzr-1+109) = xpg(1)
        zr(lzr-1+110) = xpg(3)
        zr(lzr-1+111) = xpg(5)
        zr(lzr-1+112) = xpg(7)
        zr(lzr-1+113) = xpg(9)
        zr(lzr-1+114) = xpg(11)
        zr(lzr-1+115) = xpg(13)
!
        zr(lzr-1+118) = xpg(2)
        zr(lzr-1+119) = xpg(4)
        zr(lzr-1+120) = xpg(6)
        zr(lzr-1+121) = xpg(8)
        zr(lzr-1+122) = xpg(10)
        zr(lzr-1+123) = xpg(12)
        zr(lzr-1+124) = xpg(14)
!
        zr(lzr-1+127) = poipg(1)
        zr(lzr-1+128) = poipg(2)
        zr(lzr-1+129) = poipg(3)
        zr(lzr-1+130) = poipg(4)
        zr(lzr-1+131) = poipg(5)
        zr(lzr-1+132) = poipg(6)
        zr(lzr-1+133) = poipg(7)
!
!     VALEURS AUX 7 PTS DE HAMMER NORMAL
!     DES FONCTIONS DE LAGRANGE ET DE LEURS
!     DERVIVEEES
!
        do l = 1, 7
            i1 = 108+l
            i2 = 108+9+l
            xi1 = zr(lzr-1+i1)
            xi2 = zr(lzr-1+i2)
            ll = 8*(l-1)
            ff(1) = (1.d0-xi1-xi2)*(2.d0*(1.d0-xi1-xi2)-1.d0)
            ff(2) = xi1*(2.d0*xi1-1.d0)
            ff(3) = xi2*(2.d0*xi2-1.d0)
            ff(4) = 4.d0*xi1*(1.d0-xi1-xi2)
            ff(5) = 4.d0*xi1*xi2
            ff(6) = 4.d0*xi2*(1.d0-xi1-xi2)
            dff(1, 1) = -4.d0*(1.d0-xi1-xi2)+1.d0
            dff(1, 2) = 4.d0*xi1-1.d0
            dff(1, 3) = 0.d0
            dff(1, 4) = 4.d0*(1.d0-xi1-xi2)-4.d0*xi1
            dff(1, 5) = 4.d0*xi2
            dff(1, 6) = -4.d0*xi2
            dff(2, 1) = -4.d0*(1.d0-xi1-xi2)+1.d0
            dff(2, 2) = 0.d0
            dff(2, 3) = 4.d0*xi2-1.d0
            dff(2, 4) = -4.d0*xi1
            dff(2, 5) = 4.d0*xi1
            dff(2, 6) = 4.d0*(1.d0-xi1-xi2)-4.d0*xi2
            do l1 = 1, 6
                i3 = 135+ll+l1
                i4 = 207+ll+l1
                i5 = 279+ll+l1
                zr(lzr-1+i3) = ff(l1)
                zr(lzr-1+i4) = dff(1, l1)
                zr(lzr-1+i5) = dff(2, l1)
            end do
        end do
!
!     7 FONCTIONS CUBIQUES (6 + 1 :  LA DERNIERE EST LA 10EME
!                                    FONCTION DE P3)
!
!     VALEURS AUX 3 PTS DE HAMMER REDUITS
!
        do l = 1, 3
            i1 = l
            i2 = 4+l
            xi1 = zr(lzr-1+i1)
            xi2 = zr(lzr-1+i2)
            ll = 9*(l-1)
            ff(1) = (1.d0-xi1-xi2)*(2.d0*(1.d0-xi1-xi2)-1.d0)+3.d0*xi1* &
                    xi2*(1.d0-xi1-xi2)
            ff(2) = xi1*(2.d0*xi1-1.d0)+3.d0*xi1*xi2*(1.d0-xi1-xi2)
            ff(3) = xi2*(2.d0*xi2-1.d0)+3.d0*xi1*xi2*(1.d0-xi1-xi2)
            ff(4) = 4.d0*xi1*(1.d0-xi1-xi2)-12.d0*xi1*xi2*(1.d0-xi1- &
                                                           xi2)
            ff(5) = 4.d0*xi1*xi2-12.d0*xi1*xi2*(1.d0-xi1-xi2)
            ff(6) = 4.d0*xi2*(1.d0-xi1-xi2)-12.d0*xi1*xi2*(1.d0-xi1- &
                                                           xi2)
            ff(7) = 27.d0*xi1*xi2*(1.d0-xi1-xi2)
            dff(1, 1) = -4.d0*(1.d0-xi1-xi2)+1.d0+3.d0*xi2*(1.d0-xi1- &
                                                            xi2)-3.d0*xi1*xi2
            dff(1, 2) = 4.d0*xi1-1.d0+3.d0*xi2*(1.d0-xi1-xi2)-3.d0*xi1* &
                        xi2
            dff(1, 3) = 3.d0*xi2*(1.d0-xi1-xi2)-3.d0*xi1*xi2
            dff(1, 4) = 4.d0*(1.d0-xi1-xi2)-4.d0*xi1-12.d0*xi2*(1.d0- &
                                                                xi1-xi2)+12.d0*xi1*xi2
            dff(1, 5) = 4.d0*xi2-12.d0*xi2*(1.d0-xi1-xi2)+12.d0*xi1*xi2
            dff(1, 6) = -4.d0*xi2-12.d0*xi2*(1.d0-xi1-xi2)+12.d0*xi1*xi2
            dff(1, 7) = 27.d0*xi2*(1.d0-xi1-xi2)-27.d0*xi1*xi2
            dff(2, 1) = -4.d0*(1.d0-xi1-xi2)+1.d0+3.d0*xi1*(1.d0-xi1- &
                                                            xi2)-3.d0*xi1*xi2
            dff(2, 2) = 3.d0*xi1*(1.d0-xi1-xi2)-3.d0*xi1*xi2
            dff(2, 3) = 4.d0*xi2-1.d0+3.d0*xi1*(1.d0-xi1-xi2)-3.d0*xi1* &
                        xi2
            dff(2, 4) = -4.d0*xi1-12.d0*xi1*(1.d0-xi1-xi2)+12.d0*xi1*xi2
            dff(2, 5) = 4.d0*xi1-12.d0*xi1*(1.d0-xi1-xi2)+12.d0*xi1*xi2
            dff(2, 6) = 4.d0*(1.d0-xi1-xi2)-4.d0*xi2-12.d0*xi1*(1.d0- &
                                                                xi1-xi2)+12.d0*xi1*xi2
            dff(2, 7) = 27.d0*xi1*(1.d0-xi1-xi2)-27.d0*xi1*xi2
            do l1 = 1, 7
                i3 = 351+ll+l1
                i4 = 387+ll+l1
                i5 = 423+ll+l1
                zr(lzr-1+i3) = ff(l1)
                zr(lzr-1+i4) = dff(1, l1)
                zr(lzr-1+i5) = dff(2, l1)
            end do
        end do
!
!     VALEURS AUX 7 PTS DE HAMMER NORMAL
!
        do l = 1, 7
            i1 = 108+l
            i2 = 108+9+l
            xi1 = zr(lzr-1+i1)
            xi2 = zr(lzr-1+i2)
            ll = 9*(l-1)
            ff(1) = (1.d0-xi1-xi2)*(2.d0*(1.d0-xi1-xi2)-1.d0)+3.d0*xi1*xi2*(1.d0-xi1-xi2)
            ff(2) = xi1*(2.d0*xi1-1.d0)+3.d0*xi1*xi2*(1.d0-xi1-xi2)
            ff(3) = xi2*(2.d0*xi2-1.d0)+3.d0*xi1*xi2*(1.d0-xi1-xi2)
            ff(4) = 4.d0*xi1*(1.d0-xi1-xi2)-12.d0*xi1*xi2*(1.d0-xi1-xi2)
            ff(5) = 4.d0*xi1*xi2-12.d0*xi1*xi2*(1.d0-xi1-xi2)
            ff(6) = 4.d0*xi2*(1.d0-xi1-xi2)-12.d0*xi1*xi2*(1.d0-xi1-xi2)
            ff(7) = 27.d0*xi1*xi2*(1.d0-xi1-xi2)
            dff(1, 1) = -4.d0*(1.d0-xi1-xi2)+1.d0+3.d0*xi2*(1.d0-xi1-xi2)-3.d0*xi1*xi2
            dff(1, 2) = 4.d0*xi1-1.d0+3.d0*xi2*(1.d0-xi1-xi2)-3.d0*xi1*xi2
            dff(1, 3) = 3.d0*xi2*(1.d0-xi1-xi2)-3.d0*xi1*xi2
            dff(1, 4) = 4.d0*(1.d0-xi1-xi2)-4.d0*xi1-12.d0*xi2*(1.d0-xi1-xi2)+12.d0*xi1*xi2
            dff(1, 5) = 4.d0*xi2-12.d0*xi2*(1.d0-xi1-xi2)+12.d0*xi1*xi2
            dff(1, 6) = -4.d0*xi2-12.d0*xi2*(1.d0-xi1-xi2)+12.d0*xi1*xi2
            dff(1, 7) = 27.d0*xi2*(1.d0-xi1-xi2)-27.d0*xi1*xi2
            dff(2, 1) = -4.d0*(1.d0-xi1-xi2)+1.d0+3.d0*xi1*(1.d0-xi1-xi2)-3.d0*xi1*xi2
            dff(2, 2) = 3.d0*xi1*(1.d0-xi1-xi2)-3.d0*xi1*xi2
            dff(2, 3) = 4.d0*xi2-1.d0+3.d0*xi1*(1.d0-xi1-xi2)-3.d0*xi1*xi2
            dff(2, 4) = -4.d0*xi1-12.d0*xi1*(1.d0-xi1-xi2)+12.d0*xi1*xi2
            dff(2, 5) = 4.d0*xi1-12.d0*xi1*(1.d0-xi1-xi2)+12.d0*xi1*xi2
            dff(2, 6) = 4.d0*(1.d0-xi1-xi2)-4.d0*xi2-12.d0*xi1*(1.d0-xi1-xi2)+12.d0*xi1*xi2
            dff(2, 7) = 27.d0*xi1*(1.d0-xi1-xi2)-27.d0*xi1*xi2
            do l1 = 1, 7
                i3 = 459+ll+l1
                i4 = 540+ll+l1
                i5 = 621+ll+l1
                zr(lzr-1+i3) = ff(l1)
                zr(lzr-1+i4) = dff(1, l1)
                zr(lzr-1+i5) = dff(2, l1)
            end do
        end do
!
!    FONCTIONS ASSOCIEES AUX 3 PTS DE HAMMER REDUITS
!
!    VALEURS AUX 7 PTS DE HAMMER NORMAL
!
        do l = 1, 7
            i1 = 108+l
            i2 = 108+9+l
            x(1) = zr(lzr-1+i1)
            x(2) = zr(lzr-1+i2)
            ll = 4*(l-1)
            a = 0.166666666666667d0
            b = 0.666666666666667d0
            ff(1) = 2.d0*(-(x(1)-b)-(x(2)-a))
            ff(2) = 2.d0*(x(1)-a)
            ff(3) = 2.d0*(x(2)-a)
            dff(1, 1) = -2.d0
            dff(1, 2) = 2.d0
            dff(1, 3) = 0.d0
            dff(2, 1) = -2.d0
            dff(2, 2) = 0.d0
            dff(2, 3) = 2.d0
            do l1 = 1, 3
                i3 = 702+ll+l1
                i4 = 738+ll+l1
                i5 = 774+ll+l1
                zr(lzr-1+i3) = ff(l1)
                zr(lzr-1+i4) = dff(1, l1)
                zr(lzr-1+i5) = dff(2, l1)
            end do
        end do
!
!     VALEURS DES FONCTIONS DE LAGRANGE AUX 7 NOEUDS
!     POSITION DES 7 NOEUDS
!
        zr(lzr-1+811) = 0.000000000000000d0
        zr(lzr-1+812) = 1.000000000000000d0
        zr(lzr-1+813) = 0.000000000000000d0
        zr(lzr-1+814) = 0.500000000000000d0
        zr(lzr-1+815) = 0.500000000000000d0
        zr(lzr-1+816) = 0.000000000000000d0
        zr(lzr-1+817) = 0.333333333333333d0
!
        zr(lzr-1+820) = 0.000000000000000d0
        zr(lzr-1+821) = 0.000000000000000d0
        zr(lzr-1+822) = 1.000000000000000d0
        zr(lzr-1+823) = 0.000000000000000d0
        zr(lzr-1+824) = 0.500000000000000d0
        zr(lzr-1+825) = 0.500000000000000d0
        zr(lzr-1+826) = 0.333333333333333d0
!
        do l = 1, 7
            i1 = 810+l
            i2 = 810+9+l
            xi1 = zr(lzr-1+i1)
            xi2 = zr(lzr-1+i2)
            ll = 8*(l-1)
            dff(1, 1) = -4.d0*(1.d0-xi1-xi2)+1.d0
            dff(1, 2) = 4.d0*xi1-1.d0
            dff(1, 3) = 0.d0
            dff(1, 4) = 4.d0*(1.d0-xi1-xi2)-4.d0*xi1
            dff(1, 5) = 4.d0*xi2
            dff(1, 6) = -4.d0*xi2
            dff(2, 1) = -4.d0*(1.d0-xi1-xi2)+1.d0
            dff(2, 2) = 0.d0
            dff(2, 3) = 4.d0*xi2-1.d0
            dff(2, 4) = -4.d0*xi1
            dff(2, 5) = 4.d0*xi1
            dff(2, 6) = 4.d0*(1.d0-xi1-xi2)-4.d0*xi2
            do l1 = 1, 6
                i3 = 828+ll+l1
                i4 = 900+ll+l1
                zr(lzr-1+i3) = dff(1, l1)
                zr(lzr-1+i4) = dff(2, l1)
            end do
        end do
!
!     DEFINITION DES 6 FONCTIONS D'INTERPOLATION QUI PERMETTENT
!     D'EXTRAPOLER
!     LES DEFORMATIONS OU CONTRAINTES AUX NOEUDS A PARTIR DE
!     LEURS VALEURS AUX POINTS DE HAMMER
!
        zr(lzr-1+1251) = -1.d0
        zr(lzr-1+1252) = 0.d0
        zr(lzr-1+1253) = 1.d0
!
        kompt = 0
        do m = 1, 3
            i3 = 1250+m
            xi3 = zr(lzr-1+i3)
            do l = 1, 6
                i1 = 810+l
                i2 = 810+9+l
                xi1 = zr(lzr-1+i1)
                xi2 = zr(lzr-1+i2)
                kompt = kompt+1
                ll = 6*(kompt-1)
                call fcesnd(elrefe, 1, xi1, xi2, xi3, 'LI', vfesnd)
                do k = 1, 6
                    i4 = 1260+ll+k
                    zr(lzr-1+i4) = vfesnd(k)
                end do
            end do
        end do
!
        xi3 = 0.d0
        do l = 1, 6
            i1 = 810+l
            i2 = 810+9+l
            xi1 = zr(lzr-1+i1)
            xi2 = zr(lzr-1+i2)
            call fcesnd(elrefe, 0, xi1, xi2, xi3, 'LI', vfesnd)
            ll = 4*(l-1)
            do k = 1, 3
                i3 = 1452+ll+k
                zr(lzr-1+i3) = vfesnd(k)
            end do
        end do
!
!     EN NON LINEAIRE, CREATION DE LA MATRICE MAGIQUE DE PASSAGE
!     PTS DE HAMMER AUX NOEUDS PAR MOINDRES CARRES
!
        ljmas = npgcou*nso*npgcou*npgsn
        call wkvect('&INEL.'//elrefl//'.B', 'V V R', ljmas, jmas)
!
        call mamagi(elrefe, zr(lzr), zr(jmas))
!
!     DEFINITION DES PTS DE GAUSS DANS L'EPAISSEUR, DE
!     LEURS POIDS
!
        call fcepai(zr(lzr))
!
!
!
!     DEFINITION DES 7 FONCTIONS D'INTERPOLATION QUI PERMETTENT
!     D'EXTRAPOLER
!     LES DEFORMATIONS GENERALISEES AUX NOEUDS A PARTIR DE
!     LEURS VALEURS AUX POINTS DE HAMMER
!
        xi3 = 0.d0
        do l = 1, 7
            i1 = 810+l
            i2 = 810+9+l
            xi1 = zr(lzr-1+i1)
            xi2 = zr(lzr-1+i2)
            call fcesnd(elrefe, 0, xi1, xi2, xi3, 'NL', vfesnd)
            ll = 7*(l-1)
            do k = 1, 7
                i3 = 1600+ll+k
                zr(lzr-1+i3) = vfesnd(k)
            end do
        end do

    end if
!
999 continue
end subroutine
