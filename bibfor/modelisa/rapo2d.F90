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
! aslint: disable=W0413
! aslint: disable=I0413
!
subroutine rapo2d(numdlz, iocc, fonrez, lisrez, chargz)
    implicit none
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterfort/afrela.h"
#include "asterfort/assert.h"
#include "asterfort/assvec.h"
#include "asterfort/calcul.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exlim1.h"
#include "asterfort/getvtx.h"
#include "asterfort/imprel.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/malin1.h"
#include "asterfort/mecact.h"
#include "asterfort/memare.h"
#include "asterfort/mesomm.h"
#include "asterfort/reajre.h"
#include "asterfort/reliem.h"
#include "asterfort/utmess.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "asterfort/char8_to_int.h"
!
    character(len=*) :: numdlz, chargz, fonrez, lisrez
! -------------------------------------------------------
!     RACCORD POUTRE-2D PAR DES RELATIONS LINEAIRES
!     ENTRE LE NOEUD DES MAILLES DE BORD DE LA STRUCTURE 2D
!     ET LE NOEUD DE LA POUTRE DONNE PAR L'UTILISATEUR
! -------------------------------------------------------
!  NUMDLZ        - IN    - K14  - : NOM DU NUMEDDL DU LIGREL DU MODELE
!                                     (IL A ETE CREE SUR LA VOLATILE)
!  IOCC          - IN    - I    - : NUMERO D'OCCURENCE DU MOT-FACTEUR
!  FONREZ        - IN    - K4   - : 'REEL'
!  LISREZ        - IN    - K19  - : NOM DE LA SD LISTE_RELA
!  CHARGE        - IN    - K8   - : NOM DE LA SD CHARGE
!                - JXVAR -      -
! -------------------------------------------------------
!
!
! --------- VARIABLES LOCALES ---------------------------
    integer(kind=8) :: nmocl
    parameter(nmocl=300)
    character(len=4) :: typval, typcoe
    character(len=8) :: betaf, mod, nomg
    character(len=8) :: noma, nomcmp(nmocl)
    character(len=8) :: noepou, nocmp(3), kcmp(3)
    character(len=8) :: lpain(2), lpaout(2)
    character(len=9) :: nomte
    character(len=16) :: motfac, motcle(4), typmcl(4), option
    character(len=19) :: ligrmo, ligrel
    character(len=24) :: lchin(2), lchout(2), nolili, lismai, valk(2)
    character(len=24) :: lisnoe, vale1, vale2, grnoma
    character(len=8) :: charge
    character(len=14) :: numddl
    character(len=19) :: lisrel
    integer(kind=8) :: ntypel(nmocl), icmp(6), niv, ifm, vali(2)
    integer(kind=8) :: iop, nliai, i, inom
    integer(kind=8) :: nbcmp, nddla, nbec, jprnm, nlili, k, iaprno, lonlis, ilisno
    integer(kind=8) :: jlisma, nbma, nbno, numnop
    integer(kind=8) :: ino, idch1, idch2, nbterm, jno2
    integer(kind=8) ::       ival
    integer(kind=8) :: iocc
    real(kind=8) :: igzz, coorig(3), beta, eps, un
    real(kind=8) :: xpou, ypou, s, s1, xg, yg, dnorme
    real(kind=8) :: ax, ay, axx, ayy, valr(9)
    complex(kind=8) :: betac, ccmp(3)
    complex(kind=8), pointer :: coec(:) => null()
    real(kind=8), pointer :: coer(:) => null()
    integer(kind=8), pointer :: dime(:) => null()
    real(kind=8), pointer :: direct(:) => null()
    real(kind=8), pointer :: inertie_raccord(:) => null()
    character(len=8), pointer :: lisddl(:) => null()
    character(len=8), pointer :: lisno(:) => null()
    character(len=8), pointer :: lgrf(:) => null()
    real(kind=8), pointer :: vale(:) => null()
! --------- FIN  DECLARATIONS  VARIABLES LOCALES --------
!
    call jemarq()
! --- RECUPERATION DES PARAMETRE D IMPRESSION
    call infniv(ifm, niv)
! -------------------------------------------------------
    numddl = numdlz
    charge = chargz
    lisrel = lisrez
!
    motfac = 'LIAISON_ELEM'
    call getvtx(motfac, 'OPTION', iocc=iocc, scal=option, nbret=iop)
    if (option .ne. '2D_POU') then
        call utmess('F', 'MODELISA6_39', sk=option)
    end if
!
    call getfac(motfac, nliai)
    if (nliai .eq. 0) goto 999
!
! --- INITIALISATIONS
!     ---------------
!
! --- TYPE DES VALEURS AU SECOND MEMBRE DES RELATIONS
    typval = fonrez
! --- TYPE DES VALEURS DES COEFFICIENTS DES RELATIONS
    typcoe = 'REEL'
! --- VALEUR DU SECOND MEMBRE DES RELATIONS QUAND C'EST UNE FONCTION
    betaf = '&FOZERO'
! --- VALEUR DU SECOND MEMBRE DES RELATIONS QUAND C'EST UN REEL
    beta = 0.0d0
! --- VALEUR DU SECOND MEMBRE DES RELATIONS QUAND C'EST UN COMPLEXE
    betac = (0.0d0, 0.0d0)
    eps = 1.0d-02
    un = 1.0d0
    kcmp(1) = ' '
    kcmp(2) = ' '
    kcmp(3) = ' '
    ccmp(1) = (0.0d0, 0.0d0)
    ccmp(2) = (0.0d0, 0.0d0)
    ccmp(3) = (0.0d0, 0.0d0)
    do i = 1, 6
        icmp(i) = 0
    end do
!
    ligrel = '&&RAPO2D'
    lisnoe = '&&RAPO2D.LISTE_NOEUDS'
    lismai = '&&RAPO2D.LISTE_MAILLES'
    motcle(1) = 'GROUP_MA_1'
    motcle(2) = 'MAILLE_1'
    typmcl(1) = 'GROUP_MA'
    typmcl(2) = 'MAILLE'
!
! --- -----------------------------------------------------------------
! --- MODELE ASSOCIE AU LIGREL DE CHARGE
    call dismoi('NOM_MODELE', charge(1:8), 'CHARGE', repk=mod)
!
! --- -----------------------------------------------------------------
!     LIGREL DU MODELE
    ligrmo = mod(1:8)//'.MODELE'
!
! --- -----------------------------------------------------------------
! --- MAILLAGE ASSOCIE AU MODELE
    call jeveuo(ligrmo//'.LGRF', 'L', vk8=lgrf)
    noma = lgrf(1)
    grnoma = noma//'.GROUPENO'
!
! --- -----------------------------------------------------------------
! --- RECUPERATION DU TABLEAU DES COORDONNEES
    call jeveuo(noma//'.COORDO    .VALE', 'L', vr=vale)
!
! --- -----------------------------------------------------------------
! --- RECUPERATION DES NOMS DES DDLS
    nomg = 'DEPL_R'
    nomte = 'D_DEPL_R_'
!
    call jeveuo(jexnom('&CATA.GD.NOMCMP', nomg), 'L', inom)
    call jelira(jexnom('&CATA.GD.NOMCMP', nomg), 'LONMAX', nbcmp)
    nddla = nbcmp-1
    if (nddla .gt. nmocl) then
        vali(1) = nmocl
        vali(2) = nddla
        call utmess('F', 'MODELISA8_29', ni=2, vali=vali)
    end if
    do i = 1, nddla
        nomcmp(i) = zk8(inom-1+i)
        call jenonu(jexnom('&CATA.TE.NOMTE', nomte//nomcmp(i) (1:7)), ntypel(i))
    end do
    call dismoi('NB_EC', nomg, 'GRANDEUR', repi=nbec)
!
! --- -----------------------------------------------------------------
! --- ACCES A L'OBJET .PRNM
    if (nbec .gt. 10) then
        call utmess('F', 'MODELISA_94')
    else
        call jeveuo(ligrmo//'.PRNM', 'L', jprnm)
    end if
!
! --- -----------------------------------------------------------------
! --- RECUPERATION DU .PRNO ASSOCIE AU MAILLAGE
    call jelira(numddl//'.NUME.PRNO', 'NMAXOC', nlili)
    k = 0
    do i = 1, nlili
        call jenuno(jexnum(numddl//'.NUME.LILI', i), nolili)
        if (nolili(1:8) .ne. '&MAILLA') goto 30
        k = i
30      continue
    end do
    ASSERT(k .ne. 0)
    call jeveuo(jexnum(numddl//'.NUME.PRNO', k), 'L', iaprno)
!
! --- -----------------------------------------------------------------
! --- ACQUISITION DE LA LISTE DES NOEUDS A LIER
!     (CETTE LISTE EST NON REDONDANTE)
    call malin1(motfac, charge, iocc, 1, lisnoe, &
                lonlis)
    call jeveuo(lisnoe, 'L', ilisno)
!
! --- -----------------------------------------------------------------
! --- CONSTITUTION DU LIGREL FORME DES MAILLES DE BORD DE LA SURFACE 2D
!
! --- CREATION ET AFFECTATION DU VECTEUR DE K8 DE NOM LISMAI
!     CONTENANT LES NOMS DES MAILLES FORMANT LE LIGREL A CREER
    call reliem(' ', noma, 'NU_MAILLE', motfac, iocc, &
                2, motcle(1), typmcl(1), lismai, nbma)
    call jeveuo(lismai, 'L', jlisma)
!
!     CREATION ET AFFECTATION DU LIGREL
    call exlim1(zi(jlisma), nbma, mod, 'V', ligrel)
!
! --- -----------------------------------------------------------------
! --- Recuperation du noeud "poutre" (_2) :
    motcle(1) = 'GROUP_NO_2'
    motcle(2) = 'NOEUD_2'
    motcle(3) = 'GROUP_MA_2'
    motcle(4) = 'MAILLE_2'
    typmcl(1) = 'GROUP_NO'
    typmcl(2) = 'NOEUD'
    typmcl(3) = 'GROUP_MA'
    typmcl(4) = 'MAILLE'
    call reliem(' ', noma, 'NO_NOEUD', motfac, iocc, &
                4, motcle, typmcl, '&&RAPO2D.NO2', nbno)

    if (nbno .ne. 1) then
        call utmess('F', 'MODELISA6_40', si=nbno)
    end if
    call jeveuo('&&RAPO2D.NO2', 'L', jno2)
    noepou = zk8(jno2)
    numnop = char8_to_int(noepou)
    call jedetr('&&RAPO2D.NO2')
    numnop = char8_to_int(noepou)
! --- COORDONNEES DU NOEUD POUTRE
    xpou = vale(3*(numnop-1)+1)
    ypou = vale(3*(numnop-1)+2)
!
! --- -----------------------------------------------------------------
! --- CALCUL SUR CHAQUE ELEMENT DE BORD A RELIER A LA POUTRE
!     DES CARACTERISTIQUES GEOMETRIQUES SUIVANTES :
!        SOMME/B_ELEMENT(1,X,Y,X*X,Y*Y,X*Y)DS
    lpain(1) = 'PGEOMER'
    lchin(1) = noma//'.COORDO'
    lpaout(1) = 'PCASECT'
    lchout(1) = '&&RAPO2D.PSECT'
!
    call calcul('S', 'CARA_SECT_POUT3', ligrel, 1, lchin, &
                lpain, 1, lchout, lpaout, 'V', &
                'OUI')
!
! --- -----------------------------------------------------------------
! --- VECTEUR DES QUANTITES GEOMETRIQUES PRECITEES SOMMEES
!     SUR LA SURFACE DE RACCORD, CES QUANTITES SERONT NOTEES :
!        A1 = S,AX,AY,AXX,AYY
    AS_ALLOCATE(vr=inertie_raccord, size=6)
! --- -----------------------------------------------------------------
!     SOMMATION DES QUANTITES GEOMETRIQUES ELEMENTAIRES
!     DANS LE VECTEUR &&RAPO2D.INERTIE_RACCORD :
    call mesomm(lchout(1), 6, vr=inertie_raccord)
!
    s = inertie_raccord(1)
    ax = inertie_raccord(2)
    ay = inertie_raccord(3)
    axx = inertie_raccord(4)
    ayy = inertie_raccord(5)
!
    if (abs(s) .eq. 0.d0) then
        call utmess('F', 'MODELISA6_46')
    end if
    s1 = 1.0d0/s
!
! --- -----------------------------------------------------------------
! --- COORDONNEES DU CENTRE GEOMETRIQUE G DE LA SECTION DE RACCORD
!     XG = AX/S, YG = AY/S
    xg = s1*ax
    yg = s1*ay
!
! --- -----------------------------------------------------------------
!     VERIFICATION DE L'IDENTITE GEOMETRIQUE DE G AVEC LE
!     NOEUD POUTRE A RACCORDER :
    dnorme = sqrt((xpou-xg)*(xpou-xg)+(ypou-yg)*(ypou-yg))/(s/2.)
    if (dnorme .gt. eps) then
        valr(1) = xg
        valr(2) = yg
        valr(3) = 0.
        valr(4) = xpou
        valr(5) = ypou
        valr(6) = 0.
        valr(7) = eps*100.0d0
        valr(8) = s
        valr(9) = dnorme
        valk(1) = option
        vali(1) = iocc
        call utmess('A', 'CALCULEL3_80', sk=valk(1), si=vali(1), nr=9, &
                    valr=valr)
    end if
!
!
! --- -----------------------------------------------------------------
!     CALCUL DE LA COMPOSANTE IZZ DU TENSEUR D'INERTIE EN G
    igzz = axx+ayy-s*(xg*xg+yg*yg)
!
! --- -----------------------------------------------------------------
!     NOTATION DANS LA CARTE DE NOM '&&RAPO2D.CAORIGE' DES
!     COORDONNEES DU CENTRE GEOMETRIQUE G DE LA SECTION DE RACCORD
    nocmp(1) = 'X'
    nocmp(2) = 'Y'
!
    coorig(1) = xg
    coorig(2) = yg
!
    call mecact('V', '&&RAPO2D.CAORIGE', 'LIGREL', ligrel, 'GEOM_R', &
                ncmp=2, lnomcmp=nocmp, vr=coorig)
!
! --- DETERMINATION DE 2 LISTES  DE VECTEURS PAR ELEMENT PRENANT
!     LEURS VALEURS AUX NOEUDS DES ELEMENTS.
! --- LA PREMIERE LISTE DE NOM 'VECT_NI' A POUR VALEURS AU NOEUD
!     I D'UN ELEMENT : SOMME/S_ELEMENT(NI,0)DS
! --- LA SECONDE LISTE DE NOM 'VECT_XYNI' A POUR VALEURS AU NOEUD
!     I D'UN ELEMENT : SOMME/S_ELEMENT(X*NI,Y*NI)DS
!        AVEC  X = XM - XG = NJ*XJ - XG
!              Y = YM - YG = NJ*YJ - YG
!
    lpain(1) = 'PGEOMER'
    lpain(2) = 'PORIGIN'
    lchin(1) = noma//'.COORDO'
    lchin(2) = '&&RAPO2D.CAORIGE'
    lpaout(1) = 'PVECTU1'
    lpaout(2) = 'PVECTU2'
    lchout(1) = '&&RAPO2D.VECT_NI'
    lchout(2) = '&&RAPO2D.VECT_XYZNI'
!
    call calcul('S', 'CARA_SECT_POUT4', ligrel, 2, lchin, &
                lpain, 2, lchout, lpaout, 'V', &
                'OUI')
!
! --- -----------------------------------------------------------------
! --- CREATION DES .RERR DES VECTEURS EN SORTIE DE CALCUL
    call memare('V', '&&RAPO2D', mod, 'CHAR_MECA')
!
! --- -----------------------------------------------------------------
!     ASSEMBLAGE DE LCHOUT(1) DANS LE CHAMNO DE NOM 'CH_DEPL_1'
    call jedetr('&&RAPO2D           .RELR')
    call reajre('&&RAPO2D', lchout(1), 'V')
    call assvec('V', 'CH_DEPL_1', 1, '&&RAPO2D           .RELR', [1.d0], numddl)
!
! --- -----------------------------------------------------------------
!     ASSEMBLAGE DE LCHOUT(2) DANS LE CHAMNO DE NOM 'CH_DEPL_2'
    call jedetr('&&RAPO2D           .RELR')
    call reajre('&&RAPO2D', lchout(2), 'V')
    call assvec('V', 'CH_DEPL_2', 1, '&&RAPO2D           .RELR', [1.d0], numddl)
!
    vale1 = 'CH_DEPL_1          .VALE'
    vale2 = 'CH_DEPL_2          .VALE'
    call jeveuo(vale1, 'L', idch1)
    call jeveuo(vale2, 'L', idch2)
!
! --- -----------------------------------------------------------------
! --- CREATION DES TABLEAUX NECESSAIRES A L'AFFECTATION DE LISREL
! --- MAJORANT DU NOMBRE DE TERMES DANS UNE RELATION
    nbterm = 2*lonlis+2
! --- VECTEUR DU NOM DES NOEUDS
    AS_ALLOCATE(vk8=lisno, size=nbterm)
! --- VECTEUR DU NOM DES DDLS
    AS_ALLOCATE(vk8=lisddl, size=nbterm)
! --- VECTEUR DES COEFFICIENTS REELS
    AS_ALLOCATE(vr=coer, size=nbterm)
! --- VECTEUR DES COEFFICIENTS COMPLEXES
    AS_ALLOCATE(vc=coec, size=nbterm)
! --- VECTEUR DES DIRECTIONS DES DDLS A CONTRAINDRE
    AS_ALLOCATE(vr=direct, size=2*nbterm)
! --- VECTEUR DES DIMENSIONS DE CES DIRECTIONS
    AS_ALLOCATE(vi=dime, size=nbterm)
!
! --- -----------------------------------------------------------------
! --- RELATIONS ENTRE LES NOEUDS DU BORD ET LE NOEUD POUTRE
!
! --- PREMIER GROUPE DE RELATIONS TRADUISANT :
!        SOMME/S_RACCORD(U_3D) = S_RACCORD*U_NOEUD_POUTRE
!
! --- -----------------------------------------------------------------
! --- PREMIERE RELATION :
!     -S.DX(NOEUD_POUTRE) + (SOMME/S_RACCORD(NI.DS)).DX(NOEUD_I) = 0
    nbterm = lonlis+1
!     BOUCLE SUR LES NOEUDS DES MAILLES DE LA TRACE DE LA POUTRE
    do i = 1, lonlis
        ino = char8_to_int(zk8(ilisno+i-1))
!        ADRESSE DE LA PREMIERE COMPOSANTE DU NOEUD INO DANS LES CHAMNO
        ival = zi(iaprno+(ino-1)*(nbec+2))
!
        lisno(i) = zk8(ilisno+i-1)
        lisddl(i) = 'DX'
        coer(i) = zr(idch1+ival-1)
    end do
!
    lisno(1+lonlis+1-1) = noepou
    lisddl(1+lonlis+1-1) = 'DX'
    coer(1+lonlis+1-1) = -s
!
    call afrela(coer, coec, lisddl, lisno, dime, &
                direct, nbterm, beta, betac, betaf, &
                typcoe, typval, 0.d0, lisrel)
    call imprel(motfac, nbterm, coer, lisddl, lisno, &
                beta)
!
! --- DEUXIEME RELATION :
!     -S.DY(NOEUD_POUTRE) + (SOMME/S_RACCORD(NI.DS)).DY(NOEUD_I) = 0
    nbterm = lonlis+1
!     BOUCLE SUR LES NOEUDS DES MAILLES DE LA TRACE DE LA POUTRE
    do i = 1, lonlis
        ino = char8_to_int(zk8(ilisno+i-1))
!        ADRESSE DE LA PREMIERE COMPOSANTE DU NOEUD INO DANS LES CHAMNO
        ival = zi(iaprno+(ino-1)*(nbec+2))
!
        lisno(i) = zk8(ilisno+i-1)
        lisddl(i) = 'DY'
        coer(i) = zr(idch1+ival-1)
    end do
!
    lisno(1+lonlis+1-1) = noepou
    lisddl(1+lonlis+1-1) = 'DY'
    coer(1+lonlis+1-1) = -s
!
    call afrela(coer, coec, lisddl, lisno, dime, &
                direct, nbterm, beta, betac, betaf, &
                typcoe, typval, 0.d0, lisrel)
    call imprel(motfac, nbterm, coer, lisddl, lisno, &
                beta)
!
!
! --- -----------------------------------------------------------------
! --- DEUXIEME GROUPE DE RELATIONS TRADUISANT :
!        SOMME/S_RACCORD(GM X U_3D) = I.OMEGA(NOEUD_POUTRE)
!
!
! --- TOISIEME RELATION :
!        (SOMME/S_RACCORD(X*NI.DS)).DY(NOEUD_I) -
!        (SOMME/S_RACCORD(Y*NI.DS)).DX(NOEUD_I) -
!        IZZ.DRZ(NOEUD_POUTRE)                            = 0
    nbterm = 2*lonlis+1
!     BOUCLE SUR LES NOEUDS DES MAILLES DE SURFACE DU MASSIF
    do i = 1, lonlis
        ino = char8_to_int(zk8(ilisno+i-1))
!        ADRESSE DE LA PREMIERE COMPOSANTE DU NOEUD INO DANS LES CHAMNO
        ival = zi(iaprno+(ino-1)*(nbec+2))
!
        lisno(1+2*(i-1)) = zk8(ilisno+i-1)
        lisno(1+2*(i-1)+1) = zk8(ilisno+i-1)
        lisddl(1+2*(i-1)) = 'DY'
        lisddl(1+2*(i-1)+1) = 'DX'
        coer(1+2*(i-1)) = zr(idch2+ival-1)
        coer(1+2*(i-1)+1) = -zr(idch2+ival-1+1)
    end do
!
    lisno(1+2*lonlis+1-1) = noepou
    lisddl(1+2*lonlis+1-1) = 'DRZ'
    coer(1+2*lonlis+1-1) = -igzz
!
    call afrela(coer, coec, lisddl, lisno, dime, &
                direct, nbterm, beta, betac, betaf, &
                typcoe, typval, 0.d0, lisrel)
    call imprel(motfac, nbterm, coer, lisddl, lisno, &
                beta)
!
!
! --- -----------------------------------------------------------------
! --- DESTRUCTION DES OBJETS DE TRAVAIL
!
    call jedetr('&&RAPO2D.LISTE_NOEUDS')
    call jedetr('&&RAPO2D.LISTE_MAILLES')
    call detrsd('CHAMP_GD', '&&RAPO2D.PSECT')
    AS_DEALLOCATE(vr=inertie_raccord)
    call detrsd('CARTE', '&&RAPO2D.CAORIGE')
    call detrsd('RESUELEM', '&&RAPO2D.VECT_NI')
    call detrsd('RESUELEM', '&&RAPO2D.VECT_XYZNI')
    call jedetr('&&RAPO2D           .RELR')
    call jedetr('&&RAPO2D           .RERR')
    AS_DEALLOCATE(vk8=lisno)
    AS_DEALLOCATE(vk8=lisddl)
    AS_DEALLOCATE(vr=coer)
    AS_DEALLOCATE(vc=coec)
    AS_DEALLOCATE(vr=direct)
    AS_DEALLOCATE(vi=dime)
    call detrsd('LIGREL', ligrel)
    call detrsd('CHAMP_GD', 'CH_DEPL_1')
    call detrsd('CHAMP_GD', 'CH_DEPL_2')
!
999 continue
    call jedema()
end subroutine
