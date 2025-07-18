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
!
subroutine rapoco(numdlz, iocc, fonrez, lisrez, chargz)
!
    implicit none
!
!    ATTENTION CETTE PROGRAMMATION SUPPOSE QUE L'OBJET NUEQ EST UN
!    VECTEUR IDENTITE. A MODIFIER
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterc/indik8.h"
#include "asterc/r8prem.h"
#include "asterfort/afrela.h"
#include "asterfort/assert.h"
#include "asterfort/assvec.h"
#include "asterfort/calcul.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisdg.h"
#include "asterfort/exlim1.h"
#include "asterfort/getvem.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
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
#include "asterfort/racotu.h"
#include "asterfort/reajre.h"
#include "asterfort/reliem.h"
#include "asterfort/utmess.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "asterfort/char8_to_int.h"
#include "asterfort/int_to_char8.h"
#include "asterc/r8pi.h"
!
    integer(kind=8) :: iocc
    character(len=8) :: charge
    character(len=14) :: numddl
    character(len=19) :: lisrel
    character(len=*) :: numdlz, chargz, fonrez, lisrez
! -------------------------------------------------------
!     RACCORD POUTRE-COQUE PAR DES RELATIONS LINEAIRES
!     ENTRE LES NOEUDS DES MAILLES DE BORD MODELISANT
!     LA TRACE DE LA SECTION DE LA POUTRE SUR LA COQUE
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
    character(len=8) :: betaf, mod, nomg, k8bid, cara
    character(len=8) :: noma, nomcmp(nmocl), option
    character(len=8) :: noepou, nocmp(3), kcmp(3), cmp(6)
    character(len=8) :: lpain(4), lpaout(2)
    character(len=9) :: nomte
    character(len=16) :: motfac, motcle(2), typmcl(2)
    character(len=19) :: ligrmo, ligrel
    character(len=24) :: lchin(4), lchout(2), nolili, lismai, valk(2)
    character(len=24) :: lisnoe, vale1, grnoma, vale2, nogrno
    integer(kind=8) :: ntypel(nmocl), dg, icmp(6), niv, ifm, iop, numnop, nliai
    integer(kind=8) :: vali(2), nlili, nbterm, ncara, nddla, nbma, nbno, nno, nbec
    integer(kind=8) :: nbcmp
    integer(kind=8) :: naxe, lonlis, k, j, in, ino, i, ival, n1
    integer(kind=8) :: nbgno
    integer(kind=8) ::   jlisma, jgro
    integer(kind=8) ::  iaprno, idch2, idch1, ilisno, inom
    real(kind=8) :: ig(6), coorig(3), axepou(3), valr(9)
    real(kind=8) :: ayz, axx, ax, ay, axz, axy, ayy, azz, az, beta, dnorme, eps
    real(kind=8) :: un, pi
    real(kind=8) :: xg, yg, zg, xpou, ypou, zpou, xnorm, s1, s
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
    integer(kind=8), pointer :: prnm(:) => null()
! --------- FIN  DECLARATIONS  VARIABLES LOCALES --------
!
    call jemarq()
! --- RECUPERATION DES PARAMETRE D IMPRESSION
    call infniv(ifm, niv)
! -------------------------------------------------------
    numddl = numdlz
    charge = chargz
    lisrel = lisrez
    pi = r8pi()
!
    motfac = 'LIAISON_ELEM'
    call getvtx(motfac, 'OPTION', iocc=iocc, scal=option, nbret=iop)
    if ((option .ne. 'COQ_POU') .and. (option .ne. 'COQ_TUYA')) then
        call utmess('F', 'MODELISA6_39', sk=option)
    end if
!
    call getfac(motfac, nliai)
    if (nliai .eq. 0) goto 130
!
! --- INITIALISATIONS :
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
    eps = 1.d-2
    un = 1.0d0
    kcmp(1) = ' '
    kcmp(2) = ' '
    kcmp(3) = ' '
    cmp(1) = 'DX'
    cmp(2) = 'DY'
    cmp(3) = 'DZ'
    cmp(4) = 'DRX'
    cmp(5) = 'DRY'
    cmp(6) = 'DRZ'
    ccmp(1) = (0.0d0, 0.0d0)
    ccmp(2) = (0.0d0, 0.0d0)
    ccmp(3) = (0.0d0, 0.0d0)
    do i = 1, 6
        icmp(i) = 0
    end do
!
    ligrel = '&&RAPOCO'
    lisnoe = '&&RAPOCO.LISTE_NOEUDS'
    lismai = '&&RAPOCO.LISTE_MAILLES'
    motcle(1) = 'GROUP_MA_1'
    motcle(2) = 'MAILLE_1'
    typmcl(1) = 'GROUP_MA'
    typmcl(2) = 'MAILLE'
!
! --- MODELE ASSOCIE AU LIGREL DE CHARGE :
!     ----------------------------------
    call dismoi('NOM_MODELE', charge(1:8), 'CHARGE', repk=mod)
!
! ---  LIGREL DU MODELE :
!      ----------------
    ligrmo = mod(1:8)//'.MODELE'
!
! --- MAILLAGE ASSOCIE AU MODELE :
!     --------------------------
    call jeveuo(ligrmo//'.LGRF', 'L', vk8=lgrf)
    noma = lgrf(1)
!
    grnoma = noma//'.GROUPENO'
!
! --- RECUPERATION DU TABLEAU DES COORDONNEES :
!     ---------------------------------------
    call jeveuo(noma//'.COORDO    .VALE', 'L', vr=vale)
!
! --- RECUPERATION DES NOMS DES DDLS :
!     ------------------------------
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
! --- ACCES A L'OBJET .PRNM :
!     ----------------------
    if (nbec .gt. 10) then
        call utmess('F', 'MODELISA_94')
    else
        call jeveuo(ligrmo//'.PRNM', 'L', vi=prnm)
    end if
!
! --- RECUPERATION DU .PRNO ASSOCIE AU MAILLAGE :
!     -----------------------------------------
    call jelira(numddl//'.NUME.PRNO', 'NMAXOC', nlili)
!
    k = 0
    do i = 1, nlili
        call jenuno(jexnum(numddl//'.NUME.LILI', i), nolili)
        if (nolili(1:8) .ne. '&MAILLA') goto 30
        k = i
30      continue
    end do
!
    ASSERT(k .ne. 0)
!
    call jeveuo(jexnum(numddl//'.NUME.PRNO', k), 'L', iaprno)
!
! --- ACQUISITION DE LA LISTE DES NOEUDS A LIER
! --- (CETTE LISTE EST NON REDONDANTE) :
!     -------------------------------
    call malin1(motfac, charge, iocc, 1, lisnoe, &
                lonlis)
    call jeveuo(lisnoe, 'L', ilisno)
!
! --- CONSTITUTION DU LIGREL FORME DES MAILLES DE BORD MODELISANT
! --- LA TRACE DE LA SECTION DE LA POUTRE SUR LA COQUE :
!     ------------------------------------------------
!
! --- CREATION ET AFFECTATION DU VECTEUR DE K8 DE NOM LISMAI
! --- CONTENANT LES NOMS DES MAILLES FORMANT LE LIGREL A CREER :
!     --------------------------------------------------------
    call reliem(' ', noma, 'NU_MAILLE', motfac, iocc, &
                2, motcle(1), typmcl(1), lismai, nbma)
    call jeveuo(lismai, 'L', jlisma)
!
! ---   CREATION ET AFFECTATION DU LIGREL
    call exlim1(zi(jlisma), nbma, mod, 'V', ligrel)
!
! --- ACQUISITION DES MOTS-CLES NOEUD_2 OU GROUP_NO_2 :
!     -----------------------------------------------
    nbno = 0
    nbgno = 0
!
    call getvem(noma, 'NOEUD', motfac, 'NOEUD_2', iocc, &
                0, k8bid, nbno)
!
    if (nbno .eq. 0) then
        call getvem(noma, 'GROUP_NO', motfac, 'GROUP_NO_2', iocc, &
                    0, k8bid, nbgno)
        if (nbgno .eq. 0) then
            valk(1) = motfac
            valk(2) = option
            call utmess('F', 'MODELISA6_48', nk=2, valk=valk)
        end if
    end if
!
    if (nbno .ne. 0) then
        nbno = -nbno
        if (nbno .ne. 1) then
            call utmess('F', 'MODELISA6_49')
        end if
        call getvem(noma, 'NOEUD', motfac, 'NOEUD_2', iocc, &
                    nbno, noepou, nno)
    end if
!
    if (nbgno .ne. 0) then
        nbgno = -nbgno
        if (nbgno .ne. 1) then
            call utmess('F', 'MODELISA6_50')
        end if
        call getvem(noma, 'GROUP_NO', motfac, 'GROUP_NO_2', iocc, &
                    nbgno, nogrno, nno)
        call jelira(jexnom(grnoma, nogrno), 'LONUTI', n1)
        if (n1 .ne. 1) then
            call utmess('F', 'MODELISA6_43', sk=nogrno)
        else
            call jeveuo(jexnom(grnoma, nogrno), 'L', jgro)
            in = zi(jgro+1-1)
            noepou = int_to_char8(in)
        end if
    end if
!
! --- RECUPERATION DU VECTEUR ORIENTANT LA POUTRE ET DIRIGE
! --- DE LA PARTIE COQUE VERS LA PARTIE POUTRE :
!     ----------------------------------------
    call getvr8(motfac, 'AXE_POUTRE', iocc=iocc, nbval=3, vect=axepou, &
                nbret=naxe)
    if (naxe .eq. 0) then
        call utmess('F', 'MODELISA6_51')
    end if
!
    xnorm = sqrt(axepou(1)*axepou(1)+axepou(2)*axepou(2)+axepou(3)*axepou(3))
!
    if (xnorm .le. r8prem()) then
        call utmess('F', 'MODELISA6_52')
    end if
!
    axepou(1) = axepou(1)/xnorm
    axepou(2) = axepou(2)/xnorm
    axepou(3) = axepou(3)/xnorm
!
! --- NOTATION DANS LA CARTE DE NOM '&&RAPOCO.CAXE_POU' DES
! --- COORDONNEES DU VECTEUR UNITAIRE ORIENTANT LA POUTRE :
!     ---------------------------------------------------
!
    nocmp(1) = 'X'
    nocmp(2) = 'Y'
    nocmp(3) = 'Z'
!
    call mecact('V', '&&RAPOCO.CAXE_POU', 'LIGREL', ligrel, 'GEOM_R', &
                ncmp=3, lnomcmp=nocmp, vr=axepou)
!
! --- RECUPERATION DES CARACTERISTIQUES ELEMENTAIRES :
!     ----------------------------------------------
    call getvid(motfac, 'CARA_ELEM', iocc=iocc, scal=cara, nbret=ncara)
    if (ncara .eq. 0) then
        call utmess('F', 'MODELISA6_53')
    end if
!
! ---  NUMERO DU NOEUD POUTRE A LIER :
!      -----------------------------
    numnop = char8_to_int(noepou)
!
! ---  COORDONNEES DU NOEUD POUTRE :
!      ---------------------------
    xpou = vale(3*(numnop-1)+1)
    ypou = vale(3*(numnop-1)+2)
    zpou = vale(3*(numnop-1)+3)
!
! --- VERIFICATION DU FAIT QUE LES NOEUDS DE LISNOE (DONC
! --- APPARTENANT A LA COQUE)  PORTENT LES DDLS DE ROTATION :
!     -----------------------------------------------------
    do i = 1, lonlis
! ---     NUMERO DU NOEUD COURANT DE LA LISTE
        ino = char8_to_int(zk8(ilisno+i-1))
!
        dg = prnm((ino-1)*nbec+1)
        do j = 4, 6
            icmp(j) = indik8(nomcmp, cmp(j), 1, nddla)
            if (.not. exisdg([dg], icmp(j))) then
                valk(1) = zk8(ilisno+i-1)
                valk(2) = cmp(j)
                call utmess('F', 'MODELISA6_54', nk=2, valk=valk)
            end if
        end do
    end do
!
! --- VERIFICATION DU FAIT QUE LE NOEUD POUTRE A RACCORDER PORTE
! --- LES 3 DDLS DE TRANSLATION ET LES 3 DDLS DE ROTATION :
!     ---------------------------------------------------
    dg = prnm((numnop-1)*nbec+1)
    do j = 1, 6
        icmp(j) = indik8(nomcmp, cmp(j), 1, nddla)
        if (.not. exisdg([dg], icmp(j))) then
            valk(1) = noepou
            valk(2) = cmp(j)
            call utmess('F', 'MODELISA6_45', nk=2, valk=valk)
        end if
    end do
!
! --- CALCUL SUR CHAQUE ELEMENT DE BORD A RELIER A LA POUTRE
! --- DES CARACTERISTIQUES GEOMETRIQUES SUIVANTES :
! ---  SOMME/S_ELEMENT(1,X,Y,Z,X*X,Y*Y,Z*Z,X*Y,X*Z,Y*Z)DS
!     ---------------------------------------------------
    lpain(1) = 'PGEOMER'
    lchin(1) = noma//'.COORDO'
    lpain(2) = 'PCACOQU'
    lchin(2) = cara//'.CARCOQUE'
    lpain(3) = 'PCAORIE'
    lchin(3) = '&&RAPOCO.CAXE_POU'
    lpaout(1) = 'PCASECT'
    lchout(1) = '&&RAPOCO.PSECT'
!
    call calcul('S', 'CARA_SECT_POUT3', ligrel, 3, lchin, &
                lpain, 1, lchout, lpaout, 'V', &
                'OUI')
!
! --- VECTEUR DES QUANTITES GEOMETRIQUES PRECITEES SOMMEES
! --- SUR LA SURFACE DE RACCORD, CES QUANTITES SERONT NOTEES :
! ---  A1 = S,AX,AY,AZ,AXX,AYY,AZZ,AXY,AXZ,AYZ
!     ----------------------------------------
    AS_ALLOCATE(vr=inertie_raccord, size=10)
!
! --- SOMMATION DES QUANTITES GEOMETRIQUES ELEMENTAIRES
! --- DANS LE VECTEUR &&RAPOCO.INERTIE_RACCORD :
!     ----------------------------------------
    call mesomm(lchout(1), 10, vr=inertie_raccord)
!
    s = inertie_raccord(1)
    ax = inertie_raccord(2)
    ay = inertie_raccord(3)
    az = inertie_raccord(4)
    axx = inertie_raccord(5)
    ayy = inertie_raccord(6)
    azz = inertie_raccord(7)
    axy = inertie_raccord(8)
    axz = inertie_raccord(9)
    ayz = inertie_raccord(10)
!
    if (abs(s) .lt. r8prem()) then
        call utmess('F', 'MODELISA6_55')
    end if
    s1 = 1.0d0/s
!
! --- COORDONNEES DU CENTRE GEOMETRIQUE G DE LA SECTION DE RACCORD
! --- XG = AX/S, YG = AY/S, ZG = AZ/S :
!     -------------------------------
    xg = s1*ax
    yg = s1*ay
    zg = s1*az
! --- VERIFICATION DE L'IDENTITE GEOMETRIQUE DE G AVEC LE
! --- NOEUD POUTRE A RACCORDER :
!     ------------------------
    dnorme = sqrt((xpou-xg)*(xpou-xg)+(ypou-yg)*(ypou-yg)+(zpou-zg)*(zpou-zg))/sqrt(s/pi)
    if (dnorme .gt. eps) then
        valr(1) = xg
        valr(2) = yg
        valr(3) = zg
        valr(4) = xpou
        valr(5) = ypou
        valr(6) = zpou
        valr(7) = eps*100.0d0
        valr(8) = sqrt(s/pi)
        valr(9) = dnorme
        valk(1) = option
        vali(1) = iocc
        call utmess('A', 'CALCULEL3_80', sk=valk(1), si=vali(1), nr=9, &
                    valr=valr)
    end if
!
! --- CALCUL DU TENSEUR D'INERTIE EN G, CE TENSEUR EST SYMETRIQUE :
! --- ON CALCULE LES COMPOSANTES DE LA PARTIE SUPERIEURE PAR LIGNE
!     ------------------------------------------------------------
!
! ---    IGXX = AYY + AZZ -S*(YG*YG+ZG*ZG)
    ig(1) = ayy+azz-s*(yg*yg+zg*zg)
! ---    IGXY = -AXY + S*XG*YG
    ig(2) = -axy+s*xg*yg
! ---    IGXZ = -AXZ + S*XG*ZG
    ig(3) = -axz+s*xg*zg
! ---    IGYY = AZZ + AXX -S*(ZG*ZG+XG*XG)
    ig(4) = azz+axx-s*(zg*zg+xg*xg)
! ---    IGYZ = -AYZ + S*YG*ZG
    ig(5) = -ayz+s*yg*zg
! ---    IGZZ = AXX + AYY -S*(XG*XG+YG*YG)
    ig(6) = axx+ayy-s*(xg*xg+yg*yg)
!
! --- NOTATION DANS LA CARTE DE NOM '&&RAPOCO.CAORIGE' DES
! --- COORDONNEES DU CENTRE GEOMETRIQUE G DE LA SECTION DE RACCORD
!     ------------------------------------------------------------
!
    nocmp(1) = 'X'
    nocmp(2) = 'Y'
    nocmp(3) = 'Z'
!
    coorig(1) = xg
    coorig(2) = yg
    coorig(3) = zg
!
    call mecact('V', '&&RAPOCO.CAORIGE', 'LIGREL', ligrel, 'GEOM_R', &
                ncmp=3, lnomcmp=nocmp, vr=coorig)
!
! --- DETERMINATION DE 2 LISTES  DE VECTEURS PAR ELEMENT PRENANT
! --- LEURS VALEURS AUX NOEUDS DES ELEMENTS.
! --- LA PREMIERE LISTE DE NOM 'VECT_EINI' A POUR VALEURS AU NOEUD
! --- I D'UN ELEMENT :
! ---  SOMME/S_ELEMENT(E1(1)*NI,E1(2)*NI,E1(3)*NI,
! ---                  E2(1)*NI,E2(2)*NI,E2(3)*NI)DS
! --- OU E1 EST UN VECTEUR UNITAIRE PERPENDICULAIRE A L'ELEMENT
! --- DE BORD ORIENTE DE LA COQUE VERS LA POUTRE ET
! --- E2 EST LE VECTEUR TANGENT A LA FIBRE MOYENNE DE L'ELEMENT DE BORD
! --- LA SECONDE LISTE DE NOM 'VECT_XYZNI' A POUR VALEURS AU NOEUD
! --- I D'UN ELEMENT :
! ---  SOMME/S_ELEMENT(X*NI,Y*NI,Z*NI,NI,0,0)DS
! ---  AVEC X = XM - XG = NJ*XJ - XG
! ---       Y = YM - YG = NJ*YJ - YG
! ---       Z = ZM - ZG = NJ*ZJ - ZG
!     ------------------------------
    lpain(1) = 'PGEOMER'
    lchin(1) = noma//'.COORDO'
    lpain(2) = 'PORIGIN'
    lchin(2) = '&&RAPOCO.CAORIGE'
    lpain(3) = 'PCACOQU'
    lchin(3) = cara//'.CARCOQUE'
    lpain(4) = 'PCAORIE'
    lchin(4) = '&&RAPOCO.CAXE_POU'
    lpaout(1) = 'PVECTU1'
    lpaout(2) = 'PVECTU2'
    lchout(1) = '&&RAPOCO.VECT_XYZNI'
    lchout(2) = '&&RAPOCO.VECT2'
!
    call calcul('S', 'CARA_SECT_POUT4', ligrel, 4, lchin, &
                lpain, 2, lchout, lpaout, 'V', &
                'OUI')
!
! --- CREATION DES .RERR DES VECTEURS EN SORTIE DE CALCUL
!     --------------------------------------------------------
!
    call memare('V', '&&RAPOCO', mod, 'CHAR_MECA')
!
! --- ASSEMBLAGE DE LCHOUT(1) DANS LE CHAMNO DE NOM 'CH_DEPL_1'
!     ---------------------------------------------------------
    call jedetr('&&RAPOCO           .RELR')
    call reajre('&&RAPOCO', lchout(1), 'V')
    call assvec('V', '&&RAPOCO.CH_DEPL_01', 1, '&&RAPOCO           .RELR', [1.d0], &
                numddl)
!
    vale1 = '&&RAPOCO.CH_DEPL_01.VALE'
    call jeveuo(vale1, 'L', idch1)
!
!
! --- ASSEMBLAGE DE LCHOUT(2) DANS LE CHAMNO DE NOM 'CH_DEPL_1'
!     ---------------------------------------------------------
    call jedetr('&&RAPOCO           .RELR')
    call reajre('&&RAPOCO', lchout(2), 'V')
!
    call assvec('V', '&&RAPOCO.CH_DEPL_02', 1, '&&RAPOCO           .RELR', [1.d0], &
                numddl)
!
    vale2 = '&&RAPOCO.CH_DEPL_02.VALE'
    call jeveuo(vale2, 'L', idch2)
!
! --- CREATION DES TABLEAUX DE TRAVAIL NECESSAIRES A L'AFFECTATION
! --- DE LISREL
!     ------------------------------------------------------------
!
! ---     MAJORANT DU NOMBRE DE TERMES DANS UNE RELATION
    nbterm = 5*lonlis+3
! ---     VECTEUR DU NOM DES NOEUDS
    AS_ALLOCATE(vk8=lisno, size=nbterm)
! ---     VECTEUR DU NOM DES DDLS
    AS_ALLOCATE(vk8=lisddl, size=nbterm)
! ---     VECTEUR DES COEFFICIENTS REELS
    AS_ALLOCATE(vr=coer, size=nbterm)
! ---     VECTEUR DES COEFFICIENTS COMPLEXES
    AS_ALLOCATE(vc=coec, size=nbterm)
! ---     VECTEUR DES DIRECTIONS DES DDLS A CONTRAINDRE
    AS_ALLOCATE(vr=direct, size=6*nbterm)
! ---     VECTEUR DES DIMENSIONS DE CES DIRECTIONS
    AS_ALLOCATE(vi=dime, size=nbterm)
!
! ---    RELATIONS ENTRE LES NOEUDS DU MASSIF ET LE NOEUD POUTRE
!        -------------------------------------------------------
! ---    PREMIER GROUPE DE RELATIONS TRADUISANT :
! ---      SOMME/S_RACCORD(U_COQUE) = S_RACCORD*U_NOEUD_POUTRE
!         ----------------------------------------------------
!
! ---    PREMIERE RELATION :
!     -S.DX(NOEUD_POUTRE) + (SOMME/S_RACCORD(NI.DS)).DX(NOEUD_I) = 0
!     --------------------------------------------------------------
!
    nbterm = lonlis+1
! ---    BOUCLE SUR LES NOEUDS DES MAILLES DE BORD DE LA PARTIE COQUE
    do i = 1, lonlis
        ino = char8_to_int(zk8(ilisno+i-1))
! ---    ADRESSE DE LA PREMIERE COMPOSANTE DU NOEUD INO DANS LES
! ---    CHAMNO (SOMME/S_RACCORD(NI.DS))
        ival = zi(iaprno+(ino-1)*(nbec+2)+1-1)-1
!
        lisno(i) = zk8(ilisno+i-1)
        lisddl(i) = 'DX'
        coer(i) = zr(idch1-1+ival+4)
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
! ---    DEUXIEME RELATION :
!     -S.DY(NOEUD_POUTRE) + (SOMME/S_RACCORD(NI.DS)).DY(NOEUD_I) = 0
!     --------------------------------------------------------------
!
    nbterm = lonlis+1
! ---   BOUCLE SUR LES NOEUDS DES MAILLES DE BORD DE LA PARTIE COQUE
    do i = 1, lonlis
        ino = char8_to_int(zk8(ilisno+i-1))
! ---    ADRESSE DE LA PREMIERE COMPOSANTE DU NOEUD INO DANS LES
! ---    CHAMNO (SOMME/S_RACCORD(NI.DS))
        ival = zi(iaprno+(ino-1)*(nbec+2)+1-1)-1
!
        lisno(i) = zk8(ilisno+i-1)
        lisddl(i) = 'DY'
        coer(i) = zr(idch1-1+ival+4)
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
! ---    TROISIEME RELATION :
!     -S.DZ(NOEUD_POUTRE) + (SOMME/S_RACCORD(NI.DS)).DZ(NOEUD_I) = 0
!     --------------------------------------------------------------
!
    nbterm = lonlis+1
! --- BOUCLE SUR LES NOEUDS DES MAILLES DE BORD DE LA PARTIE COQUE
    do i = 1, lonlis
        ino = char8_to_int(zk8(ilisno+i-1))
! ---    ADRESSE DE LA PREMIERE COMPOSANTE DU NOEUD INO DANS LES
! ---    CHAMNO (SOMME/S_RACCORD(NI.DS))
        ival = zi(iaprno+(ino-1)*(nbec+2)+1-1)-1
!
        lisno(i) = zk8(ilisno+i-1)
        lisddl(i) = 'DZ'
        coer(i) = zr(idch1-1+ival+4)
    end do
!
    lisno(1+lonlis+1-1) = noepou
    lisddl(1+lonlis+1-1) = 'DZ'
    coer(1+lonlis+1-1) = -s
!
    call afrela(coer, coec, lisddl, lisno, dime, &
                direct, nbterm, beta, betac, betaf, &
                typcoe, typval, 0.d0, lisrel)
    call imprel(motfac, nbterm, coer, lisddl, lisno, &
                beta)
!
! ---    DEUXIEME GROUPE DE RELATIONS TRADUISANT :
! ---      SOMME/S_RACCORD(GM X U_COQUE) = I.OMEGA(NOEUD_POUTRE)
!         ------------------------------------------------------
! ---    QUATRIEME RELATION :
!
    nbterm = 5*lonlis+3
! ---    BOUCLE SUR LES NOEUDS DES MAILLES DE BORD DE LA PARTIE COQUE
    do i = 1, lonlis
        ino = char8_to_int(zk8(ilisno+i-1))
! ---    ADRESSE DE LA PREMIERE COMPOSANTE DU NOEUD INO DANS LES CHAMNO
        ival = zi(iaprno+(ino-1)*(nbec+2)+1-1)-1
!
        lisno(1+5*(i-1)+1-1) = zk8(ilisno+i-1)
        lisno(1+5*(i-1)+2-1) = zk8(ilisno+i-1)
        lisno(1+5*(i-1)+3-1) = zk8(ilisno+i-1)
        lisno(1+5*(i-1)+4-1) = zk8(ilisno+i-1)
        lisno(1+5*(i-1)+5-1) = zk8(ilisno+i-1)
!
        lisddl(1+5*(i-1)+1-1) = 'DZ'
        lisddl(1+5*(i-1)+2-1) = 'DY'
        lisddl(1+5*(i-1)+3-1) = 'DRX'
        lisddl(1+5*(i-1)+4-1) = 'DRY'
        lisddl(1+5*(i-1)+5-1) = 'DRZ'
!
        coer(1+5*(i-1)+1-1) = zr(idch1-1+ival+2)
        coer(1+5*(i-1)+2-1) = -zr(idch1-1+ival+3)
        coer(1+5*(i-1)+3-1) = zr(idch2-1+ival+1)
        coer(1+5*(i-1)+4-1) = zr(idch2-1+ival+2)
        coer(1+5*(i-1)+5-1) = zr(idch2-1+ival+3)
    end do
!
    lisno(1+5*lonlis+1-1) = noepou
    lisno(1+5*lonlis+2-1) = noepou
    lisno(1+5*lonlis+3-1) = noepou
!
    lisddl(1+5*lonlis+1-1) = 'DRX'
    lisddl(1+5*lonlis+2-1) = 'DRY'
    lisddl(1+5*lonlis+3-1) = 'DRZ'
!
    coer(1+5*lonlis+1-1) = -ig(1)
    coer(1+5*lonlis+2-1) = -ig(2)
    coer(1+5*lonlis+3-1) = -ig(3)
!
    call afrela(coer, coec, lisddl, lisno, dime, &
                direct, nbterm, beta, betac, betaf, &
                typcoe, typval, 0.d0, lisrel)
    call imprel(motfac, nbterm, coer, lisddl, lisno, &
                beta)
!
! ---    CINQUIEME RELATION :
!
    nbterm = 5*lonlis+3
! ---    BOUCLE SUR LES NOEUDS DES MAILLES DE BORD DE LA PARTIE COQUE
    do i = 1, lonlis
        ino = char8_to_int(zk8(ilisno+i-1))
! ---    ADRESSE DE LA PREMIERE COMPOSANTE DU NOEUD INO DANS LES CHAMNO
        ival = zi(iaprno+(ino-1)*(nbec+2)+1-1)-1
!
        lisno(1+5*(i-1)+1-1) = zk8(ilisno+i-1)
        lisno(1+5*(i-1)+2-1) = zk8(ilisno+i-1)
        lisno(1+5*(i-1)+3-1) = zk8(ilisno+i-1)
        lisno(1+5*(i-1)+4-1) = zk8(ilisno+i-1)
        lisno(1+5*(i-1)+5-1) = zk8(ilisno+i-1)
!
        lisddl(1+5*(i-1)+1-1) = 'DX'
        lisddl(1+5*(i-1)+2-1) = 'DZ'
        lisddl(1+5*(i-1)+3-1) = 'DRX'
        lisddl(1+5*(i-1)+4-1) = 'DRY'
        lisddl(1+5*(i-1)+5-1) = 'DRZ'
!
        coer(1+5*(i-1)+1-1) = zr(idch1-1+ival+3)
        coer(1+5*(i-1)+2-1) = -zr(idch1-1+ival+1)
        coer(1+5*(i-1)+3-1) = zr(idch2-1+ival+2)
        coer(1+5*(i-1)+4-1) = zr(idch2-1+ival+4)
        coer(1+5*(i-1)+5-1) = zr(idch2-1+ival+5)
    end do
!
    lisno(1+5*lonlis+1-1) = noepou
    lisno(1+5*lonlis+2-1) = noepou
    lisno(1+5*lonlis+3-1) = noepou
!
    lisddl(1+5*lonlis+1-1) = 'DRX'
    lisddl(1+5*lonlis+2-1) = 'DRY'
    lisddl(1+5*lonlis+3-1) = 'DRZ'
!
    coer(1+5*lonlis+1-1) = -ig(2)
    coer(1+5*lonlis+2-1) = -ig(4)
    coer(1+5*lonlis+3-1) = -ig(5)
!
    call afrela(coer, coec, lisddl, lisno, dime, &
                direct, nbterm, beta, betac, betaf, &
                typcoe, typval, 0.d0, lisrel)
    call imprel(motfac, nbterm, coer, lisddl, lisno, &
                beta)
!
! ---    SIXIEME RELATION :
!
    nbterm = 5*lonlis+3
! ---    BOUCLE SUR LES NOEUDS DES MAILLES DE BORD DE LA PARTIE COQUE
    do i = 1, lonlis
        ino = char8_to_int(zk8(ilisno+i-1))
! ---    ADRESSE DE LA PREMIERE COMPOSANTE DU NOEUD INO DANS LES CHAMNO
        ival = zi(iaprno+(ino-1)*(nbec+2)+1-1)-1
!
        lisno(1+5*(i-1)+1-1) = zk8(ilisno+i-1)
        lisno(1+5*(i-1)+2-1) = zk8(ilisno+i-1)
        lisno(1+5*(i-1)+3-1) = zk8(ilisno+i-1)
        lisno(1+5*(i-1)+4-1) = zk8(ilisno+i-1)
        lisno(1+5*(i-1)+5-1) = zk8(ilisno+i-1)
!
        lisddl(1+5*(i-1)+1-1) = 'DY'
        lisddl(1+5*(i-1)+2-1) = 'DX'
        lisddl(1+5*(i-1)+3-1) = 'DRX'
        lisddl(1+5*(i-1)+4-1) = 'DRY'
        lisddl(1+5*(i-1)+5-1) = 'DRZ'
!
        coer(1+5*(i-1)+1-1) = zr(idch1-1+ival+1)
        coer(1+5*(i-1)+2-1) = -zr(idch1-1+ival+2)
        coer(1+5*(i-1)+3-1) = zr(idch2-1+ival+3)
        coer(1+5*(i-1)+4-1) = zr(idch2-1+ival+5)
        coer(1+5*(i-1)+5-1) = zr(idch2-1+ival+6)
    end do
!
    lisno(1+5*lonlis+1-1) = noepou
    lisno(1+5*lonlis+2-1) = noepou
    lisno(1+5*lonlis+3-1) = noepou
!
    lisddl(1+5*lonlis+1-1) = 'DRX'
    lisddl(1+5*lonlis+2-1) = 'DRY'
    lisddl(1+5*lonlis+3-1) = 'DRZ'
!
    coer(1+5*lonlis+1-1) = -ig(3)
    coer(1+5*lonlis+2-1) = -ig(5)
    coer(1+5*lonlis+3-1) = -ig(6)
!
    call afrela(coer, coec, lisddl, lisno, dime, &
                direct, nbterm, beta, betac, betaf, &
                typcoe, typval, 0.d0, lisrel)
    call imprel(motfac, nbterm, coer, lisddl, lisno, &
                beta)
!
    if ((option .eq. 'COQ_TUYA')) then
        call racotu(zi(iaprno), lonlis, zk8(ilisno), noepou, noma, &
                    ligrel, mod, cara, numddl, lisrel, coorig)
    end if
!
! --- DESTRUCTION DES OBJETS DE TRAVAIL
!     ---------------------------------
    call detrsd('LIGREL', ligrel)
    call jedetr('&&RAPOCO.LISTE_NOEUDS')
    call jedetr('&&RAPOCO.LISTE_MAILLES')
    call detrsd('CARTE', '&&RAPOCO.CAXE_POU')
    call detrsd('CHAMP_GD', '&&RAPOCO.PSECT')
    AS_DEALLOCATE(vr=inertie_raccord)
    call detrsd('CARTE', '&&RAPOCO.CAORIGE')
    call detrsd('RESUELEM', '&&RAPOCO.VECT2')
    call detrsd('RESUELEM', '&&RAPOCO.VECT_XYZNI')
    call jedetr('&&RAPOCO           .RELR')
    call jedetr('&&RAPOCO           .RERR')
    AS_DEALLOCATE(vk8=lisno)
    AS_DEALLOCATE(vk8=lisddl)
    AS_DEALLOCATE(vr=coer)
    AS_DEALLOCATE(vc=coec)
    AS_DEALLOCATE(vr=direct)
    AS_DEALLOCATE(vi=dime)
    call detrsd('CHAMP_GD', '&&RAPOCO.CH_DEPL_01')
    call detrsd('CHAMP_GD', '&&RAPOCO.CH_DEPL_02')
!
130 continue
    call jedema()
end subroutine
