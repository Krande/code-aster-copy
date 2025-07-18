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

subroutine flust3(melflu, typflu, base, nuor, amor, &
                  freq, masg, fact, vite, nbm, &
                  npv, nivpar, nivdef)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8pi.h"
#include "asterc/r8dgrd.h"
#include "asterfort/codent.h"
#include "asterfort/dismoi.h"
#include "asterfort/extmod.h"
#include "asterfort/fluimp.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/mefcen.h"
#include "asterfort/mefgmn.h"
#include "asterfort/mefgri.h"
#include "asterfort/mefint.h"
#include "asterfort/mefist.h"
#include "asterfort/mefrac.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rslipa.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    integer(kind=8) :: nbm, npv, nivpar, nivdef, nuor(*)
    real(kind=8) :: amor(*), freq(*), masg(*), vite(*), fact(*)
    character(len=8) :: typflu, base
    character(len=19) :: melflu
!  CALCUL DES PARAMETRES DE COUPLAGE FLUIDE-STRUCTURE POUR UNE
!  CONFIGURATION DE TYPE "FAISCEAU DE TUBES SOUS ECOULEMENT AXIAL"
!  OPERATEUR APPELANT : CALC_FLUI_STRU , OP0144
!-----------------------------------------------------------------------
!  IN : MELFLU : NOM DU CONCEPT DE TYPE MELASFLU PRODUIT
!  IN : TYPFLU : NOM DU CONCEPT DE TYPE TYPE_FLUI_STRU DEFINISSANT LA
!                CONFIGURATION ETUDIEE
!  IN : BASE   : NOM DU CONCEPT DE TYPE MODE_MECA DEFINISSANT LA BASE
!                MODALE DU SYSTEME AVANT PRISE EN COMPTE DU COUPLAGE
!  IN : NUOR   : LISTE DES NUMEROS D'ORDRE DES MODES SELECTIONNES POUR
!                LE COUPLAGE (PRIS DANS LE CONCEPT MODE_MECA)
!  IN : AMOR   : LISTE DES AMORTISSEMENTS REDUITS MODAUX INITIAUX
!  IN : VITE   : LISTE DES VITESSES D'ECOULEMENT ETUDIEES
!  IN : NBM    : NOMBRE DE MODES PRIS EN COMPTE POUR LE COUPLAGE
!  IN : NPV    : NOMBRE DE VITESSES D'ECOULEMENT
!  IN : NIVPAR : NIVEAU D'IMPRESSION DANS LE FICHIER RESULTAT POUR LES
!                PARAMETRES DU COUPLAGE (FREQ,AMOR)
!  IN : NIVDEF : NIVEAU D'IMPRESSION DANS LE FICHIER RESULTAT POUR LES
!                DEFORMEES MODALES
!  OUT: FREQ   : FREQUENCES ET AMORTISSEMENTS REDUITS MODAUX PERTUBES
!                PAR L'ECOULEMENT
!  OUT: MASG   : MASSES GENERALISEES DES MODES PERTURBES, SUIVANT LA
!                DIRECTION CHOISIE PAR L'UTILISATEUR
!  OUT: FACT   : PSEUDO FACTEUR DE PARTICIPATION
!-----------------------------------------------------------------------
!     ------------------------------------------------------------------
!
    character(len=1) :: k1bid
    character(len=8) :: mailla, k8b, promas, provis, numgno
    character(len=14) :: numddl
    character(len=19) :: caelem
    character(len=24) :: fsvi, fsvk, fsvr, fsgm, fscr, fsgr
    character(len=24) :: matria, nomrac
    character(len=24) :: grpno
    real(kind=8) :: som(9), ru, alpha, coorpe(3), don(5), carac(2)
    integer(kind=8) :: irot(3), iddl(6), ndim(14)
    aster_logical :: calcul(2)
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iadmag, iadnog, iaxe, icdg, icency
    integer(kind=8) :: icpg, idefm, iencei, iepaig, iequiv, ifpart, ifreqi
    integer(kind=8) :: ifscr, ifsgm, ifsgr, ifsvi, ifsvk, ifsvr, igreq
    integer(kind=8) :: ikn, ilargg, ilongg, im, imataa, imatma
    integer(kind=8) :: imatra, inbmag, inbneq, inbnog, inomcy, inomeq, inum
    integer(kind=8) :: ior, iphix, iphiy, ireq, irint
    integer(kind=8) :: irugg, ixint, iyint, iz, izg, izint, j
    integer(kind=8) :: jtypg, lfacx, lmasg, lrigg, n, n1
    integer(kind=8) :: nbcyl, nbddl, nbfin, nbgrma, nbgrmx, nbgrp, nbgtot
    integer(kind=8) :: nbnoe, nbtot, nbtron, nbz, ncoor, ndir, neq
    integer(kind=8) :: nima, nima2, nn, nt, ntypg, numno1, numno2
!
    real(kind=8) :: ang, epsit, g, pi, x1
    real(kind=8) :: x2, y1, y2, z1, z2, zmax, zmin
    real(kind=8), pointer :: vale(:) => null()
!
!-----------------------------------------------------------------------
    data iddl/1, 2, 3, 4, 5, 6/
!
!-----------------------------------------------------------------------
!
    call jemarq()
    epsit = 1.d-5
    pi = r8pi()
!
! --- 1.RECUPERATION DES INFORMATIONS APPORTEES PAR LE CONCEPT ---
! ---   TYPE_FLUI_STRU                                         ---
!
! --- 1.1.GROUPE D EQUIVALENCE  - AXE DIRECTEUR DU TUBE - TYPE
! ---     D ENCEINTE - NOMBRE DE GROUPES DE NOEUDS
!
    fsvi = typflu//'           .FSVI'
    call jeveuo(fsvi, 'L', ifsvi)
    iequiv = zi(ifsvi)
    iaxe = zi(ifsvi+1)
    iencei = zi(ifsvi+2)
    nbgrma = zi(ifsvi+3)
    ntypg = zi(ifsvi+4)
!
! --- 1.2.NOMBRE DE CYLINDRES, DE GROUPES D EQUIVALENCE, ET DE
! ---       CYLINDRES REELS PAR GROUPE D EQUIVALENCE
!
    nbgrp = nbgrma
    nbgtot = 0
    if (iequiv .eq. 1) then
        nbcyl = zi(ifsvi+5)
        if (ntypg .ne. 0) then
            nbgtot = zi(ifsvi+6+nbgrma)
        end if
        call wkvect('&&FLUST3.TMP.NBCC', 'V V I', nbgrp, inbneq)
        do i = 1, nbgrp
            zi(inbneq+i-1) = zi(ifsvi+i+5)
        end do
    else if (iequiv .eq. 0) then
        nbcyl = nbgrma
        if (ntypg .ne. 0) then
            nbgtot = zi(ifsvi+5)
        end if
    end if
!
! --- 1.3.CONCEPTS DE TYPE FONCTION DEFINISSANT LE PROFIL DE MASSE
! ---     VOLUMIQUE ET LE PROFIL DE VISCOSITE CINEMATIQUE DU FLUIDE
!
    fsvk = typflu//'           .FSVK'
    call jeveuo(fsvk, 'L', ifsvk)
    promas = zk8(ifsvk)
    provis = zk8(ifsvk+1)
!
! --- 1.4.NOM DU CONCEPT DE TYPE CARA_ELEM
!
    if (iequiv .eq. 0) caelem = zk8(ifsvk+2)
!
! --- 1.5.NORME ET ORIENTATION DU VECTEUR PESANTEUR - RUGOSITE
!
    fsvr = typflu//'           .FSVR'
    call jeveuo(fsvr, 'L', ifsvr)
    g = zr(ifsvr)
    coorpe(1) = zr(ifsvr+1)
    coorpe(2) = zr(ifsvr+2)
    coorpe(3) = zr(ifsvr+3)
    ru = zr(ifsvr+4)
!
! --- 1.6.DONNEES GEOMETRIQUES DES ENCEINTES
!
    if (iencei .eq. 1) then
        som(1) = zr(ifsvr+5)
        som(2) = zr(ifsvr+6)
        som(3) = zr(ifsvr+7)
        ikn = 8
    else if (iencei .eq. 2) then
        do i = 1, 5
            don(i) = zr(ifsvr+i+4)
        end do
        ikn = 10
!
! ---    CALCUL DES COORDONNEES DES QUATRES SOMMETS DE L'ENCEINTE
!
        ang = don(5)*r8dgrd()
        som(1) = don(1)+(don(3)*cos(ang)-don(4)*sin(ang))/2
        som(2) = don(2)+(don(3)*sin(ang)+don(4)*cos(ang))/2
        som(3) = don(1)-(don(3)*cos(ang)+don(4)*sin(ang))/2
        som(4) = don(2)-(don(3)*sin(ang)-don(4)*cos(ang))/2
        som(5) = don(1)-(don(3)*cos(ang)-don(4)*sin(ang))/2
        som(6) = don(2)-(don(3)*sin(ang)+don(4)*cos(ang))/2
        som(7) = don(1)+(don(3)*cos(ang)+don(4)*sin(ang))/2
        som(8) = don(2)+(don(3)*sin(ang)-don(4)*cos(ang))/2
    else
        call utmess('F', 'ALGELINE_44')
    end if
!
! --- 1.7.RAYONS DES TUBES DU FAISCEAU REEL ASSOCIES AUX TUBES
! ---     EQUIVALENTS
!
    if (iequiv .eq. 1) then
        call wkvect('&&FLUST3.TMP.REQ', 'V V R', nbgrp, ireq)
        do i = 1, nbgrp
            zr(ireq+i-1) = zr(ifsvr+ikn+i-1)
        end do
    else
        ireq = 1
    end if
!
! --- 1.8.NOM DES GROUPES DE MAILLES ASSOCIES AUX TUBES EQUIVALENTS,
! ---     OU NOM DE LA RACINE COMMUNE DES GROUPES DE NOEUDS,
! ---     OU NOM DES GROUPES DE NOEUDS DES TUBES REELS.
!
    fsgm = typflu//'           .FSGM'
    call jeveuo(fsgm, 'L', ifsgm)
    if (iequiv .eq. 1) then
! --     NOMS DES GROUPES CLASSES D EQUIVALENCE
        call wkvect('&&FLUST3.TMP.NOEQ', 'V V K24', nbcyl, inomeq)
        do i = 1, nbgrp
            zk24(inomeq+i-1) = zk24(ifsgm+i-1)
        end do
    else if (iequiv .eq. 0 .and. nbgrma .eq. 0) then
! --     NOM DE LA RACINE COMMUNE
        nomrac = zk24(ifsgm)
    else if (iequiv .eq. 0 .and. nbgrma .ne. 0) then
! --     NOM POUR CHAQUE CYLINDRE DE LA CLASSE D EQUIVALENCE A
! --     LAQUELLE IL APPARTIENT
        call wkvect('&&FLUST3.TMP.NOCY', 'V V K24', nbcyl, inomcy)
        do i = 1, nbcyl
            zk24(inomcy+i-1) = zk24(ifsgm+i-1)
        end do
    end if
!
! --- 1.9.COORDONNEES DES CENTRES DES TUBES DU FAISCEAU REEL
!
    if (iequiv .eq. 1) then
        fscr = typflu//'           .FSCR'
        call jeveuo(fscr, 'L', ifscr)
        call wkvect('&&FLUST3.TMP.CENT', 'V V R', 2*nbcyl, icency)
        do i = 1, nbcyl*2
            zr(icency+i-1) = zr(ifscr+i-1)
        end do
    else
        icency = 1
    end if
!
!--- 1.10. CARACTERISTIQUES DES GRILLES
!
!     --- INITIALISATION DES POINTEURS
    izg = 1
    ilongg = 1
    ilargg = 1
    iepaig = 1
    icdg = 1
    icpg = 1
    irugg = 1
    if (ntypg .ne. 0) then
        call wkvect('&&FLUST3.TMP.TYPG', 'V V I', nbgtot, jtypg)
        if (iequiv .eq. 1) then
            do i = 1, nbgtot
                zi(jtypg+i-1) = zi(ifsvi+i+6+nbgrma)
            end do
        else if (iequiv .eq. 0) then
            do i = 1, nbgtot
                zi(jtypg+i-1) = zi(ifsvi+i+5)
            end do
        end if
!
        fsgr = typflu//'           .FSGR'
        call jeveuo(fsgr, 'L', ifsgr)
        call wkvect('&&FLUST3.TMP.GRIL', 'V V R', nbgtot+6*ntypg, izg)
        do i = 1, nbgtot+6*ntypg
            zr(izg+i-1) = zr(ifsgr+i-1)
        end do
        ilongg = izg+nbgtot
        ilargg = ilongg+ntypg
        iepaig = ilargg+ntypg
        icdg = iepaig+ntypg
        icpg = icdg+ntypg
        irugg = icpg+ntypg
    else
        jtypg = 1
    end if
!
!
!-----RECUPERATION D'INFORMATIONS POUR CREATION DES OBJETS DE TRAVAIL
!
!
    call dismoi('REF_RIGI_PREM', base, 'RESU_DYNA', repk=matria)
    call dismoi('NOM_NUME_DDL', matria, 'MATR_ASSE', repk=numddl)
    call dismoi('NB_EQUA', matria, 'MATR_ASSE', repi=neq)
    call dismoi('NOM_MAILLA', matria, 'MATR_ASSE', repk=mailla)
    call jelira(mailla//'.COORDO    .VALE', 'LONMAX', nbnoe)
    nbnoe = nbnoe/3
    call jeveuo(mailla//'.COORDO    .VALE', 'L', vr=vale)
!
!
! --- 2.RECHERCHE DES NOMS DES GROUPES DE NOEUDS DONNES SOUS UNE
! ---   RACINE COMMUNE
!
    if (nbgrma .eq. 0) then
        call jelira(mailla//'.GROUPEMA', 'NOMUTI', nbgrmx)
        call wkvect('&&FLUST3.TMP.NGX', 'V V K24', nbgrmx, inomcy)
        call mefrac(mailla, nbgrmx, nomrac, nbgrma, zk24(inomcy))
        nbcyl = nbgrma
    end if
!
!
! --- 3.CONSTITUTION DES GROUPES D EQUIVALENCE
!
    call wkvect('&&FLUST3.TMP.GREQ', 'V V I', nbcyl, igreq)
    if (iequiv .eq. 0) then
        nbgrp = nbcyl
        inomeq = inomcy
        do i = 1, nbcyl
            zi(igreq+i-1) = i
        end do
    else if (iequiv .eq. 1) then
        nt = 0
        nn = 0
        do i = 1, nbgrp
            nt = nt+zi(inbneq+i-1)
            if (nt .gt. nbcyl) then
                call utmess('F', 'ALGELINE_45')
            end if
            do j = 1, zi(inbneq+i-1)
                nn = nn+1
                zi(igreq+nn-1) = i
            end do
        end do
        if (nt .ne. nbcyl) then
            call utmess('F', 'ALGELINE_45')
        end if
    end if
!
!
! --- 4.CREATION DE GROUPES DE NOEUDS A PARTIR DES GROUPES DE MAILLES
! ---   LES GROUPES DE NOEUDS ET DE MAILLES PORTENT LE MEME NOM
!
    if (iequiv .eq. 1) then
        call mefgmn(mailla, nbgrp, zk24(inomeq))
    else
        call mefgmn(mailla, nbcyl, zk24(inomcy))
    end if
!
!
! --- 5.DETERMINATION DU NOMBRE MAXIMUM DE NOEUDS PAR CYLINDRE
!
    nbz = 0
! --- CREATION DE TABLEAUX POUR LES ADRESSES DES NUMEROS DES NOEUDS
! --- DES CYLINDRES, ET POUR LE NOMBRE DE NOEUDS DE CHAQUE CYLINDRE
    call wkvect('&&FLUST3.TMP.ADR', 'V V I', nbgrp*4, iadnog)
    inbnog = iadnog+nbgrp
    iadmag = inbnog+nbgrp
    inbmag = iadmag+nbgrp
    do i = 1, nbgrp
        call codent(i, 'D0', numgno)
        grpno = '&&MEFGMN.'//numgno//'       '
        call jeveuo(grpno, 'L', zi(iadnog+i-1))
        call jelira(grpno, 'LONMAX', zi(inbnog+i-1), k1bid)
        call jeveuo(jexnom(mailla//'.GROUPEMA', zk24(inomeq+i-1)), 'L', zi(iadmag+i-1))
        call jelira(jexnom(mailla//'.GROUPEMA', zk24(inomeq+i-1)), 'LONMAX', zi(inbmag+i-1), &
                    k1bid)
        if (zi(inbnog+i-1) .gt. nbz) nbz = zi(inbnog+i-1)
    end do
!
!
! --- 6.VERIFICATION DE L AXE DIRECTEUR DU FAISCEAU. CREATION DU
! ---   TABLEAU IROT QUI DEFINIT UNE PERMUTATION CIRCULAIRE
! ---   PERMETTANT DE PASSER DU REPERE INITIAL AU REPERE AXIAL, DONT
! ---   L AXE Z A LA DIRECTION DE L AXE DIRECTEUR DU FAISCEAU.
!
    numno1 = zi(zi(iadnog))
    numno2 = zi(zi(iadnog)+1)
    x1 = vale(1+(numno1-1)*3)
    y1 = vale(1+(numno1-1)*3+1)
    z1 = vale(1+(numno1-1)*3+2)
    x2 = vale(1+(numno2-1)*3)
    y2 = vale(1+(numno2-1)*3+1)
    z2 = vale(1+(numno2-1)*3+2)
    if (abs(x1-x2) .lt. epsit .and. abs(y1-y2) .lt. epsit) then
        ndir = 3
        irot(1) = 1
        irot(2) = 2
        irot(3) = 3
    else if (abs(y1-y2) .lt. epsit .and. abs(z1-z2) .lt. epsit) then
        ndir = 1
        irot(1) = 2
        irot(2) = 3
        irot(3) = 1
    else if (abs(z1-z2) .lt. epsit .and. abs(x1-x2) .lt. epsit) then
        ndir = 2
        irot(1) = 3
        irot(2) = 1
        irot(3) = 2
    else
        call utmess('F', 'ALGELINE_46')
    end if
    if (ndir .ne. iaxe) then
        call utmess('F', 'ALGELINE_47')
    end if
!
!
! --- 7.COEFFICIENT DE PROPORTIONALITE DE LA PESENTEUR PAR
! ---   RAPPORT A LA VALEUR STANDARD (9.81) PROJETE SUR L AXE Z DU
! ---   REPERE AXIAL
!
    alpha = g*coorpe(irot(3))/9.81d0/sqrt(coorpe(1)**2+coorpe(2)**2+coorpe(3)**2)
!
!
! --- 8.RECHERCHE DES COORDONNEES DES CENTRES ET DES RAYONS DE CHAQUE
! ---   CYLINDRE
!
    if (iencei .eq. 2) then
! ---    PARAMETRE POUR LES CONDITIONS AUX LIMITES, TRAITEES PAR LA
! ---    METHODE DES IMAGES
! ---    - NOMBRE DE COURONNES D IMAGES COMPLETES   --> NIMA
! ---    - NOMBRE DE COURONNES D IMAGES SIMPLIFIEES --> NIMA2
        nima = 2
        nima2 = 8
    else
        nima = 0
        nima2 = 0
    end if
    nbtot = nbcyl*(2*nima+1)*(2*nima+1)
    nbfin = nbtot+4*(nima2)*(nima2+2*nima+1)
    ncoor = nbfin*(3+nbz)
    call wkvect('&&FLUST3.TMP.XYZR', 'V V R', ncoor, ixint)
    iyint = ixint+nbfin
    izint = iyint+nbfin
    irint = izint+nbz*nbfin
!
    call mefcen(caelem, iequiv, nbcyl, nbz, irot, &
                zi(iadnog), zi(inbnog), zi(iadmag), zi(igreq), vale, &
                zr(icency), zr(ireq), zr(ixint), zr(iyint), zr(izint), &
                zr(irint), nbgrp)
!
!
! --- 9.RECUPERATION DES RESULTATS DU CALCUL MODAL EN AIR
!
    call rslipa(base, 'FREQ', '&&FLUST3.LIFREQ', ifreqi, n1)
!
!
    call wkvect('&&FLUST3.TMP.IFAC', 'V V R', 3*nbm, ifpart)
    call wkvect('&&FLUST3.TMP.IMAT', 'V V R', 3*nbm, imatma)
    call wkvect('&&FLUST3.TMP.DMOD', 'V V R', 6*nbm*nbnoe, idefm)
    imatra = imatma+nbm
    imataa = imatra+nbm
    do im = 1, nbm
        ior = nuor(im)
        call rsadpa(base, 'L', 1, 'MASS_GENE', ior, &
                    0, sjv=lmasg, styp=k8b)
        call rsadpa(base, 'L', 1, 'RIGI_GENE', ior, &
                    0, sjv=lrigg, styp=k8b)
        call rsadpa(base, 'L', 1, 'FACT_PARTICI_DX', ior, &
                    0, sjv=lfacx, styp=k8b)
        zr(imatma+im-1) = zr(lmasg)
        zr(imatra+im-1) = zr(lrigg)
        zr(ifpart+im-1) = zr(lfacx)
        zr(ifpart+nbm+im-1) = zr(lfacx+1)
        zr(ifpart+2*nbm+im-1) = zr(lfacx+2)
        zr(imataa+im-1) = 4.d0*pi*zr(imatma+im-1)*amor(im)*zr(ifreqi+ior-1)
    end do
    call extmod(base, numddl, nuor, nbm, zr(idefm), &
                neq, nbnoe, iddl, 6)
!
!
! --- 10.INTERPOLATION DES DEFORMEES MODALES AUX POINTS DE
! ---    DISCRETISATION
!
    n = nbz*nbgrp*nbm
    call wkvect('&&FLUST3.TMP.PHI', 'V V R', 2*n, iphix)
    call wkvect('&&FLUST3.TMP.IZ', 'V V R', nbz, iz)
    call wkvect('&&FLUST3.TMP.NUM', 'V V I', nbz, inum)
    iphiy = iphix+n
    nbddl = 6
!
    call mefint(nbz, nbgrp, nbm, nbnoe, nbddl, &
                irot, zi(iadnog), zi(inbnog), zr(izint), zr(idefm), &
                zr(iphix), zr(iphiy), zr(iz), zi(inum))
!
!
! --- 11.VERIFICATIONS GEOMETRIQUES SI PRESENCE DE GRILLES
!
    if (ntypg .ne. 0) then
        zmin = zr(iz)
        zmax = zr(iz+nbz-1)
        call mefgri(ntypg, nbgtot, zr(izg), zr(ilongg), zi(jtypg), &
                    zmin, zmax)
    end if
!
!
! --- 12.TABLEAU DES DIMENSIONS
!
! --- ORDRE DE TRONCATURE DES SERIES DE LAURENT DANS LA BASE MODALE
    nbtron = 3
    ndim(1) = nbz
    ndim(2) = nbm
    ndim(3) = nbcyl
    ndim(4) = nbgrp
    ndim(5) = nbtron
    ndim(6) = iencei
    ndim(7) = nima
    ndim(8) = nima2
    ndim(9) = npv
    ndim(10) = ndir
    ndim(11) = nbnoe
    ndim(12) = neq
    ndim(13) = ntypg
    ndim(14) = nbgtot
!
!
! --- 13.APPEL DE LA PROCEDURE DE RESOLUTION
    call mefist(melflu, ndim, som, alpha, ru, &
                promas, provis, zr(imatma), zi(igreq), nuor, &
                freq, masg, fact, zr(ifpart), vite, &
                zr(ixint), zr(iyint), zr(irint), zr(iz), zr(iphix), &
                zr(iphiy), zr(idefm), zi(jtypg), zr(izg), zr(ilongg), &
                zr(ilargg), zr(iepaig), zr(icdg), zr(icpg), zr(irugg), &
                base)
!
!
! --- 14.IMPRESSIONS DANS LE FICHIER RESULTAT SI DEMANDEES ---
!
    if (nivpar .eq. 1 .or. nivdef .eq. 1) then
        carac(1) = som(9)
        carac(2) = 0.d0
        calcul(1) = .true.
        calcul(2) = .false.
        call fluimp(3, nivpar, nivdef, melflu, typflu, &
                    nuor, freq, zr(ifreqi), nbm, vite, &
                    npv, carac, calcul, [0.d0])
    end if
!
!
! --- MENAGE
!
    do i = 1, nbgrp
        call codent(i, 'D0', numgno)
        grpno = '&&MEFGMN.'//numgno
        call jedetr(grpno)
    end do
    call jedetr('&&FLUST3.TMP.NBCC')
    call jedetr('&&FLUST3.TMP.REQ')
    call jedetr('&&FLUST3.TMP.NOEQ')
    call jedetr('&&FLUST3.TMP.NOCY')
    call jedetr('&&FLUST3.TMP.CENT')
    call jedetr('&&FLUST3.TMP.TYPG')
    call jedetr('&&FLUST3.TMP.GRIL')
    call jedetr('&&FLUST3.TMP.NGX')
    call jedetr('&&FLUST3.TMP.GREQ')
    call jedetr('&&FLUST3.TMP.ADR')
    call jedetr('&&FLUST3.TMP.XYZR')
    call jedetr('&&FLUST3.LIFREQ')
    call jedetr('&&FLUST3.TMP.IFAC')
    call jedetr('&&FLUST3.TMP.IMAT')
    call jedetr('&&FLUST3.TMP.DMOD')
    call jedetr('&&FLUST3.TMP.PHI')
    call jedetr('&&FLUST3.TMP.IZ')
    call jedetr('&&FLUST3.TMP.NUM')
!
    call jedema()
end subroutine
