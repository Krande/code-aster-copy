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
! person_in_charge: daniele.colombo at ifpen.fr
! aslint: disable=W1306
!
subroutine xddlimf(modele, ino, cnxinv, jnoxfv, motcle, &
                   ch2, ndim, lsn, lst, valimr, valimf, valimc, &
                   fonree, lisrel, nomn, direct, class, mesh, &
                   hea_no)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/afrela.h"
#include "asterfort/cesexi.h"
#include "asterfort/elelin.h"
#include "asterfort/elrfvf.h"
#include "asterfort/focste.h"
#include "asterfort/fointe.h"
#include "asterc/getres.h"
#include "asterfort/gnomsd.h"
#include "asterfort/ismali.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnum.h"
#include "asterfort/reeref.h"
#include "asterfort/xcalc_code.h"
#include "asterfort/xcalc_heav.h"
#include "asterfort/xdvois.h"
#include "asterfort/jedetr.h"
#include "asterfort/int_to_char8.h"
!
    integer(kind=8) :: ino, jnoxfv, ndim
    real(kind=8) :: lsn(4), lst(4), valimr, direct(3)
    complex(kind=8) :: valimc
    character(len=4) :: fonree
    character(len=8) :: modele, motcle, valimf, nomn
    character(len=8), intent(in) :: mesh
    character(len=19) :: cnxinv, ch2, lisrel, hea_no
    aster_logical :: class
!
!      TRAITEMENT DE DDL_IMPO SUR UN NOEUD X-FEM
!             (POUR MOTCLE = DX, DY ,DZ)
!      TRAITEMENT DE DDL_IMPO SUR UN NOEUD HM-XFEM
!             (POUR MOTCLE = DX, DY, DZ ET/OU PRE1)
!      TRAITEMENT DE FACE_IMPO SUR UN NOEUD X-FEM
!             (POUR DNOR OU DTAN : MOTCLE = DEPL )
!
! IN  INO    : NUMERO DU NOEUD
! IN  MESH   : NOM DU MAILLAGE
! IN  LSN    : LSN DU NOEUD COURANT
! IN  NOMN   : NOM DU NOEU COURANT
!
! OUT CLASS  : ON AFFECTE AUSSI LA RELATION CINEMATIQUE "TOTALE"
!
    integer(kind=8) :: nbxcmp, nbnomax
    parameter(nbxcmp=60, nbnomax=20)
    integer(kind=8) :: ier, nbno, jconx2, nbmano, jma, adrma, numa, voisin(3), dimens(nbxcmp)
    integer(kind=8) :: itypma, ibid, nbnoma, nno, i, jlsnd, jlsnl, hea_pt, heavm(135), jheavnl
    integer(kind=8) :: iad, ima, j, nuno, nuno2, iadrco, icode, numac, nbnomac
    integer(kind=8) :: nterm, ncompn, jheavnd
    real(kind=8) :: lsno(3), lsn2, coor(4*ndim), param(1), alpha(1), geom(nbnomax*ndim)
    real(kind=8) :: ff(nbnomax), eps, ptm(ndim), ptp(ndim), xe(3)
    real(kind=8) :: valpar(ndim), deplm, deplp, deplun, deplde, depltr, coef(nbxcmp)
    real(kind=8) :: valh, valc, sign, deplmi, ffb(nbnomax), a1, b1, c1, a2, b2, c2
    character(len=8) :: noma, typma, elp, elpq, arete, nompar(ndim), nomres
    character(len=8) :: name_node, name_ma(20), ddl(nbxcmp), noeud(nbxcmp)
    character(len=1) :: axes(3)
    character(len=19) :: mai, fclas, fenri
    character(len=24) :: coorn, noojb
    character(len=16) :: typres, nomcmd
    parameter(eps=2.d-6)
    integer(kind=8), pointer :: nunotmp(:) => null()
    integer(kind=8), pointer :: connex(:) => null()
    character(len=8), pointer :: lgrf(:) => null()
    integer(kind=8), pointer :: ihea_no(:) => null()
    real(kind=8), pointer :: lsnv(:) => null()
    complex(kind=8) :: cbid
    aster_logical :: milieu, passe, coupee
    cbid = dcmplx(0.d0, 0.d0)
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
    call jeveuo(ch2//'.CESV', 'L', vr=lsnv)
    call jeveuo(ch2//'.CESL', 'L', jlsnl)
    call jeveuo(ch2//'.CESD', 'L', jlsnd)
!
    do i = 1, nbxcmp
        dimens(i) = 0
        noeud(i) = nomn
    end do
    axes(1) = 'X'
    axes(2) = 'Y'
    axes(3) = 'Z'
!
! --- RECUPERATION DE LA LISTE DES NOEUDS AFFECTES PAR LA CONDITION
    call jeexin('&&CADDLI.LIST_NODE', ier)
    if (ier .ne. 0) then
        call jeveuo('&&CADDLI.LIST_NODE', 'L', vi=nunotmp)
        call jelira('&&CADDLI.LIST_NODE', 'LONMAX', nbno)
    end if
    call jeexin('&&CAFACI.LIST_NODE', ier)
    if (ier .ne. 0) then
        call jeveuo('&&CAFACI.LIST_NODE', 'L', vi=nunotmp)
        call jelira('&&CAFACI.LIST_NODE', 'LONMAX', nbno)
    end if
    call jeexin('&&CAAREI.LIST_NODE', ier)
    if (ier .ne. 0) then
        call jeveuo('&&CAAREI.LIST_NODE', 'L', vi=nunotmp)
        call jelira('&&CAAREI.LIST_NODE', 'LONMAX', nbno)
    end if
! --- RECUPERATION DU NOM DU MAILLAGE :
    call jeveuo(modele//'.MODELE    .LGRF', 'L', vk8=lgrf)
    noma = lgrf(1)
! --- RECUPERATION DES MAILLES CONTENANT LE NOEUD
    call jeveuo(noma//'.CONNEX', 'L', vi=connex)
    call jeveuo(jexatr(noma//'.CONNEX', 'LONCUM'), 'L', jconx2)
    call jelira(jexnum(cnxinv, ino), 'LONMAX', nbmano)
    call jeveuo(jexnum(cnxinv, ino), 'L', adrma)
! --- ON RECUPERE LE TYPE DE MAILLE
    mai = noma//'.TYPMAIL'
    call jeveuo(mai, 'L', jma)
    numa = zi(jnoxfv-1+2*(ino-1)+1)
    itypma = zi(jma-1+numa)
    call jenuno(jexnum('&CATA.TM.NOMTM', itypma), typma)
! --- INITIALISATION
!     INDICATEUR POUR IMPOSER LA RELATION CINEMATIQUE "TOTALE" OU NON
    class = .true.
!     INDICATEUR POUR LES NOEUDS QUE L'ON DOIT TRAITER UNIQUEMENT AVEC LA
!     RELATION CINEMATIQUE "TOTALE"
    passe = .false.
    call xdvois(typma, ino, noma, numa, jlsnd, jlsnl, jconx2, &
                ch2, lsn, nbmano, jma, adrma, ndim, coupee, &
                nno, arete, milieu, lsno, voisin)
! --- SI DEUX VOISINS N'ONT PAS ETE TROUVES, CE NOEUD N'EST PAS UN NOEUD XFEM ACTIF
    if (voisin(1)*voisin(2) .eq. 0 .and. .not. milieu) then
        passe = .true.
    end if
! --- ON IMPOSERA TOUJOURS LES CONDITIONS AUX LIMITES TOTALES SAUF POUR LES
!     NOEUDS SITUES SUR LA FISSURE
    if (lsn(1) .eq. 0.d0) then
        class = .false.
    end if
! --- RECUPERATION DES COORDONNEES DES NOEUDS
    if (.not. passe) then
        coorn = noma//'.COORDO    .VALE'
        call jeveuo(coorn, 'L', iadrco)
        do i = 1, ndim
            coor(i) = zr(iadrco-1+3*(ino-1)+i)
            coor(ndim+i) = zr(iadrco-1+3*(voisin(1)-1)+i)
            coor(2*ndim+i) = zr(iadrco-1+3*(voisin(2)-1)+i)
            if (.not. ismali(typma)) then
                coor(3*ndim+i) = zr(iadrco-1+3*(voisin(3)-1)+i)
            end if
        end do
!     L'ARETE EST COUPEE
        if (coupee) then
! --- ON EVALUE ALPHA
            if (lsn(1) .ne. 0.d0) then
                param(1) = 0.d0
                call reeref(arete, nno, lsno, param, 1, alpha, &
                            ff(1:nno))
            end if
! --- ON RECUPERE LES COORDONNES DE PTM ET PTP
            if (lsn(1) .eq. 0.d0) then
                do i = 1, ndim
                    ptm(i) = (1.d0-eps)*coor(i)+eps*coor(i+ndim)
                    ptp(i) = (1.d0-eps)*coor(i)+eps*coor(i+2*ndim)
                end do
            else
                do i = 1, ndim
                    ptm(i) = (1.d0-(1.d0-eps)*(5.d-1+alpha(1)/2.d0))*coor(ndim+i)+(1.d0-eps)* &
                             (alpha(1)/2.d0+5.d-1)*coor(i+2*ndim)
                    ptp(i) = (1.d0-eps)*(5.d-1-alpha(1)/2.d0)*coor(ndim+i)+ &
                             (1.d0-(1.d0-eps)*(-alpha(1)/2.d0+5.d-1))*coor(i+2*ndim)
                end do
            end if
! --- ON EVALUE LA FONCTION EN PTM ET PTP AINSI QU'AUX NOEUDS
            nompar(1) = 'X'
            nompar(2) = 'Y'
            if (ndim .eq. 3) then
                nompar(3) = 'Z'
            end if
            do i = 1, ndim
                valpar(i) = ptm(i)
            end do
            call fointe('FM', valimf, ndim, nompar, valpar, &
                        deplm, icode)
            do i = 1, ndim
                valpar(i) = ptp(i)
            end do
            call fointe('FM', valimf, ndim, nompar, valpar, &
                        deplp, icode)
            do i = 1, ndim
                valpar(i) = coor(ndim+i)
            end do
            call fointe('FM', valimf, ndim, nompar, valpar, &
                        deplun, icode)
            do i = 1, ndim
                valpar(i) = coor(2*ndim+i)
            end do
            call fointe('FM', valimf, ndim, nompar, valpar, &
                        deplde, icode)
            if (.not. ismali(typma) .and. lsn(1) .ne. 0.d0 .and. lsno(3) .ne. 0.d0) then
                do i = 1, ndim
                    valpar(i) = coor(3*ndim+i)
                end do
                call fointe('FM', valimf, ndim, nompar, valpar, &
                            depltr, icode)
            end if
! --- PREMIER CAS, LE NOEUD EST SUR LA FISSURE
            if (lsn(1) .eq. 0.d0 .and. lst(1) .lt. 0.d0) then
                valh = (deplp-deplm)/2.d0
                valc = deplp
            end if
! --- DEUXIEME CAS: NOEUD HORS FISSURE, MAILLAGE LINEAIRE
!     ON EVALUE LA FONCTION AUX NOEUDS VOISINS 1 ET 2
            if ((ismali(typma) .and. lsn(1) .ne. 0.d0)) then
                if (lsn(1) .lt. 0.d0) then
                    valh = (deplp-ff(2)*deplde)/(2.d0*ff(1))-deplun/2.d0
                elseif (lsn(1) .gt. 0.d0) then
                    valh = (ff(1)*deplun-deplm)/(2.d0*ff(2))+deplde/2.d0
                else
                    ASSERT(.false.)
                end if
            end if
! --- TROISIEME CAS: MAILLLAGE QUADRATIQUE ET NOEUD HORS FISSURE, ET NOEUD
!     MILIEU SUR LA FISSURE
            if (.not. ismali(typma) .and. lsn(1) .ne. 0.d0 .and. lsno(3) .eq. 0.d0) then
                param(1) = -sign(1.d0, lsn(1))*5.d-1
!     ON EVALUE LA FONCTION EN PARAM
                call elrfvf(arete, param, ff)
                do i = 1, ndim
                    valpar(i) = coor(ndim+i)*ff(1)+coor(2*ndim+i)*ff(2)+coor(3*ndim+i)*ff(3)
                end do
                call fointe('FM', valimf, ndim, nompar, valpar, &
                            deplmi, icode)
                if (lsn(1) .gt. 0.d0) then
                    valh = (ff(1)*deplun+ff(3)*deplm-deplmi)/(2.d0*ff(2))+deplde/2.d0
                else
                    valh = (deplmi-ff(2)*deplde-ff(3)*deplp)/(2.d0*ff(1))-deplun/2.d0
                end if
            end if
! --- QUATRIEME CAS: MAILLLAGE QUADRATIQUE ET NOEUD HORS FISSURE, ET NOEUD
!     MILIEU PAS SUR LA FISSURE
            if (.not. ismali(typma) .and. lsn(1) .ne. 0.d0 .and. lsno(3) .ne. 0.d0) then
                call elrfvf(arete, alpha, ff)
                param(1) = (sign(1.d0, alpha(1))*1.d0+alpha(1))/2.d0
!     ON EVALUE LA FONCTION EN PARAM
                call elrfvf(arete, param, ffb)
                do i = 1, ndim
                    valpar(i) = coor(ndim+i)*ffb(1)+coor(2*ndim+i)*ffb(2)+coor(3*ndim+i)*ffb(3)
                end do
                call fointe('FM', valimf, ndim, nompar, valpar, &
                            deplmi, icode)
                if (alpha(1) .gt. 0.d0) then
                    a1 = 2.d0*ff(1)
                    b1 = 2.d0*ff(3)
                    c1 = deplp-deplun*ff(1)-deplde*ff(2)-depltr*ff(3)
                    a2 = 2.d0*ffb(1)
                    b2 = 2.d0*ffb(3)
                    c2 = deplmi-deplun*ffb(1)-deplde*ffb(2)-depltr*ffb(3)
                else
                    a1 = -2.d0*ff(2)
                    b1 = -2.d0*ff(3)
                    c1 = deplm-deplun*ff(1)-deplde*ff(2)-depltr*ff(3)
                    a2 = -2.d0*ffb(2)
                    b2 = -2.d0*ffb(3)
                    c2 = deplmi-deplun*ffb(1)-deplde*ffb(2)-depltr*ffb(3)
                end if
                ASSERT(a1*b2 .ne. a2*b1)
!     SI NOEUD MILIEU
                if (milieu) then
                    valh = (a1*c2-c1*a2)/(a1*b2-a2*b1)
!     SI NOEUD SOMMET LE PLUS PROCHE DE LA FISSURE
                elseif (alpha(1)*lsn(1) .ge. 0.d0) then
                    if (lsn(1) .gt. 0.d0) then
                        valh = (deplun*ff(1)+depltr*ff(3)-deplm)/(2.d0*ff(2))+deplde/2.d0
                    else
                        valh = (deplp-ff(3)*depltr-deplde*ff(2))/(2.d0*ff(1))-deplun/2.d0
                    end if
!     SI NOEUD SOMMET LE LOIN PROCHE DE LA FISSURE
                else
                    valh = (c1*b2-b1*c2)/(a1*b2-a2*b1)
                end if
            end if
! --- CREATION DES FONCTIONS CONSTANTES POUR LA VALEUR A IMPOSER SUR LES
!     DDLS CLASSIQUES ET ENRICHIS
            call getres(nomres, typres, nomcmd)
            noojb = '12345678.FXF.123456.PROL'
            if (lsn(1) .eq. 0.d0) then
                call gnomsd(nomres, noojb, 14, 19)
                fclas = noojb(1:19)
                call focste(fclas, 'XXX', valc, 'G')
            end if
            call gnomsd(nomres, noojb, 14, 19)
            fenri = noojb(1:19)
            call focste(fenri, 'XXX', valh, 'G')
        elseif (.not. coupee) then
            passe = .true.
!     ON BOUCLE SUR LES MAILLES AUXQUELLES LE NOEUD APPARTIENT
            lsn2 = 0.d0
            nbnomac = 0
            numac = 0
            do ima = 1, nbmano
                numa = zi(adrma-1+ima)
                itypma = zi(jma-1+numa)
                call jenuno(jexnum('&CATA.TM.NOMTM', itypma), typma)
                if (typma .eq. 'SEG3') then
                    goto 80
                elseif (typma .eq. 'SEG2') then
                    goto 80
                elseif (typma .eq. 'TRIA6') then
                    if (ndim .eq. 3) goto 80
                    elp = 'TR6'
                elseif (typma .eq. 'TRIA3') then
                    if (ndim .eq. 3) goto 80
                    elp = 'TR3'
                elseif (typma .eq. 'QUAD8') then
                    if (ndim .eq. 3) goto 80
                    elp = 'QU8'
                elseif (typma .eq. 'QUAD4') then
                    if (ndim .eq. 3) goto 80
                    elp = 'QU4'
                elseif (typma .eq. 'TETRA4') then
                    elp = 'TE4'
                elseif (typma .eq. 'TETRA10') then
                    elp = 'T10'
                elseif (typma .eq. 'PYRAM5') then
                    elp = 'PY5'
                elseif (typma .eq. 'PYRAM13') then
                    elp = 'P13'
                elseif (typma .eq. 'PENTA6') then
                    elp = 'PE6'
                elseif (typma .eq. 'PENTA15') then
                    elp = 'P15'
                elseif (typma .eq. 'HEXA8') then
                    elp = 'HE8'
                elseif (typma .eq. 'HEXA20') then
                    elp = 'H20'
                end if
!     SI LE DDL EST PRE1, ON PREND LA MAILLE LINEAIRE
                if (motcle .eq. 'PRE1') then
                    elpq = elp
                    call elelin(3, elpq, elp, ibid, nbnoma)
                end if
!     ON CHERCHE UN NOEUD DE LA MAILLE TEL QUE LSN SOIT DE SIGNE OPPOSE
                ptp(:) = 0.d0
                call elrfvf(elp, ptp, ff, nbnoma)
                coorn = noma//'.COORDO    .VALE'
                call jeveuo(coorn, 'L', iadrco)
                do i = 1, nbnoma
                    nuno = connex(zi(jconx2+numa-1)+i-1)
! ---         ON REGARDE SI LE NOEUD APPARTIENT AU GRP DE NOEUD AFFECTE
                    do j = 1, nbno
                        nuno2 = nunotmp(j)
                        if (nuno2 .eq. nuno) then
                            call cesexi('C', jlsnd, jlsnl, numa, i, &
                                        1, 1, iad)
! ---         ON CHERCHE LE NOEUD DE SIGNE LSN OPPOSE ET SITUE LE PLUS LOIN
!             POSSIBLE DE LA FISSURE
                            if (lsnv(iad)*lsn(1) .le. 0.d0) then
                                if (abs(lsnv(iad)) .gt. abs(lsn2)) then
                                    voisin(1) = connex(zi(jconx2+numa-1)+i-1)
                                    lsn2 = lsnv(iad)
                                    numac = numa
                                    nbnomac = nbnoma
                                end if
                            end if
                        end if
                    end do
                end do
80              continue
            end do
            if (numac .eq. 0) goto 78
!     ON RECUPERE LES COORDONNEES DE LA MAILLE SELECTIONNEE
            heavm(1:135) = 0
            call jeveuo(hea_no//'.CESV', 'L', vi=ihea_no)
            call jeveuo(hea_no//'.CESL', 'L', jheavnl)
            call jeveuo(hea_no//'.CESD', 'L', jheavnd)
            ncompn = zi(jheavnd-1+5+4*(numac-1)+3)
            do i = 1, nbnomac
                do j = 1, ndim
                    geom(ndim*(i-1)+j) = zr(iadrco-1+3*(connex(zi(jconx2+numac-1)+i-1)-1)+j)
                end do
                name_node = int_to_char8(connex(zi(jconx2+numac-1)+i-1))
                name_ma(i) = name_node
!     RECUPERATION DE LA DEFINITION DES DDLS HEAVISIDES
!            ASSERT(ncompn.eq.5)
                do j = 1, ncompn
                    call cesexi('S', jheavnd, jheavnl, numac, i, &
                                1, j, iad)
                    heavm(ncompn*(i-1)+j) = ihea_no(iad)
                end do
            end do
            hea_pt = xcalc_code(1, he_real=[-sign(1.d0, lsn(1))])
!     ON CHERCHE LES COORDONNEES DE CE NOEUD
            do i = 1, ndim
                coor(i) = zr(iadrco-1+3*(ino-1)+i)
                coor(ndim+i) = zr(iadrco-1+3*(voisin(1)-1)+i)
            end do
!     ON CHERCHE LES COORDONNES DU POINT OU L'ON VA EVALUER LA FONCTION
            nompar(1) = 'X'
            nompar(2) = 'Y'
            if (ndim .eq. 3) then
                nompar(3) = 'Z'
            end if
            do i = 1, ndim
                valpar(i) = ((abs(lsn(1))+abs(lsn2)/2.d0)*coor(ndim+i)+ &
                             abs(lsn2)/2.d0*coor(i))/(abs(lsn(1))+abs(lsn2))
            end do
            call fointe('FM', valimf, ndim, nompar, valpar, &
                        deplmi, icode)
            call getres(nomres, typres, nomcmd)
            noojb = '12345678.FXFEM.1234.PROL'
            call gnomsd(nomres, noojb, 16, 19)
            fenri = noojb(1:19)
            call focste(fenri, 'XXX', deplmi, 'G')
            call reeref(elp, nbnomac, geom, valpar, ndim, xe, ff(1:nbnomac))
            if (motcle .eq. 'DX' .or. motcle .eq. 'DY' .or. motcle .eq. 'DZ') then
                do i = 1, nbnomac
                    ddl(2*i-1) = 'D'//motcle(2:2)
                    ddl(2*i) = 'H1'//motcle(2:2)
                    coef(2*i-1) = ff(i)
                    coef(2*i) = xcalc_heav(heavm(1+ncompn*(i-1)), hea_pt, &
                                           heavm(5+ncompn*(i-1)))*ff(i)
                    noeud(2*i-1) = name_ma(i)
                    noeud(2*i) = name_ma(i)
                end do
                nterm = 2*nbnomac
                call afrela(coef, [cbid], ddl, noeud, dimens, &
                            [0.d0], nterm, valimr, valimc, fenri, &
                            'REEL', fonree, 0.d0, lisrel)
            elseif (motcle(1:8) .eq. 'DEPL    ') then
                do j = 1, ndim
                    do i = 1, nbnomac
                        ddl(2*nbnoma*(j-1)+2*i-1) = 'D'//axes(j)
                        coef(2*nbnoma*(j-1)+2*i-1) = direct(j)*ff(i)
                        noeud(2*nbnoma*(j-1)+2*i-1) = name_ma(i)
                        ddl(2*nbnoma*(j-1)+2*i) = 'H1'//axes(j)
                        coef(2*nbnoma*(j-1)+2*i) = direct(j)*ff(i)* &
                                                   xcalc_heav(heavm(1+ncompn*(i-1)), &
                                                              hea_pt, heavm(5+ncompn*(i-1)))
                        noeud(2*nbnoma*(j-1)+2*i) = name_ma(i)
                    end do
                end do
                nterm = nbnomac*2*ndim
                call afrela(coef, [cbid], ddl, noeud, dimens, &
                            [0.d0], nterm, deplmi, valimc, valimf, &
                            'REEL', 'REEL', 0.d0, lisrel)
            elseif (motcle .eq. 'PRE1') then
                do i = 1, nbnomac
                    ddl(2*i-1) = 'PRE1'
                    ddl(2*i) = 'H1PRE1'
                    coef(2*i-1) = ff(i)
                    coef(2*i) = xcalc_heav(heavm(1+ncompn*(i-1)), hea_pt, &
                                           heavm(5+ncompn*(i-1)))*ff(i)
                    noeud(2*i-1) = name_ma(i)
                    noeud(2*i) = name_ma(i)
                end do
                nterm = 2*nbnomac
                call afrela(coef, [cbid], ddl, noeud, dimens, &
                            [0.d0], nterm, deplmi, valimc, valimf, &
                            'REEL', 'REEL', 0.d0, lisrel)
            end if
        end if
    end if
!       CAS FACE_IMPO DNOR OU DTAN
78  continue
    if (motcle(1:8) .eq. 'DEPL    ' .and. .not. passe) then
        if (lsn(1) .eq. 0.d0) then
            i = 0
            do j = 1, ndim
!          COEFFICIENTS ET DDLS DE LA RELATION
                i = i+1
                ddl(i) = 'D'//axes(j)
                coef(i) = direct(j)
            end do
            nterm = i
            call afrela(coef, [cbid], ddl, noeud, dimens, &
                        [0.d0], nterm, valimr, valimc, fclas, &
                        'REEL', fonree, 0.d0, lisrel)
        end if
        i = 0
        do j = 1, ndim
            i = i+1
            ddl(i) = 'H1'//axes(j)
            coef(i) = direct(j)
        end do
        nterm = i
        call afrela(coef, [cbid], ddl, noeud, dimens, &
                    [0.d0], nterm, valimr, valimc, fenri, &
                    'REEL', fonree, 0.d0, lisrel)
!      CAS DDL_IMPO DX DY DZ (ET/OU PRE1 => POUR HM-XFEM ONLY)
    elseif ((motcle .eq. 'DX' .or. motcle .eq. 'DY' .or. motcle .eq. 'DZ') .and. &
            .not. passe) then
!        COEFFICIENTS ET DDLS DE LA RELATION
        if (lsn(1) .eq. 0.d0) then
            ddl(1) = 'D'//motcle(2:2)
            coef(1) = 1.d0
            nterm = 1
            call afrela(coef, [cbid], ddl, noeud, dimens, &
                        [0.d0], nterm, valimr, valimc, fclas, &
                        'REEL', fonree, 0.d0, lisrel)
        end if
        nterm = 1
        ddl(1) = 'H1'//motcle(2:2)
        coef(1) = 1.d0
        call afrela(coef, [cbid], ddl, noeud, dimens, &
                    [0.d0], nterm, valimr, valimc, fenri, &
                    'REEL', fonree, 0.d0, lisrel)
    elseif (motcle .eq. 'PRE1' .and. .not. passe) then
!          COEFFICIENTS ET DDLS DE LA RELATION
        if (lsn(1) .eq. 0.d0) then
            ddl(1) = 'PRE1'
            coef(1) = 1.d0
            nterm = 1
            call afrela(coef, [cbid], ddl, noeud, dimens, &
                        [0.d0], nterm, valimr, valimc, fclas, &
                        'REEL', fonree, 0.d0, lisrel)
        end if
        nterm = 1
        ddl(1) = 'H1'//motcle(1:4)
        coef(1) = 1.d0
        call afrela(coef, [cbid], ddl, noeud, dimens, &
                    [0.d0], nterm, valimr, valimc, fenri, &
                    'REEL', fonree, 0.d0, lisrel)
    end if
!
!    call jedetr(fclas//'.PROL')
!    call jedetr(fclas//'.VALE')
!    call jedetr(fenri//'.PROL')
!    call jedetr(fenri//'.VALE')
!
    call jedema()
end subroutine
