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

subroutine op0184()
    implicit none
!
!     LECTURE D'UN RESULTAT PLEXUS (PRESSION) SUR FICHIER IDEAS
!     LA PRESSION EST CALCULEE PAR PLEXUS SUR DES SEG2 (CONSTANTE)
!     LE MAILLAGE ASTER PEUT COMPORTER DES SEG2, DES SEG3, DES COQUES
!
!     -----------------------------------------------------------------
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterfort/calcul.h"
#include "asterfort/copisd.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/infmaj.h"
#include "asterfort/iunifi.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/nbelem.h"
#include "asterfort/nbgrel.h"
#include "asterfort/pj3da4.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsagsd.h"
#include "asterfort/rscrsd.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsinfo.h"
#include "asterfort/rsnoch.h"
#include "asterfort/rsorac.h"
#include "asterfort/titre.h"
#include "asterfort/typele.h"
#include "asterfort/ulisop.h"
#include "asterfort/ulopen.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
!
    integer(kind=8) :: nbv, nbmapl, ibid, ntout, nnume, np, iul
    integer(kind=8) :: vali(2)
    integer(kind=8) :: nbordr, jnume, n1, nlinst, jlist, nbinst, nis, nc, nfor
    integer(kind=8) :: iret, nbtitr, ifsig, l, ipas, k, numpas, ino
    integer(kind=8) :: ivar(6), iord, jpres, ima, nbnoas, jnoma, i, nnu
    integer(kind=8) :: iadrno, imp, jdme, ntseg, imamin, jdco, jdno, no1, no2, nutyel
    integer(kind=8) :: jcelv, idec, nbelgr, liel, iel, iadno, nno
    integer(kind=8) :: jinst, nbtrou, lordr, ntpoi, itest, iad, imapl
    integer(kind=8) :: nbordt, te, nbgr, igr, tord(1)
    real(kind=8) :: rbid, pres, epsi, temps, tref, cm(3), a(3), b(3), la, lb, d2
    real(kind=8) :: d2min
    real(kind=8) :: valr
    complex(kind=8) :: cbid
    character(len=6) :: kar
    character(len=8) :: resu, nomapl, nomast, k8b, listr8, form, crit
    character(len=8) :: nomo, lpain(1), lpaout(1)
    character(len=16) :: nomcmd, concep, nsymb, nomte, k16nom
    character(len=19) :: nomch, ligrmo, chpres, capres
    character(len=24) :: coorn, typma, coorp, mlgcnx, lchin(1), lchout(1)
    character(len=24) :: noliel, chgeom, option, connex
    character(len=80) :: k80b, k80bm, k80bid
    integer(kind=8), pointer :: celd(:) => null()
    character(len=80), pointer :: titr(:) => null()
!     -----------------------------------------------------------------
!
    call jemarq()
    call infmaj()
    k80bm = ' '
    call getres(resu, concep, nomcmd)
!
!     LECTURE DE LA NUMEROTATION DES DDL
!
    call getvid(' ', 'MAIL_PLEXUS', scal=nomapl, nbret=nbv)
    call dismoi('NB_MA_MAILLA', nomapl, 'MAILLAGE', repi=nbmapl)
    call getvid(' ', 'MAILLAGE', scal=nomast, nbret=nbv)
    call dismoi('NB_NO_MAILLA', nomast, 'MAILLAGE', repi=nbnoas)
    call wkvect('&&OP0184.NOAST_MAPLEX', 'V V I', nbnoas, jnoma)
    coorn = nomast//'.COORDO    .VALE'
    call jeveuo(coorn, 'L', iadrno)
    coorp = nomapl//'.COORDO    .VALE'
    call jeveuo(coorp, 'L', jdco)
    typma = nomapl//'.TYPMAIL'
    call jeveuo(typma, 'L', jdme)
    mlgcnx = nomapl//'.CONNEX'
    call jenonu(jexnom('&CATA.TM.NOMTM', 'SEG2'), ntseg)
!
! IL PEUT Y AVOIR DES POI1 DANS LES ELEMENTS PLEXUS ON LES IGNORE
!
    call jenonu(jexnom('&CATA.TM.NOMTM', 'POI1'), ntpoi)
    itest = 0
!
    do ino = 1, nbnoas
        do i = 1, 3
            cm(i) = zr(iadrno+3*ino-3+i-1)
        end do
        imamin = 0
        d2min = 1.d10
        do imp = 1, nbmapl
            nutyel = zi(jdme+imp-1)
            if (nutyel .eq. ntpoi) then
                itest = -1
                goto 40
            else if (nutyel .eq. ntseg) then
                goto 20
            else
                call utmess('F', 'UTILITAI3_9')
            end if
20          continue
            call jeveuo(jexnum(mlgcnx, imp), 'L', jdno)
            no1 = zi(jdno)
            no2 = zi(jdno+1)
            do i = 1, 3
                a(i) = zr(jdco+(no1-1)*3+i-1)
                b(i) = zr(jdco+(no2-1)*3+i-1)
            end do
            call pj3da4(cm, a, b, la, lb, &
                        d2)
            if (d2 .lt. d2min) then
                imamin = imp
                d2min = d2
            end if
40          continue
        end do
        zi(jnoma-1+ino) = imamin
    end do
! TEST SUR LA PRESENCE DE MAILLE PONCTUELLE
    if (itest .ne. 0) then
        call utmess('I', 'UTILITAI3_10')
    end if
!
!     --- QUELS SONT LES INSTANTS A RELIRE ---
!
    nbordr = 0
    call getvtx(' ', 'TOUT_ORDRE', scal=k8b, nbret=ntout)
    if (ntout .eq. 0) then
        call getvis(' ', 'NUME_ORDRE', nbval=0, nbret=nnume)
        if (nnume .ne. 0) then
            nbordr = -nnume
            call wkvect('&&OP0184.NUME_ORDRE', 'V V I', nbordr, jnume)
            call getvis(' ', 'NUME_ORDRE', nbval=nbordr, vect=zi(jnume), nbret=n1)
        else
            call getvid(' ', 'LIST_ORDRE', scal=listr8, nbret=nnu)
            if (nnu .ne. 0) then
                call jeveuo(listr8//'.VALE', 'L', jnume)
                call jelira(listr8//'.VALE', 'LONMAX', nbordr)
            else
                call getvid(' ', 'LIST_INST', scal=listr8, nbret=nlinst)
                if (nlinst .ne. 0) then
                    call jeveuo(listr8//'.VALE', 'L', jlist)
                    call jelira(listr8//'.VALE', 'LONMAX', nbordr)
                else
                    call getvr8(' ', 'INST', nbval=0, nbret=nis)
                    if (nis .ne. 0) then
                        nbinst = -nis
                        call wkvect('&&OP0184.INST', 'V V R', nbordr, jlist)
                        call getvr8(' ', 'INST', nbval=nbordr, vect=zr(jlist), nbret=n1)
                    end if
                end if
            end if
        end if
    end if
!
!     --- LECTURE DE LA PRECISION ET DU CRITERE ---
!
    call getvr8(' ', 'PRECISION', scal=epsi, nbret=np)
    call getvtx(' ', 'CRITERE', scal=crit, nbret=nc)
!
!     FORMAT IDEAS OBLIGATOIRE
!
    call getvtx(' ', 'FORMAT', scal=form, nbret=nfor)
!
    if (form .ne. 'IDEAS') then
        call utmess('F', 'UTILITAI3_11')
    end if
    call jeexin(nomapl//'           .TITR', iret)
    if (iret .eq. 0) then
        call utmess('F', 'UTILITAI3_12')
    else
        call jeveuo(nomapl//'           .TITR', 'L', vk80=titr)
        call jelira(nomapl//'           .TITR', 'LONMAX', nbtitr)
        if (nbtitr .ge. 1) then
            if (titr(1) (10:31) .ne. 'AUTEUR=INTERFACE_IDEAS') then
                call utmess('F', 'UTILITAI3_12')
            end if
        else
            call utmess('A', 'UTILITAI3_13')
        end if
    end if
!
!     CREATION DE LA SD RESULTAT
!
    nbordt = max(1, nbordr)
    call rscrsd('G', resu, 'EVOL_CHAR', nbordt)
!
!     CREATION DU CHAMP DE PRESSION
!
    call getvid(' ', 'MODELE', scal=nomo, nbret=nbv)
    chgeom = nomast//'.COORDO'
    lpain(1) = 'PGEOMER'
    lchin(1) = chgeom
    lpaout(1) = 'PPRES_R'
    chpres = '&&OP0184.CHPRES'
    lchout(1) = chpres
    ligrmo = nomo//'.MODELE'
    option = 'TOU_INI_ELNO'
    call calcul('S', option, ligrmo, 1, lchin, &
                lpain, 1, lchout, lpaout, 'V', &
                'OUI')
    call jeveuo(chpres//'.CELD', 'L', vi=celd)
    call jeveuo(chpres//'.CELV', 'E', jcelv)
    connex = nomast//'.CONNEX'
    noliel = ligrmo//'.LIEL'
    nbgr = nbgrel(ligrmo)
!
    capres = '&&OP0184.PRES'
    call wkvect(capres, 'V V R', nbmapl, jpres)
!
    call getvis(' ', 'UNITE', scal=ifsig, nbret=l)
    k16nom = ' '
    if (ulisop(ifsig, k16nom) .eq. 0) then
        call ulopen(ifsig, ' ', ' ', 'NEW', 'O')
    end if
    ipas = 0
!
!     LECTURE DES DATASET DU FICHIER IDEAS
!     ON NE LIT QUE DES CHAMPS CONSTANTS PAR ELEMENTS
!
60  continue
!
    read (ifsig, '(A6)', end=160, err=180) kar
    if (kar .eq. '    56') then
        read (ifsig, '(40A2)', end=180) k80b
        read (ifsig, '(40A2)', end=180) k80b
        read (ifsig, '(40A2)', end=180) k80b
        read (ifsig, '(A80)', end=180) k80b
! LECTURE DE LA VARIABLE PRESSION EN FONCTION DU TYPE
! DE MATERIAU PLEXUS UTILISE
!
!
        if (.not. ( &
            k80b(49:52) .eq. 'MULT' .or. k80b(49:52) .eq. 'EAU ' .or. k80b(49:52) .eq. 'FLUI')) &
            goto 60
        read (ifsig, '(40A2)', end=180) k80bid
        k8b = k80bid(1:8)
    else if (kar .eq. '  2414') then
        read (ifsig, '(1I10)', end=180) ibid
        read (ifsig, '(40A2)', end=180) k80b
        read (ifsig, '(1I10)', end=180) ibid
        if (ibid .ne. 2) goto 60
        read (ifsig, '(40A2)', end=180) k80b
        read (ifsig, '(40A2)', end=180) k80b
        read (ifsig, '(40A2)', end=180) k80b
        read (ifsig, '(A80)', end=180) k80bm
        read (ifsig, '(40A2)', end=180) k80b
    else
        goto 60
    end if
    read (ifsig, '(6I10)', end=180) (ivar(k), k=1, 6)
!
!        IVAR(3) : TYPE DE DONNEE =0 POUR UNKNOWN
!        IVAR(5) : TYPE DE DONNEE =2 POUR REELLE
!        IVAR(6) : NOMBRE DE VALEURS PAR ELEMENT =20 ICI POUR MULT
!
!
    if (k80b(49:52) .eq. 'MULT') then
        if (ivar(3) .ne. 0) goto 60
        if (ivar(5) .ne. 2) goto 60
        if (ivar(6) .ne. 20) goto 60
    else if (k80b(49:52) .eq. 'FLUI') then
        if (ivar(3) .ne. 0) goto 60
        if (ivar(5) .ne. 2) goto 60
        if (ivar(6) .ne. 2) goto 60
    else if (k80bm(49:52) .eq. 'FLUI') then
        if (ivar(3) .ne. 0) goto 60
        if (ivar(5) .ne. 2) goto 60
        if (ivar(6) .ne. 2) goto 60
    else
        if (ivar(3) .ne. 0) goto 60
        if (ivar(5) .ne. 2) goto 60
        if (ivar(6) .ne. 10) goto 60
    end if
!
!
!        VERIFICATION QUE LE DATASET EST AU BON INSTANT
!
    if (kar .eq. '    56') then
        read (ifsig, '(4I10)', end=180) ibid, ibid, ibid, numpas
        read (ifsig, '(E13.5)', end=180) temps
    else
        read (ifsig, '(8I10)', end=180) ibid
        read (ifsig, '(8I10)', end=180) ibid
        read (ifsig, '(6E13.5)', end=180) temps
        read (ifsig, '(6E13.5)', end=180) rbid
    end if
!
    if (ntout .ne. 0) then
        goto 90
    else
        if (nbordr .ne. 0) then
            if (kar .eq. '  2414') then
                call utmess('F', 'UTILITAI3_14')
            end if
            do iord = 1, nbordr
                if (zi(jnume+iord-1) .eq. numpas) goto 90
            end do
        else if (nbinst .ne. 0) then
            do iord = 1, nbinst
                tref = zr(jlist+iord-1)
                if (crit(1:4) .eq. 'RELA') then
                    if (abs(tref-temps) .le. abs(epsi*temps)) goto 90
                else if (crit(1:4) .eq. 'ABSO') then
                    if (abs(tref-temps) .le. abs(epsi)) goto 90
                end if
            end do
        end if
        goto 60
    end if
90  continue
    ipas = ipas+1
!
!        LECTURE DES PRESSIONS
!
100 continue
    read (ifsig, '(I10)', end=180) ima
    if (ima .eq. -1) goto 110
    read (ifsig, '(E13.5)', end=180) pres
    zr(jpres-1+ima) = pres
    if (k80bm(49:52) .eq. 'FLUI') then
        goto 100
    else if (k80b(49:52) .eq. 'MULT') then
        read (ifsig, '(E13.5)', end=180) rbid
        read (ifsig, '(E13.5)', end=180) rbid
        read (ifsig, '(E13.5)', end=180) rbid
    else
        read (ifsig, '(E13.5)', end=180) rbid
    end if
!
    goto 100
110 continue
!
    do igr = 1, nbgr
        idec = celd(celd(4+igr)+8)
        te = typele(ligrmo, igr)
        call jenuno(jexnum('&CATA.TE.NOMTE', te), nomte)
        if (nomte .eq. 'MEDKTR3' .or. nomte .eq. 'MEDKQU4' .or. nomte .eq. 'MET3SEG3' .or. &
            nomte .eq. 'MEC3TR7H' .or. nomte .eq. 'MEC3QU9H') then
            nbelgr = nbelem(ligrmo, igr)
            call jeveuo(jexnum(noliel, igr), 'L', liel)
            do iel = 1, nbelgr
                ima = zi(liel-1+iel)
                call jeveuo(jexnum(connex, ima), 'L', iadno)
                call jelira(jexnum(connex, ima), 'LONMAX', nno)
                iad = jcelv-1+idec-1+nno*(iel-1)
                do i = 1, nno
                    ino = zi(iadno-1+i)
                    imapl = zi(jnoma+ino-1)
                    pres = zr(jpres-1+imapl)
!  SUITE AUX CORRECTIONS SUR LE SIGNE DE LA PRESSION
!  ON NE MODIFIE SURTOUT PAS LA PRESSION LUE
!  ==> LES TE SAVENT CE QU4ILS ONT A FAIRE
                    zr(iad+i) = pres
                end do
            end do
        end if
    end do
!
    nsymb = 'PRES'
    call rsexch(' ', resu, nsymb, ipas, nomch, &
                iret)
    if (iret .eq. 100) then
    else if (iret .eq. 110) then
        call rsagsd(resu, 0)
        call rsexch(' ', resu, nsymb, ipas, nomch, &
                    iret)
    else
        vali(1) = ipas
        vali(2) = iret
        call utmess('F', 'UTILITAI8_7', ni=2, vali=vali)
    end if
    call copisd('CHAMP_GD', 'G', chpres, nomch)
    call rsnoch(resu, nsymb, ipas)
    call rsadpa(resu, 'E', 1, 'INST', ipas, &
                0, sjv=jinst, styp=k8b)
    zr(jinst) = temps
    do i = 1, nbmapl
        zr(jpres-1+i) = 0.d0
    end do
!
    goto 60
!
160 continue
!
    call utmess('I', 'UTILITAI8_8')
    call rsorac(resu, 'LONUTI', ibid, rbid, k8b, &
                cbid, epsi, crit, tord, 1, &
                nbtrou)
    nbordr = tord(1)
    if (nbordr .le. 0) then
        call utmess('F', 'UTILITAI2_97')
    end if
    call wkvect('&&OP0184.NUME_ORDR', 'V V I', nbordr, lordr)
    call rsorac(resu, 'TOUT_ORDRE', ibid, rbid, k8b, &
                cbid, epsi, crit, zi(lordr), nbordr, &
                nbtrou)
    do iord = 1, nbordr
        call rsadpa(resu, 'L', 1, 'INST', zi(lordr+iord-1), &
                    0, sjv=jinst, styp=k8b)
        vali(1) = zi(lordr+iord-1)
        valr = zr(jinst)
        call utmess('I', 'UTILITAI8_9', si=vali(1), sr=valr)
    end do
    call utmess('I', 'VIDE_1')
!
    call titre()
!
    iul = iunifi('RESULTAT')
    call rsinfo(resu, iul)
!
    goto 190
!
180 continue
    call utmess('F', 'UTILITAI3_15')
!
190 continue
    call jedema()
!
end subroutine
