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

subroutine crmema(promes, iampee)
!
! CALCUL DE LA PSEUDO-MATRICE MASSE POUR LA MODIFICATION STRUCTURALE
!
!   IN  : PROMES : NOM DU CONCEPT PROJ_MESU_MODAL (POUR BASE DE PROJ)
!   IN  : IAMPEE : ADRESSE DE LA MATRICE MACR_ELEM MASSE
!
    implicit none
!
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterc/r8prem.h"
#include "asterfort/cnocns.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/remome.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsexch.h"
#include "asterfort/rslsvd.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    character(len=8) :: promes
    integer(kind=8) :: iampee
!
!
!
!
    character(len=8) :: basemo, nomres, k8bid, typ, trav, modmes, sol, modlms
    character(len=8) :: nomgd, vals, u, v, wks, modele, ma
    character(len=16) :: nomchp, noresu, maelm, k16bid
    character(len=19) :: chamno, ch1s, nu
    character(len=24) :: vnoeud, vrange, basepr, noeums, baseit, vsu, mesint
    character(len=24) :: modid, vref, refms
!
    integer(kind=8) :: nbmesu, nbvecb, nbord, isol, affici(2)
    integer(kind=8) ::  lred, lrange, lint, ier, iposi, ipuls
    integer(kind=8) :: imod, jmod, iret, llncmp, iddl, lmesu, jddl, iexist
    integer(kind=8) :: iposj, ino, nddle, nddli, ico, ipos
    integer(kind=8) :: lnoeud, ltrav, lredi, lwks, lrefms, lref
    integer(kind=8) ::  jcnsc
    integer(kind=8) :: ibid, nbcmpi, numgd, lmaelm
    integer(kind=8) :: lu, lvals, lv, lvsu
!
    real(kind=8) :: masg, eps
    character(len=8), pointer :: cnsk(:) => null()
    integer(kind=8), pointer :: desm(:) => null()
    integer(kind=8), pointer :: ordr(:) => null()
    integer(kind=8), pointer :: deeq(:) => null()
    real(kind=8), pointer :: cnsv(:) => null()
    integer(kind=8), pointer :: cnsd(:) => null()
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
    call getres(nomres, k8bid, k8bid)
!
! RECUPERATION DES ELEMENTS RELATIFS A L'EXPANSION
    noeums = promes//'.PROJM    .PJMNO'
    call jelira(noeums, 'LONUTI', nbmesu)
!
    refms = promes//'.PROJM    .PJMRF'
    call jeveuo(refms, 'L', lrefms)
    k16bid = zk16(lrefms-1+1)
    modlms = k16bid(1:8)
    nomchp = zk16(lrefms-1+2)
    k16bid = zk16(lrefms-1+3)
    basemo = k16bid(1:8)
!
    basepr = promes//'.PROJM    .PJMBP'
    call jeveuo(basepr, 'L', lred)
    call jelira(basepr, 'LONUTI', nbvecb)
    nbvecb = nbvecb/nbmesu
! NBVECB : NOMBRE DE VECTEURS DE BASE
!
! RECUPERATION DES ELEMENTS RELATIFS AU MACRO ELEMENT
    noresu = nomres
    nu = noresu(1:14)//'.NUME'
    call jeveuo(nomres//'.DESM', 'L', vi=desm)
    nddle = desm(4)
    nddli = desm(5)
!
    baseit = nomres//'.PROJM    .PJMBP'
    vsu = nomres//'.PROJM    .PJMIG'
!
! RECUPERATION DES ELEMENTS RELATIFS AU MODES MESURES
    call getvid('DEFINITION', 'MODE_MESURE', iocc=1, scal=modmes, nbret=ibid)
    vref = nomres//'.PROJM    .PJMRF'
    call jeexin(vref, iexist)
    if (iexist .eq. 0) then
        call wkvect(vref, 'G V K16', 5, lref)
        zk16(lref-1+1) = modlms
        zk16(lref-1+2) = nomchp
        zk16(lref-1+3) = basemo
        zk16(lref-1+4) = modmes
        zk16(lref-1+5) = promes
    end if
!
    modid = nomres//'.PROJM    .PJMMM'
    call jeexin(modid, iexist)
    if (iexist .eq. 0) then
        call remome(promes, modmes, nomres)
    end if
    call jeveuo(modid, 'L', lmesu)
    call jelira(modid, 'LONUTI', nbord)
    nbord = nbord/nbmesu
! NBORD : NOMBRE DE NUMERO D'ORDRE
!
    trav = '&TRAV'
    wks = '&WKS'
!
! ===============================
!  TEST : CALCUL INVERSE MATRICE DE PASSAGE DEJA REALISE
! ===============================
    call jeexin(vsu, iexist)
    if (iexist .eq. 0) then
!
! CREATION DE LA BASE RESTREINTE AUX DDL EXTERIEUR
        call getvid('DEFINITION', 'MODELE', iocc=1, scal=modele, nbret=ibid)
        call dismoi('NOM_MAILLA', modele, 'MODELE', repk=ma)
!
        call wkvect(baseit, 'G V R', nddle*nbvecb, lredi)
        call jeveuo(basemo//'           .ORDR', 'L', vi=ordr)
        ch1s = '&BASEIT.CH1S'
!
        do imod = 1, nbvecb
            call rsexch('F', basemo, nomchp, ordr(imod), chamno, &
                        iret)
!
! TRANSFORMATION DE CHAMNO EN CHAM_NO_S : CH1S
            call detrsd('CHAM_NO_S', ch1s)
            call cnocns(chamno, 'V', ch1s)
!
! RECUPERATION DU CHAMP AU NOEUD
            call jeveuo(ch1s//'.CNSK', 'L', vk8=cnsk)
            call jeveuo(ch1s//'.CNSD', 'L', vi=cnsd)
            call jeveuo(ch1s//'.CNSC', 'L', jcnsc)
            call jeveuo(ch1s//'.CNSV', 'L', vr=cnsv)
!
            nbcmpi = cnsd(2)
            nomgd = cnsk(2)
!
            call jenonu(jexnom('&CATA.GD.NOMGD', nomgd), numgd)
            call jeveuo(jexnum('&CATA.GD.NOMCMP', numgd), 'L', llncmp)
!
            if (imod .eq. 1) then
! CREATION DE LA LISTE DES DDL EXTERIEUR
                vnoeud = nomres//'.PROJM    .PJMNO'
                vrange = nomres//'.PROJM    .PJMRG'
                call wkvect(vnoeud, 'G V I', nddle, lnoeud)
                call wkvect(vrange, 'G V K8', nddle, lrange)
            end if
!
            call jeveuo(nu//'.DEEQ', 'L', vi=deeq)
!
            do iddl = 1, nddle
                ino = deeq(nddli*2+(iddl-1)*2+1)
                ico = deeq(nddli*2+(iddl-1)*2+2)
                if (imod .eq. 1) then
                    typ = zk8(llncmp-1+ico)
                    zi(lnoeud-1+iddl) = ino
                    zk8(lrange-1+iddl) = typ
                end if
                ipos = (imod-1)*nddle+iddl
                zr(lredi-1+ipos) = cnsv((ino-1)*nbcmpi+ico)
            end do
!
        end do
!
        call jeecra(vnoeud, 'LONUTI', nddle)
        call jeecra(vrange, 'LONUTI', nddle)
        call jeecra(baseit, 'LONUTI', nddle*nbvecb)
!
! FIN CREATION DE LA BASE RESTREINTE AUX DDL EXTERIEUR (LREDI)
!
! ===============================
! CALCUL MATRICE DE PASSAGE TIT (NOTATION MC)
! ===============================
!
! MESINT : MATRICE DE PASSAGE TIT (CAPTEUR -> INTERFACE)
        mesint = nomres//'.TIT'
        call wkvect(mesint, 'V V R', nddle*nbmesu, lint)
!
! CALCUL DU PSEUDO INVERSE DE LA BASE REDUITE AUX DDL MESURE
!
        if (nbmesu .lt. nbvecb) then
            affici(1) = nbmesu
            affici(2) = nbvecb
            call utmess('F', 'SOUSTRUC_82', ni=2, vali=affici)
        end if
!
        vals = '&VALS'
        u = '&U'
        v = '&V'
!
        call wkvect(vals, 'V V R', nbvecb, lvals)
        call wkvect(u, 'V V R', nbmesu*nbmesu, lu)
        call wkvect(v, 'V V R', nbmesu*nbvecb, lv)
!
! CALCUL PSEUDO INVERSE BASE REDUITE AUX DDL MESURE (LTRAV)
!
        sol = '&SOLUT'
        call wkvect(sol, 'V V R', nbmesu*nbmesu, isol)
!
        call wkvect(wks, 'V V R', nbmesu, lwks)
!
        do iddl = 1, nbmesu
            do jddl = 1, nbmesu
                ipos = (jddl-1)*nbmesu+iddl
                if (iddl .eq. jddl) then
                    zr(isol-1+ipos) = 1.d0
                else
                    zr(isol-1+ipos) = 0.d0
                end if
            end do
        end do
!
        eps = 1.d2*r8prem()
        call wkvect(trav, 'V V R', nbvecb*nbmesu, ltrav)
        do iddl = 1, nbvecb*nbmesu
            zr(ltrav-1+iddl) = zr(lred-1+iddl)
        end do
!
        call rslsvd(nbmesu, nbmesu, nbvecb, zr(ltrav), zr(lvals), &
                    zr(lu), zr(lv), nbmesu, zr(isol), eps, &
                    ier, zr(lwks))
!
        call jedetr(trav)
        if (ier .ne. 0) then
            call utmess('F', 'UTILITAI3_8')
        end if
!
        call wkvect(trav, 'V V R', nbvecb*nbmesu, ltrav)
!
        do imod = 1, nbvecb
            do jddl = 1, nbmesu
                ipos = (jddl-1)*nbmesu+imod
                iposj = (jddl-1)*nbvecb+imod
                zr(ltrav-1+iposj) = zr(isol-1+ipos)
            end do
        end do
!
        call jedetr(vals)
        call jedetr(u)
        call jedetr(v)
!
        call jedetr(wks)
        call jedetr(sol)
!
! FIN CALCUL PSEUDO INVERSE BASE REDUITE AUX DDL MESURE (LTRAV)
!
        do iddl = 1, nddle
            do jddl = 1, nbmesu
                ipos = (jddl-1)*nddle+iddl
                zr(lint-1+ipos) = 0.d0
                do imod = 1, nbvecb
                    iposi = (imod-1)*nddle+iddl
                    iposj = (jddl-1)*nbvecb+imod
                    zr(lint-1+ipos) = zr(lint-1+ipos)+zr(lredi-1+iposi)*zr(ltrav-1+iposj)
                end do
            end do
        end do
!
        call jedetr(trav)
!
! ===============================
! FIN CALCUL MATRICE DE PASSAGE TIT (LINT)
! ===============================
!
!
! ===============================
! CALCUL DU PSEUDO INVERSE DE (TIT*PHI) PAR SVD (NOTATION MC)
! MATRICE DE PASSAGE COORDONNEES GENERALISEES -> DDL_INTERFACE
! ICI PHI : MODES PROPRES IDENTIFIES
! ===============================
!
! CALCUL DU PRODUIT TIT*PHIid : LTRAV
        call wkvect(trav, 'V V R', nddle*nbord, ltrav)
!
        do iddl = 1, nddle
            do jmod = 1, nbord
                ipos = (jmod-1)*nddle+iddl
                zr(ltrav-1+ipos) = 0.d0
                do jddl = 1, nbmesu
                    iposi = (jddl-1)*nddle+iddl
                    iposj = (jmod-1)*nbmesu+jddl
                    zr(ltrav-1+ipos) = zr(ltrav-1+ipos)+zr(lint-1+iposi)*zr(lmesu-1+iposj)
                end do
            end do
        end do
!
        call jedetr(mesint)
!
        if (nddle .lt. nbord) then
            affici(1) = nddle
            affici(2) = nbord
            call utmess('F', 'SOUSTRUC_83', ni=2, vali=affici)
        end if
!
        call wkvect(vals, 'V V R', nbord, lvals)
        call wkvect(u, 'V V R', nddle*nddle, lu)
        call wkvect(v, 'V V R', nddle*nbord, lv)
!
        call wkvect(sol, 'V V R', nddle*nddle, isol)
!
        call wkvect(wks, 'V V R', nddle, lwks)
!
        do iddl = 1, nddle
            do jddl = 1, nddle
                ipos = (jddl-1)*nddle+iddl
                if (iddl .eq. jddl) then
                    zr(isol-1+ipos) = 1.d0
                else
                    zr(isol-1+ipos) = 0.d0
                end if
            end do
        end do
!
! CALCUL DE L'INVERSE DU PRODUIT TIT*PHIid : LVSU
!
        call wkvect(vsu, 'G V R', nbord*nddle, lvsu)
!
        call rslsvd(nddle, nddle, nbord, zr(ltrav), zr(lvals), &
                    zr(lu), zr(lv), nddle, zr(isol), eps, &
                    ier, zr(lwks))
!
        if (ier .ne. 0) then
            call utmess('F', 'UTILITAI3_8')
        end if
!
        call jedetr(wks)
!
        do imod = 1, nbord
            do jddl = 1, nddle
                ipos = (jddl-1)*nbord+imod
                iposj = (jddl-1)*nddle+imod
                zr(lvsu-1+ipos) = zr(isol-1+iposj)
            end do
        end do
!
        call jeecra(vsu, 'LONUTI', nbord*nddle)
!
        call jedetr(sol)
!
        call jedetr(trav)
!
        call jedetr(vals)
        call jedetr(u)
        call jedetr(v)
!
    else
        call jeveuo(vsu, 'L', lvsu)
    end if
!
! ===============================
! FIN TEST : CALCUL INVERSE MATRICE DE PASSAGE DEJA REALISE
! ===============================
!
!  CALCUL DE VSUT : LTRAV
    call wkvect(trav, 'V V R', nddle*nbord, ltrav)
    call wkvect(wks, 'V V R', nbord*nddle, lwks)
!
    do imod = 1, nbord
        call rsadpa(modmes, 'L', 1, 'MASS_GENE', imod, &
                    0, sjv=ipuls, styp=k8bid)
        masg = zr(ipuls)
        do iddl = 1, nddle
            ipos = (iddl-1)*nbord+imod
            iposi = (imod-1)*nddle+iddl
            zr(ltrav-1+iposi) = zr(lvsu-1+ipos)*sqrt(masg)
            zr(lwks-1+ipos) = zr(lvsu-1+ipos)*sqrt(masg)
        end do
    end do
!
! ===============================
! CALCUL DES TERMES DE LA <MATRICE ELEMENTAIRE>
! ===============================
!
    maelm = nomres//'.MAEL_M'
    call wkvect(maelm, 'V V R', nddle*nddle, lmaelm)
!
    do iddl = 1, nddle
        do jddl = 1, nddle
            ipos = (jddl-1)*nddle+iddl
            zr(lmaelm-1+ipos) = 0.d0
            do imod = 1, nbord
                iposi = (imod-1)*nddle+iddl
                iposj = (jddl-1)*nbord+imod
                zr(lmaelm-1+ipos) = zr(lmaelm-1+ipos)+zr(ltrav-1+iposi)*zr(lwks-1+iposj)
            end do
        end do
    end do
!
! RANGEMENT DES RESULTATS DANS IAMPEE
!
    do iddl = 1, nddle
        do jddl = 1, iddl
            ipos = (iddl-1)*iddl/2+jddl
            iposi = (jddl-1)*nddle+iddl
            zr(iampee-1+ipos) = zr(lmaelm-1+iposi)
        end do
    end do
!
    call jedetr(trav)
    call jedetr(wks)
    call jedetr(maelm)
!
! ===============================
! FIN CALCUL DES TERMES DES <MATRICE ELEMENTAIRE>
! ===============================
!
    call jedema()
!
end subroutine
