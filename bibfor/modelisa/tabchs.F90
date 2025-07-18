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

subroutine tabchs(tabin, typchs, base, nomgd, ma, &
                  chs)
!
!
!     TRAITEMENT DE COMMANDE:   CREA_CHAMP / OPTION: 'EXTR' / TABLE
!
!     CREATION D'UN CHAMP SIMPLE A PARTIR D'UNE TABLE
!
!
!     IN : TABIN  : NOM DE LA TABLE
!     IN : TYPCHS : TYPE DU CHAMP SIMPLE (NOEU/ELEM/ELNO/ELGA)
!     IN : BASE   : BASE DE CREATION (G/V)
!     IN : NOMGD  : NOM DE LA GRANDEUR
!     IN : MA     : NOM DU MAILLAGE
!     IN/JXOUT : CHS: NOM DU CHAMP SIMPLE A CREER
!
    implicit none
!
!     ------------------------------------------------------------------
! 0.1. ==> ARGUMENT
!
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/cescre.h"
#include "asterfort/cesexi.h"
#include "asterfort/cnscre.h"
#include "asterfort/dismoi.h"
#include "asterfort/indiis.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/tbexip.h"
#include "asterfort/tbexve.h"
#include "asterfort/utmess.h"
#include "asterfort/verigd.h"
#include "asterfort/verima.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "asterfort/char8_to_int.h"
!
    character(len=1) :: base
    character(len=8) :: nomgd, ma
    character(len=16) :: typchs
    character(len=19) :: chs, tabin
!
!
! 0.2. ==> COMMUNS
!
!
!      ==> VARIABLES LOCALES
    integer(kind=8) :: ncmp, jcnsl, i
    integer(kind=8) :: vali(2), nblig, isp
    integer(kind=8) :: nbval, ier
    integer(kind=8) :: nuno, numa, nbma, jcesd, nbssp
    integer(kind=8) :: jcesl, jcesc, iad
    integer(kind=8) :: nbcol, nbno, ksp, kpt, jcon1, jcon2
    integer(kind=8) :: jcolma, jcolno, jcolpt, jcolsp, ipt
    integer(kind=8) :: icmp, ili, iret, jobj2, jobj3
    character(len=8) :: nono, tsca, noma
    character(len=24) :: objlg, objr, objtmp
    character(len=24) :: valk(3)
    aster_logical :: lmail, lnoeu, lpoin, lspoin
    character(len=24), pointer :: ncmp1(:) => null()
    integer(kind=8), pointer :: ncmp2(:) => null()
    integer(kind=8), pointer :: pg_tot(:) => null()
    integer(kind=8), pointer :: sp_tot(:) => null()
    character(len=24), pointer :: tblp(:) => null()
    real(kind=8), pointer :: cnsv(:) => null()
    integer(kind=8), pointer :: tbnp(:) => null()
    real(kind=8), pointer :: cesv(:) => null()
    aster_logical :: lcolle, lcolle2
! ---------------------------------------------------------------------
!
!
    call jemarq()
!
    call jeveuo(tabin//'.TBNP', 'L', vi=tbnp)
    nbcol = tbnp(1)
    nblig = tbnp(2)
    call jeveuo(tabin//'.TBLP', 'L', vk24=tblp)
!
    call dismoi('NB_MA_MAILLA', ma, 'MAILLAGE', repi=nbma)
    call dismoi('TYPE_SCA', nomgd, 'GRANDEUR', repk=tsca)
    ASSERT(tsca .eq. 'R')
    call jeveuo(ma//'.CONNEX', 'L', jcon1)
    call jeveuo(jexatr(ma//'.CONNEX', 'LONCUM'), 'L', jcon2)
!
    lcolle = .false.
    call jeexin(ma//'.NOMNOE', ier)
    if (ier .ne. 0) then
        lcolle = .true.
    end if
    lcolle2 = .false.
    call jeexin(ma//'.NOMMAI', ier)
    if (ier .ne. 0) then
        lcolle2 = .true.
    end if
!
!     1. VERIFICATION DES PARAMETRES DE LA TABLE :
!     --------------------------------------------
    call tbexip(tabin, 'MAILLE', lmail, tsca)
    call tbexip(tabin, 'NOEUD', lnoeu, tsca)
    call tbexip(tabin, 'POINT', lpoin, tsca)
    call tbexip(tabin, 'SOUS_POINT', lspoin, tsca)
!
    valk(1) = tabin(1:8)
    valk(2) = typchs
!
!     PRESENCE DU PARAMETRE MAILLE
    if (typchs .eq. 'EL' .and. .not. lmail) then
        call utmess('F', 'MODELISA9_1', nk=2, valk=valk)
    end if
!
!     PRESENCE/ABSENCE POUR CHAMPS NOEU :
    if (typchs .eq. 'NOEU') then
        if (.not. lnoeu) then
            call utmess('F', 'MODELISA9_1', nk=2, valk=valk)
        end if
        if (lmail .or. lpoin .or. lspoin) then
            call utmess('F', 'MODELISA9_1', nk=2, valk=valk)
        end if
    end if
!
!     PRESENCE/ABSENCE POUR CHAMPS ELGA :
    if (typchs .eq. 'ELGA') then
        if (.not. lpoin .or. lnoeu) then
            call utmess('F', 'MODELISA9_1', nk=2, valk=valk)
        end if
    end if
!
!     PRESENCE/ABSENCE POUR CHAMPS ELNO :
    if (typchs .eq. 'ELNO') then
        if (.not. lpoin .and. .not. lnoeu) then
            call utmess('F', 'MODELISA9_1', nk=2, valk=valk)
        end if
        if (lpoin .and. lnoeu) then
            call utmess('F', 'MODELISA9_1', nk=2, valk=valk)
        end if
    end if
!
!     PRESENCE/ABSENCE POUR CHAMPS ELEM :
    if (typchs .eq. 'ELEM') then
        if (lpoin .or. lnoeu) then
            call utmess('F', 'MODELISA9_1', nk=2, valk=valk)
        end if
    end if
!
!
!     2. RECUPERATION DES COLONNES DE LA TABLE :
!     -------------------------------------------
!     -- 2.1 COLONNE 'NOEUD' :
    if (lnoeu) then
!        ON VERIFIE QUE LES NOEUDS FOURNIS DANS LA TABLE
!        APPARTIENNENT AU MAILLAGE
        objtmp = '&&TABCHS.NOEUD'
        call tbexve(tabin, 'NOEUD', objtmp, 'V', nbval, tsca)
        if (tsca .eq. 'K8') then
            call jeveuo(objtmp, 'L', jcolno)
            call verima(ma, zk8(jcolno), nbval, 'NOEUD')
        else
            valk(2) = tsca
            call utmess('F', 'MODELISA9_7', nk=2, valk=valk)
        end if
    end if
!
!     -- 2.2 COLONNE 'MAILLE' :
    if (lmail) then
!        ON VERIFIE QUE LES MAILLES FOURNIES DANS LA TABLE
!        APPARTIENNENT AU MAILLAGE
        objtmp = '&&TABCHS.MAILLE'
        call tbexve(tabin, 'MAILLE', objtmp, 'V', nbval, tsca)
        if (tsca .eq. 'K8') then
            call jeveuo(objtmp, 'L', jcolma)
            call verima(ma, zk8(jcolma), nbval, 'MAILLE')
        else
            valk(2) = tsca
            call utmess('F', 'MODELISA9_7', nk=2, valk=valk)
        end if
    end if
!
!     -- 2.3 COLONNE 'POINT' :
    if (lpoin) then
        objtmp = '&&TABCHS.POINT'
        call tbexve(tabin, 'POINT', objtmp, 'V', nbval, tsca)
        call jeveuo(objtmp, 'L', jcolpt)
    end if
!
!     -- 2.4 COLONNE 'SOUS_POINT' :
    if (lspoin) then
        objtmp = '&&TABCHS.SPOINT'
        call tbexve(tabin, 'SOUS_POINT', objtmp, 'V', nbval, tsca)
        call jeveuo(objtmp, 'L', jcolsp)
    end if
!
!
!     -- ON REPERE LES COLONNES QUI CORRESPONDENT AUX CMPS.
!        CE SONT CELLES QUI NE SONT PAS : MAILLE, NOEUD, ...
    AS_ALLOCATE(vk24=ncmp1, size=nbcol)
    AS_ALLOCATE(vi=ncmp2, size=nbcol)
    ncmp = 0
    do i = 1, nbcol
        if (tblp(1+4*(i-1)) .ne. 'MAILLE' .and. tblp(1+4*(i-1)) .ne. 'NOEUD' .and. &
            tblp(1+4*(i-1)) .ne. 'POINT' .and. tblp(1+4*(i-1)) .ne. 'SOUS_POINT') then
            ncmp = ncmp+1
            ncmp1(ncmp) = tblp(1+4*(i-1))
            ncmp2(ncmp) = i
        end if
    end do
!
!     ON VERIFIE QUE LE NOM ET LE TYPE DES COMPOSANTES
!        DE LA TABLE CORRESPONDENT A LA GRANDEUR LUE
    call verigd(nomgd, ncmp1, ncmp, iret)
    if (iret .ne. 0) then
        call utmess('F', 'MODELISA9_2', nk=ncmp, valk=ncmp1)
    end if
!
!
!
    if (typchs .eq. 'NOEU') then
!     ------------------------------------
!
! ---    CREATION DU CHAM_NO_S
        call cnscre(ma, nomgd, ncmp, ncmp1, base, &
                    chs)
!
! ---    REMPLISSAGE DU CHAM_S
        call jeveuo(chs//'.CNSV', 'E', vr=cnsv)
        call jeveuo(chs//'.CNSL', 'E', jcnsl)
!
        do icmp = 1, ncmp
            objlg = tblp(1+4*(ncmp2(icmp)-1)+3)
            call jeveuo(objlg, 'L', jobj2)
            objr = tblp(1+4*(ncmp2(icmp)-1)+2)
            call jeveuo(objr, 'L', jobj3)
            do ili = 1, nblig
                if (zi(jobj2+ili-1) .eq. 1) then
                    nono = zk8(jcolno+ili-1)
                    nuno = char8_to_int(nono, lcolle, ma, "NOEUD")
                    ASSERT(nuno .gt. 0)
                    cnsv(1+(nuno-1)*ncmp+icmp-1) = zr(jobj3+ili-1)
                    zl(jcnsl+(nuno-1)*ncmp+icmp-1) = .true.
                end if
            end do
        end do
!
!
!
    else if (typchs(1:2) .eq. 'EL') then
!     ------------------------------------
!
        if (typchs .eq. 'ELNO') then
!          POUR LES CHAMPS ELNO :
!           - SI LNOEU, ON CALCULE '&&TABCHS.POINT'
!           - ON VERIFIE QUE LE NUMERO DE POINT EST POSSIBLE
            if (lnoeu) then
                ASSERT(.not. lpoin)
                call wkvect('&&TABCHS.POINT', 'V V I', nblig, jcolpt)
            end if
            do ili = 1, nblig
                noma = zk8(jcolma+ili-1)
                numa = char8_to_int(noma, lcolle2, ma, "MAILLE")
                nbno = zi(jcon2-1+numa+1)-zi(jcon2-1+numa)
                if (lpoin) then
                    ipt = zi(jcolpt-1+ili)
                else
                    ASSERT(lnoeu)
                    nono = zk8(jcolno-1+ili)
                    nuno = char8_to_int(nono, lcolle, ma, "NOEUD")
                    ipt = indiis(zi(jcon1-1+zi(jcon2-1+numa)), nuno, 1, &
                                 nbno)
                    zi(jcolpt-1+ili) = ipt
                end if
                if (ipt .eq. 0 .or. ipt .gt. nbno) then
                    valk(1) = tabin
                    valk(2) = noma
                    valk(3) = nono
                    vali(1) = ipt
                    vali(2) = nbno
                    call utmess('F', 'MODELISA9_5', nk=3, valk=valk, ni=2, &
                                vali=vali)
                end if
            end do
        end if
!
!
!       CALCUL DU NOMBRE DE SOUS_POINT PAR ELEMENT (&&TABCHS.SP_TOT):
!       CALCUL DE NBSSP : MAX DU NOMBRE DE SOUS_POINT
        if (lspoin) then
            AS_ALLOCATE(vi=sp_tot, size=nbma)
            nbssp = 1
            do ili = 1, nblig
                ksp = zi(jcolsp+ili-1)
                ASSERT(ksp .gt. 0)
                nbssp = max(nbssp, ksp)
                noma = zk8(jcolma+ili-1)
                numa = char8_to_int(noma, lcolle2, ma, "MAILLE")
                ASSERT(numa .gt. 0)
                sp_tot(numa) = max(sp_tot(numa), ksp)
            end do
        else
            nbssp = 1
        end if
!
!
!       CALCUL DU NOMBRE DE POINTS DE GAUSS PAR ELEMENT
!       (&&TABCHS.PG_TOT):
        if (typchs .eq. 'ELGA') then
            AS_ALLOCATE(vi=pg_tot, size=nbma)
            do ili = 1, nblig
                kpt = zi(jcolpt+ili-1)
                ASSERT(kpt .gt. 0)
                noma = zk8(jcolma+ili-1)
                numa = char8_to_int(noma, lcolle2, ma, "MAILLE")
                pg_tot(numa) = max(pg_tot(numa), kpt)
            end do
        end if
!
!
!       CREATION DU CHAM_ELEM_S VIERGE :
        if (nbssp .eq. 1) then
            if (typchs .eq. 'ELNO' .or. typchs .eq. 'ELEM') then
                call cescre(base, chs, typchs, ma, nomgd, &
                            ncmp, ncmp1, [-1], [-1], [-ncmp])
            else if (typchs .eq. 'ELGA') then
                call cescre(base, chs, typchs, ma, nomgd, &
                            ncmp, ncmp1, pg_tot, [-1], [-ncmp])
            end if
        else
            if (typchs .eq. 'ELNO' .or. typchs .eq. 'ELEM') then
                call cescre(base, chs, typchs, ma, nomgd, &
                            ncmp, ncmp1, [-1], sp_tot, [-ncmp])
            else if (typchs .eq. 'ELGA') then
                call cescre(base, chs, typchs, ma, nomgd, &
                            ncmp, ncmp1, pg_tot, sp_tot, [-ncmp])
            end if
        end if
!
!
! ---   REMPLISSAGE DU CHAM_S
        call jeveuo(chs//'.CESD', 'L', jcesd)
        call jeveuo(chs//'.CESV', 'E', vr=cesv)
        call jeveuo(chs//'.CESL', 'E', jcesl)
        call jeveuo(chs//'.CESC', 'L', jcesc)
!
!
        do icmp = 1, ncmp
            objlg = tblp(1+4*(ncmp2(icmp)-1)+3)
            call jeveuo(objlg, 'L', jobj2)
            objr = tblp(1+4*(ncmp2(icmp)-1)+2)
            call jeveuo(objr, 'L', jobj3)
!
            do ili = 1, nblig
                if (zi(jobj2+ili-1) .eq. 0) goto 70
!
                noma = zk8(jcolma+ili-1)
                numa = char8_to_int(noma, lcolle2, ma, "MAILLE")
!
                ipt = 1
                if (lpoin) ipt = zi(jcolpt+ili-1)
!
                isp = 1
                if (lspoin) isp = zi(jcolsp+ili-1)
!
                nono = ' '
                if (lnoeu) then
                    nono = zk8(jcolno+ili-1)
                    ipt = zi(jcolpt+ili-1)
                end if
!
                call cesexi('S', jcesd, jcesl, numa, ipt, &
                            isp, icmp, iad)
                ASSERT(iad .ne. 0)
                if (iad .lt. 0) then
                    iad = -iad
                else
                    valk(1) = tabin(1:8)
                    valk(2) = noma
                    valk(3) = nono
                    vali(1) = ipt
                    vali(2) = isp
                    call utmess('F', 'MODELISA9_6', nk=3, valk=valk, ni=2, &
                                vali=vali)
                end if
                cesv(iad) = zr(jobj3+ili-1)
                zl(jcesl+iad-1) = .true.
70              continue
            end do
        end do
    end if
!
    call jedetr('&&TABCHS.MAILLE')
    call jedetr('&&TABCHS.NOEUD')
    call jedetr('&&TABCHS.POINT')
    call jedetr('&&TABCHS.SPOINT')
    AS_DEALLOCATE(vk24=ncmp1)
    AS_DEALLOCATE(vi=ncmp2)
    AS_DEALLOCATE(vi=pg_tot)
    AS_DEALLOCATE(vi=sp_tot)
!
    call jedema()
end subroutine
