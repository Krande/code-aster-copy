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

subroutine focrr0(nomfon, interp, base, resu, nomcha, &
                  maille, noeud, cmp, npoint, nusp, &
                  ivari, nomvari, nbordr, lordr)
    implicit none
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelibe.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/lxlgut.h"
#include "asterfort/posddl.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsnopa.h"
#include "asterfort/utch19.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
    integer(kind=8), intent(in) :: nbordr, lordr(*), npoint, nusp, ivari
    character(len=1), intent(in) :: base
    character(len=8), intent(in) :: interp, maille, noeud, cmp
    character(len=16), intent(in) :: nomcha
    character(len=16), intent(in) :: nomvari
    character(len=19), intent(in) :: nomfon, resu
!     RECUPERATION D'UNE FONCTION DANS UNE STRUCTURE "RESULTAT"
!                  POUR TOUS LES NUMEROS D'ORDRE
!     ------------------------------------------------------------------
! IN : NOMFON : NOM DE LA FONCTION
! IN  : INTERP : TYPE D'INTERPOLATION DE LA FONCTION
! IN  : BASE   : BASE OU L'ON CREE LA FONCTION
! IN  : RESU   : NOM DE LA STRUCTURE RESULTAT
! IN  : NOMCHA : NOM DU CHAMP
! IN  : NOEUD  : NOEUD
! IN  : MAILLE : MAILE
! IN  : CMP    : COMPOSANTE
! IN  : NPOINT : NUMERO DU POINT ( CAS DES CHAM_ELEMS )
! IN  : NUSP   : NUMERO DU SOUS-POINT ( CAS DES CHAM_ELEMS )
! IN  : IVARI  : NUMERO DE LA CMP (POUR VARI_R)
! IN  : NOMVARI: NOM DE LA CMP (POUR VARI_R)
!     ------------------------------------------------------------------
    character(len=1) :: type
    character(len=24) :: valk(2)
    character(len=4) :: typch2
    character(len=8) :: k8b, noma, nogd, noeuz
    character(len=16) :: nomcmd, typcon, nomacc, typcha, typres
    character(len=19) :: profch, profc2, cham19
    complex(kind=8) :: valc
    integer(kind=8) :: npoinz, nuspz
    integer(kind=8) :: ibid, iddl, ie, ii, inoeud, iordr
    integer(kind=8) :: iret, lacce, lfon, lg1, lg2, lpro
    integer(kind=8) :: lvale, lvar, nbacc, vali
    real(kind=8) :: valr
    character(len=16), pointer :: acces(:) => null()
!-----------------------------------------------------------------------
    call jemarq()
    call getres(k8b, typcon, nomcmd)
! GETTCO ne fonctionne pas avec les noms compose issus de sensibilite
!      CALL GETTCO(RESU,TYPRES)
    call dismoi('TYPE_RESU', resu(1:8), 'RESULTAT', repk=typres)
!
    call rsnopa(resu, 0, '&&FOCRR0.VAR.ACCES', nbacc, ibid)
    call jeexin('&&FOCRR0.VAR.ACCES', iret)
    if (iret .gt. 0) then
        call jeveuo('&&FOCRR0.VAR.ACCES', 'E', vk16=acces)
        nomacc = acces(1)
    else
        nomacc = ' '
    end if
!
!     --- REMPLISSAGE DU .PROL ---
    ASSERT(lxlgut(nomfon) .le. 24)
    call wkvect(nomfon//'.PROL', base//' V K24', 6, lpro)
    if (typres(1:10) .eq. 'DYNA_HARMO') then
        zk24(lpro) = 'FONCT_C'
    else
        zk24(lpro) = 'FONCTION'
    end if
    zk24(lpro+1) = interp
    zk24(lpro+2) = nomacc
    if (nomvari .ne. ' ') then
        ASSERT(ivari .gt. 0)
        zk24(lpro+3) = nomvari
    else
        zk24(lpro+3) = cmp
    end if
    zk24(lpro+4) = 'EE      '
    zk24(lpro+5) = nomfon
!
    if (typres(1:10) .eq. 'DYNA_HARMO') then
        call wkvect(nomfon//'.VALE', base//' V R', 3*nbordr, lvar)
    else
        call wkvect(nomfon//'.VALE', base//' V R', 2*nbordr, lvar)
    end if
    lfon = lvar+nbordr
    call rsexch('F', resu, nomcha, lordr(1), cham19, &
                ie)
    call dismoi('TYPE_SUPERVIS', cham19, 'CHAMP', repk=typcha)
!
    if (typcha(1:7) .eq. 'CHAM_NO') then
        call dismoi('NUME_EQUA', cham19, 'CHAM_NO', repk=profch)
        call dismoi('NOM_MAILLA', cham19, 'CHAM_NO', repk=noma)
        call posddl('CHAM_NO', cham19, noeud, cmp, inoeud, &
                    iddl)
        if (inoeud .eq. 0) then
            lg1 = lxlgut(noeud)
            call utmess('F', 'UTILITAI_92', sk=noeud(1:lg1))
        else if (iddl .eq. 0) then
            lg1 = lxlgut(noeud)
            lg2 = lxlgut(cmp)
            valk(1) = cmp(1:lg2)
            valk(2) = noeud(1:lg1)
            call utmess('F', 'UTILITAI_93', nk=2, valk=valk)
        end if
        ii = 0
        do iordr = 1, nbordr
            call jemarq()
!
!           --- EXTRACTION DU CHAMP ET DE LA VALEUR DE L'ACCES ----
            call rsexch(' ', resu, nomcha, lordr(iordr), cham19, &
                        ie)
            if (ie .eq. 0) then
                call dismoi('NUME_EQUA', cham19, 'CHAM_NO', repk=profc2)
                if (profc2 .ne. profch) then
                    profch = profc2
                    call posddl('CHAM_NO', cham19, noeud, cmp, inoeud, &
                                iddl)
                    if (inoeud .eq. 0) then
                        lg1 = lxlgut(noeud)
                        call utmess('F', 'UTILITAI_92', sk=noeud(1:lg1))
                    else if (iddl .eq. 0) then
                        lg1 = lxlgut(noeud)
                        lg2 = lxlgut(cmp)
                        valk(1) = cmp(1:lg2)
                        valk(2) = noeud(1:lg1)
                        call utmess('F', 'UTILITAI_93', nk=2, valk=valk)
                    end if
                end if
                call rsadpa(resu, 'L', 1, nomacc, lordr(iordr), &
                            0, sjv=lacce, styp=k8b)
                call jeveuo(cham19//'.VALE', 'L', lvale)
                if (typres(1:10) .eq. 'DYNA_HARMO') then
                    zr(lvar+iordr-1) = zr(lacce)
                    zr(lfon+ii) = dble(zc(lvale+iddl-1))
                    ii = ii+1
                    zr(lfon+ii) = dimag(zc(lvale+iddl-1))
                    ii = ii+1
                else
                    zr(lvar+iordr-1) = zr(lacce)
                    zr(lfon+iordr-1) = zr(lvale+iddl-1)
                end if
                call jelibe(cham19//'.VALE')
            end if
            call jedema()
        end do
!
    else if (typcha(1:9) .eq. 'CHAM_ELEM') then
        noeuz = noeud
        npoinz = npoint
        nuspz = nusp
! ---    VERIFICATION DE LA PRESENCE DES MOTS CLE GROUP_MA (OU MAILLE)
! ---    ET GROUP_NO (OU NOEUD OU POINT) DANS LE CAS D'UN CHAM_ELEM
!        -------------------------------------------------------------
        call dismoi('TYPE_CHAMP', cham19, 'CHAMP', repk=typch2)
        if (typch2 .eq. 'ELEM') then
            npoinz = 1
            nuspz = 1
            noeuz = ' '
            if (maille .eq. ' ') then
                call utmess('F', 'CHAMPS_11')
            end if
        else if (typch2 .eq. 'ELNO') then
            nuspz = 1
            if (maille .eq. ' ' .or. (noeud .eq. ' ' .and. npoint .eq. 0)) then
                call utmess('F', 'CHAMPS_12')
            end if
        else
            if (maille .eq. ' ' .or. npoint .eq. 0) then
                call utmess('F', 'CHAMPS_13')
            end if
        end if
        call dismoi('NOM_MAILLA', cham19, 'CHAM_ELEM', repk=noma)
        call dismoi('NOM_GD', cham19, 'CHAM_ELEM', repk=nogd)
        call dismoi('TYPE_SCA', nogd, 'GRANDEUR', repk=type)
        ii = 0
        do iordr = 1, nbordr
            call jemarq()
!
!           --- EXTRACTION DU CHAMP ET DE LA VALEUR DE L'ACCES ----
            call rsexch(' ', resu, nomcha, lordr(iordr), cham19, &
                        ie)
            if (ie .eq. 0) then
                call rsadpa(resu, 'L', 1, nomacc, lordr(iordr), &
                            0, sjv=lacce, styp=k8b)
                call utch19(cham19, noma, maille, noeuz, npoinz, &
                            nuspz, ivari, cmp, type, valr, &
                            valc, vali, iret)
                if (iret .eq. 0) then
                    zr(lvar+iordr-1) = zr(lacce)
                    if (type .eq. 'R') then
                        zr(lfon+iordr-1) = valr
                    else
                        zr(lfon+ii) = dble(valc)
                        ii = ii+1
                        zr(lfon+ii) = dimag(valc)
                        ii = ii+1
                    end if
                end if
            end if
            call jedema()
        end do
    else
        call utmess('F', 'UTILITAI_94', sk=typcha)
    end if
    call jedetr('&&FOCRR0.VAR.ACCES')
!
    call jedema()
end subroutine
