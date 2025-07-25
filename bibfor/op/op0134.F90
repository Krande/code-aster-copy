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

subroutine op0134()
    implicit none
!     CALCUL D'UNE FONCTION INTERPRETEE
!     ------------------------------------------------------------------
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterfort/gettco.h"
#include "asterfort/assert.h"
#include "asterfort/calcfo.h"
#include "asterfort/calcna.h"
#include "asterfort/foattr.h"
#include "asterfort/foimpr.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/infmaj.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/ordonn.h"
#include "asterfort/titre.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
    integer(kind=8) :: ifm, niv, n1, nbvalp, nbvalf, lvalp, lvalf, nbnova, lprol
    aster_logical :: compl
    character(len=8) :: nopn, nopf
    character(len=16) :: nomcmd, typres
    character(len=19) :: nomfon, nomfin, listp, listf, typco
    character(len=24) :: noparp, noparf, valk(3)
    character(len=24), pointer :: nova(:) => null()
!     ------------------------------------------------------------------
!
    call jemarq()
!
    call infmaj()
    call infniv(ifm, niv)
!
    call getres(nomfon, typres, nomcmd)
!
    call getvid(' ', 'FONCTION', scal=nomfin, nbret=n1)
    call gettco(nomfin, typco)
!
! --- LISTE DES VALEURS DU PARAMETRE
!
    call getvr8(' ', 'VALE_PARA', nbval=0, nbret=n1)
    if (n1 .ne. 0) then
        nbvalp = -n1
        call wkvect('&&OP0134.VALP', 'V V R', nbvalp, lvalp)
        call getvr8(' ', 'VALE_PARA', nbval=nbvalp, vect=zr(lvalp), nbret=n1)
    else
        call getvid(' ', 'LIST_PARA', scal=listp, nbret=n1)
        call jeveuo(listp//'.VALE', 'L', lvalp)
        call jelira(listp//'.VALE', 'LONUTI', nbvalp)
    end if
!
! --- NAPPE OU FONCTION
!
    compl = .false.
    if (typco(1:7) .eq. 'FORMULE') then
        if (typco(1:9) .eq. 'FORMULE_C') compl = .true.
        call jelira(nomfin//'.NOVA', 'LONUTI', nbnova)
        call jeveuo(nomfin//'.NOVA', 'L', vk24=nova)
        if (nbnova .eq. 1) then
            noparp = nova(1)
        else if (nbnova .eq. 2) then
            noparp = nova(1)
            noparf = nova(2)
        end if
!
    else if (typco(1:8) .eq. 'FONCTION') then
        if (typco(1:10) .eq. 'FONCTION_C') compl = .true.
        nbnova = 1
        call jeveuo(nomfin//'.PROL', 'L', lprol)
        noparp = zk24(lprol+2)
!
    else if (typco(1:5) .eq. 'NAPPE') then
        nbnova = 2
        call jeveuo(nomfin//'.PROL', 'L', lprol)
        noparp = zk24(lprol+2)
        noparf = zk24(lprol+6)
!
    else
        ASSERT(.false.)
    end if
!
!
    if (nbnova .eq. 1) then
! ------------------------------------------------------------------
!                 FONCTION
! ------------------------------------------------------------------
        call calcfo(compl, nomfin, nomfon, nbvalp, zr(lvalp), &
                    noparp)
!
    else if (nbnova .eq. 2) then
! ------------------------------------------------------------------
!                 NAPPE
! ------------------------------------------------------------------
        call getvr8(' ', 'VALE_PARA_FONC', nbval=0, nbret=n1)
        if (n1 .ne. 0) then
            nbvalf = -n1
            call wkvect('&&OP0134.VALF', 'V V R', nbvalf, lvalf)
            call getvr8(' ', 'VALE_PARA_FONC', nbval=nbvalf, vect=zr(lvalf), nbret=n1)
        else
            call getvid(' ', 'LIST_PARA_FONC', scal=listf, nbret=n1)
            if (n1 .ne. 0) then
                call jeveuo(listf//'.VALE', 'L', lvalf)
                call jelira(listf//'.VALE', 'LONUTI', nbvalf)
            else
                call utmess('F', 'FONCT0_49')
            end if
        end if
!
!        VERIFIER LA COHERENCE DES NOMS DES PARAMETRES
        call getvtx(' ', 'NOM_PARA', scal=nopn, nbret=n1)
!        FACULTATIF
        if (n1 .ne. 0 .and. nopn .ne. noparp) then
            valk(1) = nomfin
            valk(2) = noparp
            valk(3) = nopn
            if (typco(1:7) .eq. 'FORMULE') then
                call utmess('F', 'FONCT0_58', nk=3, valk=valk)
            else
                call utmess('F', 'FONCT0_59', nk=3, valk=valk)
            end if
        end if
!
        call getvtx(' ', 'NOM_PARA_FONC', scal=nopf, nbret=n1)
!        OBLIGATOIRE
        ASSERT(n1 .eq. 1)
        if (nopf .ne. noparf) then
            valk(1) = nomfin
            valk(2) = noparf
            valk(3) = nopf
            if (typco(1:7) .eq. 'FORMULE') then
                call utmess('F', 'FONCT0_60', nk=3, valk=valk)
            else
                call utmess('F', 'FONCT0_61', nk=3, valk=valk)
            end if
        end if
!
        call calcna(nomfin, nomfon, nbvalp, zr(lvalp), noparp, &
                    nbvalf, zr(lvalf), noparf)
!
    else
!
        call utmess('F', 'FONCT0_48')
!
    end if
!
! --- SURCHARGE EVENTUELLE DU .PROL
!
    call foattr(' ', 1, nomfon)
!
! --- VERIFICATION QU'ON A BIEN CREER UNE FONCTION ---
!     ET REMISE DES ABSCISSES EN ORDRE CROISSANT
!
    call ordonn(nomfon, 0)
!
    call titre()
    if (niv .gt. 1) call foimpr(nomfon, niv, ifm, 0, listp)
!
    call jedema()
end subroutine
