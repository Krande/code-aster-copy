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

subroutine rfrcha()
    implicit none
!     OPERATEUR "RECU_FONCTION"
!     ------------------------------------------------------------------
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterfort/dismoi.h"
#include "asterfort/foattr.h"
#include "asterfort/focste.h"
#include "asterfort/foimpr.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/infmaj.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/lxlgut.h"
#include "asterfort/ordonn.h"
#include "asterfort/posddl.h"
#include "asterfort/titre.h"
#include "asterfort/utch19.h"
#include "asterfort/utcmp1.h"
#include "asterfort/utmess.h"
#include "asterfort/utnono.h"
    integer(kind=8) ::  lg1, lg2, iddl, inoeud, nch
    integer(kind=8) :: n1, iret, ivari
    integer(kind=8) :: nm, ngm, npoint, np, nn
    integer(kind=8) :: ngn, nc, ifm, niv, nusp
    real(kind=8) :: epsi, valr
    complex(kind=8) :: valc
    character(len=1) :: type
    character(len=24) :: valk(2)
    character(len=4) :: typch2
    character(len=8) :: k8b, crit, maille, noma, intres
    character(len=8) :: noeud, cmp, nomgd
    character(len=16) :: nomcmd, typcon, typcha, nom_vari
    character(len=19) :: nomfon, cham19
    character(len=24) :: nogno, nogma
    integer(kind=8) :: vali
    real(kind=8), pointer :: vale(:) => null()
!     ------------------------------------------------------------------
    call jemarq()
! --- RECUPERATION DU NIVEAU D'IMPRESSION
    call infmaj()
    call infniv(ifm, niv)
!
    call getres(nomfon, typcon, nomcmd)
!
    call getvtx(' ', 'CRITERE', scal=crit, nbret=n1)
    call getvr8(' ', 'PRECISION', scal=epsi, nbret=n1)
    intres = 'NON     '
    call getvtx(' ', 'INTERP_NUME', scal=intres, nbret=n1)
!
    npoint = 0
    cmp = ' '
    noeud = ' '
    maille = ' '
    nogma = ' '
    nogno = ' '
    call getvtx(' ', 'MAILLE', scal=maille, nbret=nm)
    call getvtx(' ', 'GROUP_MA', scal=nogma, nbret=ngm)
    call getvis(' ', 'SOUS_POINT', scal=nusp, nbret=np)
    if (np .eq. 0) nusp = 0
    call getvis(' ', 'POINT', scal=npoint, nbret=np)
    call getvtx(' ', 'NOEUD', scal=noeud, nbret=nn)
    call getvtx(' ', 'GROUP_NO', scal=nogno, nbret=ngn)
!
!     -----------------------------------------------------------------
!                      --- CAS D'UN CHAM_GD ---
!     -----------------------------------------------------------------
    call getvid(' ', 'CHAM_GD', scal=cham19, nbret=nch)
    if (nch .ne. 0) then
        call dismoi('TYPE_SUPERVIS', cham19, 'CHAMP', repk=typcha)
        call dismoi('NOM_MAILLA', cham19, 'CHAMP', repk=noma)
        if (typcha(1:7) .eq. 'CHAM_NO') then
!       ----------------------------------
            if (ngn .ne. 0) then
                call utnono(' ', noma, 'NOEUD', nogno, noeud, &
                            iret)
                if (iret .eq. 10) then
                    call utmess('F', 'ELEMENTS_67', sk=nogno)
                else if (iret .eq. 1) then
                    valk(1) = nogno
                    valk(2) = noeud
                    call utmess('A', 'SOUSTRUC_87', nk=2, valk=valk)
                end if
            end if
            call getvtx(' ', 'NOM_CMP', scal=cmp, nbret=nc)
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
            call jeveuo(cham19//'.VALE', 'L', vr=vale)
            call focste(nomfon, cmp, vale(iddl), 'G')
            goto 10
        else if (typcha(1:9) .eq. 'CHAM_ELEM') then
!     -----------------------------------
! ---  VERIFICATION DE LA PRESENCE DES MOTS CLE GROUP_MA (OU MAILLE)
! ---  ET GROUP_NO (OU NOEUD OU POINT) DANS LE CAS D'UN CHAM_ELEM
            if (ngm .ne. 0) then
                call utnono(' ', noma, 'MAILLE', nogma, maille, &
                            iret)
                if (iret .eq. 10) then
                    call utmess('F', 'ELEMENTS_73', sk=nogma)
                else if (iret .eq. 1) then
                    valk(1) = maille
                    call utmess('A', 'UTILITAI6_72', sk=valk(1))
                end if
            end if
            if (ngn .ne. 0) then
                call utnono(' ', noma, 'NOEUD', nogno, noeud, &
                            iret)
                if (iret .eq. 10) then
                    call utmess('F', 'ELEMENTS_67', sk=nogno)
                else if (iret .eq. 1) then
                    valk(1) = nogno
                    valk(2) = noeud
                    call utmess('A', 'SOUSTRUC_87', nk=2, valk=valk)
                end if
            end if
            call dismoi('TYPE_CHAMP', cham19, 'CHAMP', repk=typch2)
            if (typch2 .eq. 'ELEM') then
                npoint = 1
                nusp = 1
                noeud = ' '
                if (maille .eq. ' ') then
                    call utmess('F', 'CHAMPS_11')
                end if
            else if (typch2 .eq. 'ELNO') then
                nusp = 1
                if (maille .eq. ' ' .or. (noeud .eq. ' ' .and. npoint .eq. 0)) then
                    call utmess('F', 'CHAMPS_12')
                end if
            else
                if (maille .eq. ' ' .or. npoint .eq. 0) then
                    call utmess('F', 'CHAMPS_13')
                end if
            end if
            call dismoi('NOM_GD', cham19, 'CHAM_ELEM', repk=nomgd)
            call dismoi('TYPE_SCA', nomgd, 'GRANDEUR', repk=type)
            if (type .ne. 'R') then
                call utmess('F', 'UTILITAI4_19')
            end if
            call utcmp1(nomgd, ' ', 1, cmp, ivari, nom_vari)
            call utch19(cham19, noma, maille, noeud, npoint, &
                        nusp, ivari, cmp, type, valr, &
                        valc, vali, iret)
            if (iret .eq. 0) then
                call focste(nomfon, cmp, valr, 'G')
            end if
            goto 10
        else
            call utmess('F', 'UTILITAI4_20', sk=typcha)
        end if
    end if
!
!     -----------------------------------------------------------------
10  continue
    call foattr(' ', 1, nomfon)
!
!     --- VERIFICATION QU'ON A BIEN CREER UNE FONCTION ---
!         ET REMISE DES ABSCISSES EN ORDRE CROISSANT
    call ordonn(nomfon, 0)
!
    call titre()
    if (niv .gt. 1) call foimpr(nomfon, niv, ifm, 0, k8b)
!
    call jedema()
end subroutine
