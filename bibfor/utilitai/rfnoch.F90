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

subroutine rfnoch()
    implicit none
!
!     OPERATEUR "RECU_FONCTION"   MOT CLE "NOEUD_CHOC"
!     ------------------------------------------------------------------
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterfort/dismoi.h"
#include "asterfort/foattr.h"
#include "asterfort/focrch.h"
#include "asterfort/foimpr.h"
#include "asterfort/getvid.h"
#include "asterfort/getvtx.h"
#include "asterfort/infmaj.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/ordonn.h"
#include "asterfort/titre.h"
#include "asterfort/utmess.h"
#include "asterfort/utnono.h"
    integer(kind=8) :: ifm, niv
    integer(kind=8) :: n, nc, ng, int, ind, nsst, iret
    character(len=8) :: k8b, noma, sst, basemo, noeud
    character(len=24) :: valk(2), nogno, intitu
    character(len=16) :: parax, paray, nomcmd, typcon
    character(len=19) :: listr, nomfon, resu
!     ------------------------------------------------------------------
    call jemarq()
!
    call infmaj()
    call infniv(ifm, niv)
!
    call getres(nomfon, typcon, nomcmd)
!
    call getvtx(' ', 'INTITULE', scal=intitu, nbret=int)
    call getvid(' ', 'RESU_GENE', scal=resu, nbret=n)
    call getvtx(' ', 'PARA_X', scal=parax, nbret=n)
    call getvtx(' ', 'PARA_Y', scal=paray, nbret=n)
    call getvid(' ', 'LIST_PARA', scal=listr, nbret=ind)
    call getvtx(' ', 'SOUS_STRUC', scal=sst, nbret=nsst)
!
    call getvtx(' ', 'NOEUD_CHOC', scal=noeud, nbret=nc)
    call getvtx(' ', 'GROUP_NO_CHOC', scal=nogno, nbret=ng)
!
    if (nc .ne. 0) then
!
        call focrch(nomfon, resu, noeud, parax, paray, &
                    'G', int, intitu, ind, listr, &
                    sst, nsst, iret)
!
    else
!
        call dismoi('BASE_MODALE', resu, 'RESU_DYNA', repk=basemo)
        call dismoi('NOM_MAILLA', basemo, 'RESULTAT', repk=noma)
!
        call utnono(' ', noma, 'NOEUD', nogno, noeud, &
                    iret)
        if (iret .eq. 10) then
            call utmess('F', 'ELEMENTS_67', sk=nogno)
        else if (iret .eq. 1) then
            valk(1) = nogno
            valk(2) = noeud
            call utmess('A', 'SOUSTRUC_87', nk=2, valk=valk)
        end if
!
        call focrch(nomfon, resu, noeud, parax, paray, &
                    'G', int, intitu, ind, listr, &
                    sst, nsst, iret)
    end if
!
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
