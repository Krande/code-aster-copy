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

subroutine irmad0(ifc, versio, nstat, chamno, nomsym)
    implicit none
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterfort/dismoi.h"
#include "asterfort/irmad1.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/nbec.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "asterfort/int_to_char8.h"
!
    integer(kind=8) :: versio, nstat
    character(len=*) :: chamno(*), nomsym
!     IMPRESSION D'UNE LISTE DE CHAMNO A COMPOSANTES REELLES OU
!     COMPLEXES AU FORMAT IDEAS  ( COMMANDE MACRO_MADMACS )
!     ON IMPRIME LE CHAMNO A L'AIDE D'UN DATASET 252
!C
!     ------------------------------------------------------------------
!
    character(len=1) :: type, typi
    integer(kind=8) :: gd, gdi
    character(len=8) :: k8b, nomma, nomgd
    character(len=16) :: nomcmd
    character(len=19) :: chamn
    character(len=24) :: nomnu
!     ------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iad, iaec, iaprno
    integer(kind=8) :: ibid, ifc, ino, iret, itype
    integer(kind=8) ::  nbno, ncmpmx, nec
    character(len=8), pointer :: nomnoe(:) => null()
    integer(kind=8), pointer :: numnoe(:) => null()
    integer(kind=8), pointer :: nueq(:) => null()
!-----------------------------------------------------------------------
    call jemarq()
    chamn = chamno(1)
!
!
!     --- NOM DU MAILLAGE
    call dismoi("NOM_MAILLA", chamn, "CHAM_NO", repk=nomma)
!
!     --- NOM DU PROFIL AUX NOEUDS ASSOCIE S'IL EXISTE
    call dismoi("NUME_EQUA", chamn, "CHAM_NO", repk=nomnu)
!
    call jelira(chamn//'.VALE', 'TYPE', cval=type)
    if (type .eq. 'R') then
        itype = 1
    else if (type .eq. 'C') then
        itype = 2
    else
        call getres(k8b, k8b, nomcmd)
        call utmess('A', 'PREPOST_97', sk=type(1:1))
        goto 999
    end if
!
    call dismoi("NUM_GD", chamn, "CHAM_NO", repi=gd)
    call jenuno(jexnum('&CATA.GD.NOMGD', gd), nomgd)
!
!     --- ON VERIFIE QUE TOUS LES CHAMPS SONT IDENTIQUES
!
    do i = 1, nstat
        chamn = chamno(i)
        call dismoi("NUM_GD", chamn, "CHAM_NO", repi=gdi)
        if (gdi .ne. gd) then
            call utmess('F', 'PREPOST2_67')
        end if
        call jelira(chamn//'.VALE', 'TYPE', cval=typi)
        if (typi .ne. type) then
            call utmess('F', 'PREPOST2_69')
        end if
    end do
!
!     --- NOMBRE D'ENTIERS CODES POUR LA GRANDEUR NOMGD
    nec = nbec(gd)
    call jeexin('&&IRMAD0.ENT_COD', iret)
    if (iret .ne. 0) call jedetr('&&IRMAD0.ENT_COD')
    call wkvect('&&IRMAD0.ENT_COD', 'V V I', nec, iaec)
    call jelira(jexnum('&CATA.GD.NOMCMP', gd), 'LONMAX', ncmpmx)
    call jeveuo(jexnum('&CATA.GD.NOMCMP', gd), 'L', iad)
!
!     --- SI LE CHAMP EST A REPRESENTATION CONSTANTE: RIEN DE SPECIAL
!
!     --- SI LE CHAMP EST DECRIT PAR UN "PRNO":
    call jeveuo(nomnu(1:19)//'.NUEQ', 'L', vi=nueq)
    call jenonu(jexnom(nomnu(1:19)//'.LILI', '&MAILLA'), ibid)
    call jeveuo(jexnum(nomnu(1:19)//'.PRNO', ibid), 'L', iaprno)
!
!     --- NOMBRE DE NOEUDS DU MAILLAGE: NBNO
    call dismoi('NB_NO_MAILLA', nomma, 'MAILLAGE', repi=nbno)
!
!     --- CREATION LISTES DES NOMS ET DES NUMEROS DES NOEUDS A IMPRIMER
    AS_ALLOCATE(vk8=nomnoe, size=nbno)
    AS_ALLOCATE(vi=numnoe, size=nbno)
    do ino = 1, nbno
        nomnoe(ino) = int_to_char8(ino)
        numnoe(ino) = ino
    end do
!
    ncmpmx = 6
    call irmad1(ifc, versio, nbno, zi(iaprno), nueq, &
                nec, zi(iaec), ncmpmx, itype, nstat, &
                chamno, zk8(iad), nomsym, numnoe)
!
! --- MENAGE
    call jedetr('&&IRMAD0.ENT_COD')
    AS_DEALLOCATE(vk8=nomnoe)
    AS_DEALLOCATE(vi=numnoe)
!
999 continue
    call jedema()
end subroutine
