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

subroutine utmach(champz, ncmp, nocmp, typemz, litroz, &
                  nbtrou)
    implicit none
#include "jeveux.h"
#include "asterfort/carces.h"
#include "asterfort/celces.h"
#include "asterfort/cesexi.h"
#include "asterfort/cesred.h"
#include "asterfort/cnocns.h"
#include "asterfort/cnsred.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/int_to_char8.h"
!
    integer(kind=8) :: nbtrou, ncmp
    character(len=8) :: nocmp(*)
    character(len=*) :: champz, typemz, litroz
!
!     CETTE ROUTINE PERMET DE CREER UN OBJET JEVEUX CONTENANT UNE LISTE
!     DE NOMS OU NUMEROS DE MAILLES OU DE NOEUDS CORRESPONDANT AUX
!     MOTS-CLES TRANSMIS EN ARGUMENTS.
!
! IN  : CHAMP  : NOM D'UN CHAMP
! IN  : NCMP   : NOMBRE DE COMPOSANTES DE NOCMP
! IN  : NOCMP  : LISTE DES COMPOSANTES
! IN  : TYPEM  : PRECISE LE TYPE DE LISTE QUE L'ON VEUT RECUPERER
!              : 'NU'  : NUMEROS DE MAILLES OU DE NOEUDS
!              : 'NO'  : NOMS    DE MAILLES OU DE NOEUDS
! IN/JXOUT : LITROZ : NOM DE L'OBJET JEVEUX QUI CONTIENDRA LA LISTE DES
!                     ENTITES (MAILLE OU NOEUD) TROUVEES
! OUT : NBTROU : NOMBRE D'ENTITES TROUVEES
!     ------------------------------------------------------------------
!
    integer(kind=8) :: ierd, jcesd, jcesk, jcesl, nbent, jent, i, nbpt, nbsp, ipt
    integer(kind=8) :: isp, icp, iad, idlist, icmp, ncmpmx, gd, ier
    character(len=2) :: typem
    character(len=4) :: docu
    character(len=8) :: k8b, nomgd
    character(len=19) :: champ, chtra1, chtra2
    character(len=24) :: litrou
    character(len=24) :: valk(5)
! DEB ------------------------------------------------------------------
    call jemarq()
!
    litrou = litroz
    champ = champz
    typem = typemz
!
    nbtrou = 0
    if (ncmp .eq. 0) goto 999
!
    chtra1 = '&&UTMACH.CHAMP_COP'
    chtra2 = '&&UTMACH.CHAMP_RED'
!
    call dismoi('TYPE_CHAMP', champ, 'CHAMP', repk=docu)
    call dismoi('NUM_GD', champ, 'CHAMP', repi=gd)
    call jenuno(jexnum('&CATA.GD.NOMGD', gd), nomgd)
    if (nomgd(1:6) .eq. 'VARI_R') goto 999
!
! --- VERIFICATION QUE LES COMPOSANTES SONT DANS LE CHAMP
!
    call jelira(jexnum('&CATA.GD.NOMCMP', gd), 'LONMAX', ncmpmx)
    call jeveuo(jexnum('&CATA.GD.NOMCMP', gd), 'L', iad)
    ier = 0
    do icp = 1, ncmp
        do icmp = 1, ncmpmx
            if (nocmp(icp) .eq. zk8(iad+icmp-1)) goto 2
        end do
        ier = ier+1
        call utmess('E', 'UTILITAI5_48', sk=nocmp(icp))
2       continue
    end do
    if (ier .ne. 0) then
        call utmess('F', 'PREPOST_60')
    end if
!
!
    if (docu .eq. 'ELGA' .or. docu .eq. 'ELNO' .or. docu .eq. 'ELEM' .or. docu .eq. 'CART') then
!          ----------------
        if (docu .eq. 'CART') then
            call carces(champ, 'ELEM', k8b, 'V', chtra1, &
                        'A', ierd)
        else
            call celces(champ, 'V', chtra1)
        end if
        call cesred(chtra1, 0, [0], ncmp, nocmp, &
                    'V', chtra2)
        call jeveuo(chtra2//'.CESD', 'L', jcesd)
        call jeveuo(chtra2//'.CESK', 'L', jcesk)
        call jeveuo(chtra2//'.CESL', 'L', jcesl)
        nbent = zi(jcesd-1+1)
        call wkvect('&&UTMACH.LIST_ENT', 'V V I', nbent, jent)
        do i = 1, nbent
            nbpt = zi(jcesd-1+5+4*(i-1)+1)
            nbsp = zi(jcesd-1+5+4*(i-1)+2)
            do ipt = 1, nbpt
                do isp = 1, nbsp
                    do icp = 1, ncmp
                        call cesexi('C', jcesd, jcesl, i, ipt, &
                                    isp, icp, iad)
                        if (iad .gt. 0) then
                            zi(jent+i-1) = 1
                            goto 10
                        else
                        end if
                    end do
                end do
            end do
10          continue
        end do
        call detrsd('CHAM_ELEM_S', chtra1)
        call detrsd('CHAM_ELEM_S', chtra2)
!
!
    else if (docu .eq. 'NOEU') then
!              ----------------
        call cnocns(champ, 'V', chtra1)
        call cnsred(chtra1, 0, [0], ncmp, nocmp, &
                    'V', chtra2)
        call jeveuo(chtra2//'.CNSD', 'L', jcesd)
        call jeveuo(chtra2//'.CNSK', 'L', jcesk)
        call jeveuo(chtra2//'.CNSL', 'L', jcesl)
        nbent = zi(jcesd-1+1)
        call wkvect('&&UTMACH.LIST_ENT', 'V V I', nbent, jent)
        do i = 1, nbent
            do icp = 1, ncmp
                if (zl(jcesl-1+(i-1)*ncmp+icp)) then
                    zi(jent+i-1) = 1
                    goto 20
                end if
            end do
20          continue
        end do
        call detrsd('CHAM_NO_S', chtra1)
        call detrsd('CHAM_NO_S', chtra2)
!
    else
!
        call utmess('F', 'UTILITAI5_49', sk=docu)
!
    end if
!
    if (nbent .eq. 0) then
        valk(1) = champ
        valk(2) = nocmp(1)
        valk(3) = nocmp(2)
        valk(4) = nocmp(3)
        valk(5) = nocmp(4)
        if (docu .eq. 'NOEU') then
            call utmess('F', 'UTILITAI8_61', nk=5, valk=valk)
        else
            call utmess('F', 'UTILITAI8_62', nk=5, valk=valk)
        end if
    end if
!
    nbtrou = 0
    do i = 1, nbent
        if (zi(jent+i-1) .eq. 1) nbtrou = nbtrou+1
    end do
!
    if (typem .eq. 'NU') then
!          ---------------
        call wkvect(litrou, 'V V I', nbtrou, idlist)
        nbtrou = 0
        do i = 1, nbent
            if (zi(jent+i-1) .eq. 1) then
                nbtrou = nbtrou+1
                zi(idlist+nbtrou-1) = i
            end if
        end do
!
    else if (typem .eq. 'NO') then
!              ---------------
        call wkvect(litrou, 'V V K8', nbtrou, idlist)
        nbtrou = 0
        do i = 1, nbent
            if (zi(jent+i-1) .eq. 1) then
                nbtrou = nbtrou+1
                zk8(idlist+nbtrou-1) = int_to_char8(zi(jent+i-1))
            end if
        end do
!
    else
        call utmess('F', 'PREPOST3_6', sk=typem)
    end if
!
    call jedetr('&&UTMACH.LIST_ENT')
!
999 continue
!
    call jedema()
end subroutine
