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
!
subroutine recmst(graexc, grdmod, nnoeex, ilnoex, ilcpex, &
                  nmost1, modsta)
!    C. DUVAL
!-----------------------------------------------------------------------
!  BUT: RECUPERER LES INFORMATIONS DE TYPE MODE STATIQUE POUR
    implicit none
!        LE CALCUL DYNAMIQUE ALEATOIRE
!
!-----------------------------------------------------------------------
!
! GRAEXC   /IN / : GRANDEUR EXCITATION
! GRDMOD   /IN / : GRANDEUR A RECUPERE DANS LES MODES DYN ET STAT
! NNOEEX   /IN / : NOMBRE DE NOEUDS EXCITATION
! ILNOEX   /IN / : POINTEUR DANS ZK8 SUR LES NOEUDS EXCITATION
! ILCPEX   /IN / : POINTEUR DANS ZK8 SUR LES DDLS EXCITATION
! NMOST1   /OUT/ : NOMBRE D OCCURENCE DU MOT CLE MODE_STAT
! MODSTA   /OUT/ : CONCEPT MODE_STAT
!                     AUX MODES STATIQUES EN DEPLACEMENT
!                     AUX MODES STATIQUES EN GRANDEUR REPONSE (CALCUL)
!
!
!
!
!
#include "jeveux.h"
#include "asterfort/getvid.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jeveut.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/rsadpa.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    character(len=8) :: kbid, modsta
    character(len=16) :: graexc, grdmod
    character(len=24) :: k24bd1, k24bd2, k24bd3
!
!---------MODES STATIQUES
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i1, i2, i3, i3b, i4, ibid, ilamsc
    integer(kind=8) :: ilamst, ilcpex, ilnoex, ilorms, jpara, n
    integer(kind=8) :: nmost1, nmost2, nnoeex
    character(len=8) :: c_nume_noeud
    integer(kind=8), pointer :: ordr(:) => null()
!-----------------------------------------------------------------------
    call jemarq()
    call getvid(' ', 'MODE_STAT', nbval=0, nbret=n)
    nmost1 = -n
    if (n .ne. 0) then
        call getvid(' ', 'MODE_STAT', scal=modsta, nbret=ibid)
    end if
!
!---------CONSTITUTION DE LA LISTE DES ADRESSES DES MODES STATIQUES
!
    if (graexc .ne. 'DEPL_R') goto 999
!
    k24bd1 = modsta//'           .NOEU'
!     CALL JELIRA(K24BD1,'LONMAX',NMOST2,K24BD2)
!     CALL JEVEUO(K24BD1,'L',IAD1)
    call jelira(modsta//'           .ORDR', 'LONMAX', nmost2)
    call jeveuo(modsta//'           .ORDR', 'L', vi=ordr)
    call wkvect('&&OP0131.LISTADORMOSTA', 'V V I', nnoeex, ilorms)
    do i1 = 1, nnoeex
        zi(ilorms+i1-1) = 0
        c_nume_noeud(1:8) = ' '
        if (zk8(ilnoex+i1-1) (1:1) .eq. 'N') then
            c_nume_noeud(1:7) = zk8(ilnoex+i1-1) (2:8)
            c_nume_noeud(8:8) = ' '
        end if
        do i2 = 1, nmost2
            call rsadpa(modsta, 'L', 1, 'NOEUD_CMP', ordr(i2), &
                        0, sjv=jpara, styp=kbid)
            if ((zk8(ilnoex+i1-1)//zk8(ilcpex+i1-1)) .eq. zk16(jpara)) then
                zi(ilorms+i1-1) = i2
            elseif ((c_nume_noeud//zk8(ilcpex+i1-1)) .eq. zk16(jpara)) then
                zi(ilorms+i1-1) = i2
            end if
        end do
    end do
!
!---------CONSTITUTION DES ADRESSES DES MODES STATIQUES
!
    call wkvect('&&OP0131.LISTADRMODSTA', 'V V I', nnoeex, ilamst)
    call wkvect('&&OP0131.LISTADRMODSTAC', 'V V I', nnoeex, ilamsc)
    do i1 = 1, nnoeex
        i2 = zi(ilorms+i1-1)
        if (i2 .eq. 0) then
            call utmess('F', 'ALGORITH10_33')
        end if
        k24bd1 = modsta//'           .TACH'
        call jenonu(jexnom(k24bd1(1:19)//'.DESC', 'DEPL'), ibid)
        call jeveuo(jexnum(k24bd1, ibid), 'L', i3)
        call jenonu(jexnom(k24bd1(1:19)//'.DESC', grdmod), ibid)
        call jeveuo(jexnum(k24bd1, ibid), 'L', i3b)
        k24bd2 = zk24(i3+i2-1)
        k24bd3 = k24bd2(1:19)//'.VALE'
        call jeveut(k24bd3, 'L', i4)
        zi(ilamst+i1-1) = i4
        k24bd2 = zk24(i3b+i2-1)
        k24bd3 = k24bd2(1:19)//'.VALE'
        call jeveut(k24bd3, 'L', i4)
        zi(ilamsc+i1-1) = i4
    end do
999 continue
    call jedema()
end subroutine
