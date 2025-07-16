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
subroutine nxcvci(loadNameJv, loadInfoJv, loadFuncJv, numedd, tempmoi, &
                  instap, cncine)
!
!
    implicit none
!
!
! BUT : CALCULER LE CHAM_NO CNCINE QUI CONTIENT  L'INCREMENT DE
!       TEMPERATURE IMPOSE PAR LES CHARGES CINEMATIQUES.
!       POUR CELA, ON FAIT LA DIFFERENCE ENTRE LES INSTANTS "+" ET "-"
!       MAIS POUR L'INSTANT "-", IL FAUT PARTIR DU "VRAI" CHAMP
!       DE TEMPERATURE.
!----------------------------------------------------------------------
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/ascavc.h"
#include "asterfort/assert.h"
#include "asterfort/copisd.h"
#include "asterfort/detrsd.h"
#include "asterfort/exisd.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/vtcmbl.h"
#include "asterfort/vtcreb.h"

    character(len=24), intent(in) :: loadNameJv, loadInfoJv, loadFuncJv
    character(len=24) :: numedd
    character(len=19) :: tempmoi, cncine
    character(len=24) :: l2cnci(2), cncinm, cncinp, dlci
    character(len=8) :: char1
    real(kind=8) :: instap, coefr(2)
    integer(kind=8) :: neq, ieq, neq2, iret, jinfc, ichar
    integer(kind=8) :: nbchar, jlchar
    character(len=1) :: typch(2)
    aster_logical :: lvcine
    integer(kind=8), pointer :: v_dlci(:) => null()
    real(kind=8), pointer :: cncim(:) => null()
    real(kind=8), pointer :: vale(:) => null()
!----------------------------------------------------------------------
!
    call jemarq()

    if (cncine .eq. ' ') cncine = '&&ASCAVC.VCI'
!
!     -- CREATION DE CNCINE = 0. PARTOUT :
!     --------------------------------------
    call exisd('CHAMP_GD', cncine, iret)
    if (iret .eq. 0) then
        call vtcreb(cncine, 'V', 'R', nume_ddlz=numedd, nb_equa_outz=neq)
    end if
    call jelira(cncine(1:19)//'.VALE', 'LONMAX', ival=neq)
    call jelira(tempmoi(1:19)//'.VALE', 'LONMAX', ival=neq2)
    ASSERT(neq .eq. neq2)
    call jeveuo(cncine(1:19)//'.VALE', 'E', vr=vale)
!
!
!     -- Y-A-T-IL DES CHARGES CINEMATIQUES ?
!     -----------------------------------------------------------------
    lvcine = .false.
    call jeveuo(loadInfoJv, 'L', jinfc)
    do ichar = 1, zi(jinfc)
        if (zi(jinfc+ichar) .lt. 0) lvcine = .true.
    end do
!
!     -- Y-A-T-IL DES CHARGES CONTENANT DES CHARGES CINEMATIQUES ?
!     -----------------------------------------------------------------
    call jeveuo(loadNameJv, 'L', jlchar)
    call jelira(loadNameJv, 'LONMAX', ival=nbchar)
    do ichar = 1, nbchar
        char1 = zk24(jlchar-1+ichar) (1:8)
    end do
!
!     -- S'IL N'Y A PAS DE CHARGES CINEMATIQUES, IL N'Y A RIEN A FAIRE:
!     -----------------------------------------------------------------
    if (.not. lvcine) goto 999
!
!
!     -- S'IL Y A DES CHARGES CINEMATIQUES :
!     -----------------------------------------------------------------
    cncinm = '&&NMCHAR.CNCIMM'
    cncinp = '&&NMCHAR.CNCIMP'
    dlci = '&&NMCHAR.DLCI'
!
!
!     CALCUL DE TIMP+ :
!     ---------------------
    call ascavc(loadNameJv, loadInfoJv, loadFuncJv, numedd, instap, cncinp, dlci_=dlci)
    call jeveuo(dlci, 'L', vi=v_dlci)
!
!
!     CALCUL DE TIMP- : C'EST U- LA OU ON IMPOSE LE DEPLACEMENT
!                       ET 0. AILLEURS
!     ---------------------------------------------------------
    call copisd('CHAMP_GD', 'V', tempmoi, cncinm)
    call jeveuo(cncinm(1:19)//'.VALE', 'E', vr=cncim)
    do ieq = 1, neq
        if (v_dlci(ieq) .eq. 0) then
            cncim(ieq) = 0.d0
        end if
    end do
!
!     DIFFERENCE TIMP+ - TIMP- :
!     ---------------------------
    coefr(1) = -1.d0
    coefr(2) = +1.d0
    l2cnci(1) = cncinm
    l2cnci(2) = cncinp
    typch(1) = 'R'
    typch(2) = 'R'
    call vtcmbl(2, typch, coefr, typch, l2cnci, &
                typch(1), cncine)
!
!     MENAGE :
!     ---------
    call detrsd('CHAM_NO', cncinm)
    call detrsd('CHAM_NO', cncinp)
    call jedetr(dlci)
!
999 continue
    call jedema()
!
end subroutine
