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
subroutine nmcvci(model, hhoField, &
                  charge, infoch, fomult, numedd, depmoi, &
                  instap, cncine)
!
    use HHO_type
!
    implicit none
!
!
! BUT : CALCULER LE CHAM_NO CNCINE QUI CONTIENT  L'INCREMENT DE
!       DEPLACEMENT IMPOSE PAR LES CHARGES CINEMATIQUES.
!       POUR CELA, ON FAIT LA DIFFERENCE ENTRE LES INSTANTS "+" ET "-"
!       MAIS POUR L'INSTANT "-", IL FAUT PARTIR DU "VRAI" CHAMP
!       DE DEPLACEMENT.
!----------------------------------------------------------------------
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/ascavc.h"
#include "asterfort/assert.h"
#include "asterfort/copisd.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/vecdid_cine.h"
#include "asterfort/vtcmbl.h"
#include "asterfort/vtcreb.h"

    character(len=24), intent(in) :: model
    type(HHO_Field), intent(in) :: hhoField
    character(len=24) :: charge, infoch, fomult, numedd
    character(len=19) :: depmoi, cncine
    character(len=24) :: l2cnci(2), cncinm, cncinp, dlci
    character(len=8) :: char1, answer
    real(kind=8) :: instap, coefr(2)
    integer(kind=8) :: neq, ieq, neq2, iret, jinfc, ichar
    integer(kind=8) :: nbchar, jlchar
    character(len=1) :: typch(2)
    aster_logical :: lvcine, l_hho, l_didi
    integer(kind=8), pointer :: v_dlci(:) => null()
    real(kind=8), pointer :: cncim(:) => null()
    real(kind=8), pointer :: vale(:) => null()
!----------------------------------------------------------------------
!
    call jemarq()

    call dismoi('EXI_HHO', model, 'MODELE', repk=answer)
    l_hho = answer .eq. 'OUI'
!
!     -- CREATION DE CNCINE = 0. PARTOUT :
!     --------------------------------------
    call exisd('CHAMP_GD', cncine, iret)
    if (iret .eq. 0) then
        call vtcreb(cncine, 'V', 'R', nume_ddlz=numedd, nb_equa_outz=neq)
    end if
    call jelira(cncine(1:19)//'.VALE', 'LONMAX', ival=neq)
    call jelira(depmoi(1:19)//'.VALE', 'LONMAX', ival=neq2)
    ASSERT(neq .eq. neq2)
    call jeveuo(cncine(1:19)//'.VALE', 'E', vr=vale)
!
!
!     -- Y-A-T-IL DES CHARGES CINEMATIQUES ?
!     -----------------------------------------------------------------
    lvcine = .false.
    l_didi = .false.
    call jeveuo(infoch, 'L', jinfc)
    nbchar = zi(jinfc)
    do ichar = 1, nbchar
        if (zi(jinfc+ichar) .lt. 0) then
            lvcine = .true.
            if (zi(jinfc-1+1+ichar+3*nbchar+2) .eq. 1) then
                l_didi = .true.
            end if
        end if
    end do
!
!     -- Y-A-T-IL DES CHARGES CONTENANT DES CHARGES CINEMATIQUES ?
!     -----------------------------------------------------------------
    call jeveuo(charge, 'L', jlchar)
    call jelira(charge, 'LONMAX', ival=nbchar)
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
!     CALCUL DE UIMP+ :
!     ---------------------
    if (l_hho) then
        call ascavc(charge, infoch, fomult, numedd, instap, cncinp, dlci, &
                    l_hho, hhoField)
    else
        call ascavc(charge, infoch, fomult, numedd, instap, cncinp, dlci_=dlci)
    end if
    call jeveuo(dlci, 'L', vi=v_dlci)
!
!
!     CALCUL DE UIMP- : C'EST U- LA OU ON IMPOSE LE DEPLACEMENT
!                       ET 0. AILLEURS
!     ---------------------------------------------------------
    call copisd('CHAMP_GD', 'V', depmoi, cncinm)
    call jeveuo(cncinm(1:19)//'.VALE', 'E', vr=cncim)
    do ieq = 1, neq
        if (v_dlci(ieq) .eq. 0) then
            cncim(ieq) = 0.d0
        end if
    end do
!
!     DIFFERENCE UIMP+ - UIMP- :
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
!   Dirichlet differentiel
!
    if (l_didi) then
        call vecdid_cine(charge, infoch, numedd, cncine)
    end if
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
