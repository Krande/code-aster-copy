! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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
subroutine meonme(model, nbLoad, listLoadName, mate, mateco, matrElem)
!
implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/calcul.h"
#include "asterfort/codent.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/mecham.h"
#include "asterfort/memare.h"
#include "asterfort/reajre.h"
#include "asterfort/utmess.h"
!
integer, intent(in) :: nbLoad
character(len=8), intent(in) :: model, listLoadName(*)
character(len=19), intent(in) :: matrElem
character(len=*), intent(in) :: mate, mateco
!
! --------------------------------------------------------------------------------------------------
!
!     CALCUL DES MATRICES ELEMENTAIRES D 'IMPEDANCE PAR ONDE INCIDENTE
!     ACOUSTIQUE DANS LE PHENOMENE MECANIQUE
!
! --------------------------------------------------------------------------------------------------
!
! IN  : MODELE : NOM DU MODELE
! IN  : NCHAR  : NOMBRE DE CHARGES
! IN  : LCHAR  : LISTE DES CHARGES
! IN  : MATE   : CARTE DE MATERIAU
! VAR : MATEL  : NOM  DU  MATELE (N RESUELEM) PRODUIT
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: option = 'ONDE_FLUI'
    integer :: nh
    character(len=8) :: caraElem, lpain(3), lpaout(1)
    character(len=24) :: ligrmo, lchin(3), lchout(1)
    character(len=24) :: chgeom, chcara(18), chharm
    integer :: iLoad, icode, ilires, iret
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
    ASSERT(model(1:1) .ne. ' ')
    caraElem = ' '
    nh       = 0
    call mecham(option, model, caraElem, nh, chgeom,&
                chcara, chharm, icode)
!
    call jeexin(matrElem//'.RERR', iret)
    if (iret .gt. 0) then
        call jedetr(matrElem//'.RERR')
        call jedetr(matrElem//'.RELR')
    endif
    call memare('G', matrElem, model, mate, ' ', option)
!
    lpaout(1) = 'PMATUUR'
    lchout(1) = matrElem(1:8)//'.ME001'
    ilires = 0
    if (listLoadName(1) .ne. ' ') then
        ligrmo = model//'.MODELE'
        lpain(1) = 'PGEOMER'
        lchin(1) = chgeom
        lpain(3) = 'PMATERC'
        lchin(3) = mateco
        do iLoad = 1, nbLoad
            call exisd('CHAMP_GD', listLoadName(iLoad)//'.CHME.ONDE ', iret)
            if (iret .ne. 0) then
                lpain(2) = 'PONDECR'
                lchin(2) = listLoadName(iLoad)//'.CHME.ONDE .DESC'
                ilires = ilires + 1
                call codent(ilires, 'D0', lchout(1) (12:14))
                call calcul('S', option, ligrmo, 3, lchin,&
                            lpain, 1, lchout, lpaout, 'G',&
                            'OUI')
                call reajre(matrElem, lchout(1), 'G')
            endif
        end do
    endif
!
    if (ilires .eq. 0) then
        call utmess('F', 'CHARGES6_84')
    endif
    call jedema()
end subroutine
