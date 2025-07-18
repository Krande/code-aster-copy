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

subroutine fimped(modele, mateco, numedd, neq, vitini, &
                  vitent, veccor, veanec, vaanec, temps, &
                  foimpe)
    implicit none
#include "jeveux.h"
#include "asterfort/asasve.h"
#include "asterfort/calcul.h"
#include "asterfort/corich.h"
#include "asterfort/detrsd.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mecact.h"
#include "asterfort/reajre.h"
    integer(kind=8) :: i, jvaanc, neq, npain
    character(len=8) :: lpain(5), lpaout(1)
    character(len=24) :: modele, mateco, numedd, vitini, veccor
    character(len=24) :: vitent, chinst
    character(len=24) :: veanec, vaanec, lchin(5), lchout(1)
    character(len=24) :: chgeom, ligrel
    real(kind=8) :: foimpe(neq), temps
    character(len=8), pointer :: lgrf(:) => null()
    real(kind=8), pointer :: vale(:) => null()
!
!-----------------------------------------------------------------------
    call jemarq()
!
    chinst = '&&CHINST'
    call mecact('V', chinst, 'MODELE', modele(1:8)//'.MODELE', 'INST_R', &
                ncmp=1, nomcmp='INST', sr=temps)
    call jedetr(veanec(1:19)//'.RELR')
!
    ligrel = modele(1:8)//'.MODELE'
    call jeveuo(ligrel(1:19)//'.LGRF', 'L', vk8=lgrf)
    chgeom = lgrf(1)//'.COORDO'
    lpain(1) = 'PGEOMER'
    lchin(1) = chgeom
    lpain(2) = 'PMATERC'
    lchin(2) = mateco
!
    lpain(3) = 'PVITPLU'
    lchin(3) = vitini
    lpain(4) = 'PVITENT'
    lchin(4) = vitent
!
    lpain(5) = 'PINSTR'
    lchin(5) = chinst
!
    npain = 5
    lpaout(1) = 'PVECTUR'
    lchout(1) = veccor
!
    call calcul('S', 'IMPE_ABSO', ligrel, npain, lchin, &
                lpain, 1, lchout, lpaout, 'V', &
                'OUI')
!
    call corich('E', lchout(1), ichin_=-1)
!
    call reajre(veanec, lchout(1), 'V')
    call asasve(veanec, numedd, 'R', vaanec)
!
    call jeveuo(vaanec, 'L', jvaanc)
    call jeveuo(zk24(jvaanc) (1:19)//'.VALE', 'L', vr=vale)
!
    do i = 1, neq
        foimpe(i) = vale(i)
    end do
    call detrsd('CHAMP_GD', zk24(jvaanc) (1:19))
!
!
    call jedema()
end subroutine
