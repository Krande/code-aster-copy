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

subroutine memam2(option, modele, mate, mateco, &
                  cara, compor, time, chacce, &
                  vecel, basez, ligrez)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/calcul.h"
#include "asterfort/corich.h"
#include "asterfort/detrsd.h"
#include "asterfort/gcnco2.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mecham.h"
#include "asterfort/memare.h"
#include "asterfort/reajre.h"
#include "asterfort/utmess.h"
#include "asterfort/vrcins.h"
!
    real(kind=8) :: time
    character(len=*) :: option, modele, chacce, mate, mateco, cara, vecel, basez, ligrez
!     CALCULE LES VECTEURS ELEMENTAIRES ( MASSE_MECA * CHACCE )
!
! ----------------------------------------------------------------------
! IN  : OPTION : OPTION DE CALCUL
! IN  : MODELE : NOM DU MODELE (OBLIGATOIRE)
! IN  : MATE   : CARTE DE MATERIAUX
! IN  : CARA   : CHAMP DE CARAC_ELEM
! IN  : TIME   : INSTANT DE CALCUL
! IN  : CHACCE : CHAMP D'ACCELERATION
! IN  : VECEL  : NOM DU VECT_ELEM RESULTAT
! IN  : BASEZ  : NOM DE LA BASE
! IN  : LIGREZ  : (SOUS-)LIGREL DE MODELE POUR CALCUL REDUIT
!                  SI ' ', ON PREND LE LIGREL DU MODELE
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
    character(len=1) :: base
    character(len=2) :: codret
    character(len=8) :: lpain(18), lpaout(1), newnom
    character(len=19) :: chvarc
    character(len=24) :: ligrmo, lchin(18), lchout(1), compor
    character(len=24) :: chgeom, chcara(18), chharm, vecelz
    integer(kind=8) :: icode, iret, nh
    character(len=24), pointer :: rerr(:) => null()
!-----------------------------------------------------------------------
    data chvarc/'&&MEMAM2.VARC'/
    call jemarq()
    newnom = '.0000000'
    vecelz = vecel
    base = basez
    if (modele(1:1) .eq. ' ') then
        call utmess('F', 'CALCULEL2_82')
    end if
!
    nh = 0
    call mecham('MASS_MECA', modele, cara, nh, chgeom, &
                chcara, chharm, icode)
    call vrcins(modele, mate, ' ', time, chvarc, &
                codret)
!
    call memare(base, vecel, modele, option)
    call jeveuo(vecelz(1:19)//'.RERR', 'E', vk24=rerr)
    rerr(3) (1:3) = 'OUI'
!
    call jeexin(vecelz(1:19)//'.RELR', iret)
    if (iret .gt. 0) call jedetr(vecelz(1:19)//'.RELR')
    if (icode .eq. 1) goto 10
!
    ligrmo = ligrez
    if (ligrmo .eq. ' ') ligrmo = modele(1:8)//'.MODELE'
!
    lpaout(1) = 'PVECTUR'
!
    lpain(1) = 'PGEOMER'
    lchin(1) = chgeom
    lpain(2) = 'PMATERC'
    lchin(2) = mateco
    lpain(3) = 'PVARCPR'
    lchin(3) = chvarc
    lpain(4) = 'PCAGNPO'
    lchin(4) = chcara(6)
    lpain(5) = 'PCAGEPO'
    lchin(5) = chcara(5)
    lpain(6) = 'PCACOQU'
    lchin(6) = chcara(7)
    lpain(7) = 'PCADISM'
    lchin(7) = chcara(3)
    lpain(8) = 'PCAORIE'
    lchin(8) = chcara(1)
    lpain(9) = 'PCASECT'
    lchin(9) = chcara(8)
    lpain(10) = 'PCAARPO'
    lchin(10) = chcara(9)
    lpain(11) = 'PCACABL'
    lchin(11) = chcara(10)
    lpain(12) = 'PCAGNBA'
    lchin(12) = chcara(11)
    lpain(13) = 'PCAPOUF'
    lchin(13) = chcara(13)
    lpain(14) = 'PACCELR'
    lchin(14) = chacce
    lpain(15) = 'PNBSP_I'
    lchin(15) = chcara(16)
    lpain(16) = 'PFIBRES'
    lchin(16) = chcara(17)
    lpain(17) = 'PCOMPOR'
    lchin(17) = compor
    lpain(18) = 'PCINFDI'
    lchin(18) = chcara(15)
!
    lchout(1) = '&&MEMAM2.???????'
    call gcnco2(newnom)
    lchout(1) (10:16) = newnom(2:8)
    call corich('E', lchout(1), ichin_=-1)
    call calcul('S', option, ligrmo, 18, lchin, &
                lpain, 1, lchout, lpaout, base, &
                'OUI')
!
    call reajre(vecelz, lchout(1), base)
!
10  continue
    call detrsd('CHAMP_GD', chvarc)
!
    call jedema()
end subroutine
