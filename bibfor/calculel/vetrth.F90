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

subroutine vetrth(model, loadNameJv, loadInfoJv, caraElem, mateco, &
                  inst, chtn, chti, chlapm, chlapp, &
                  veres)
!
!
    implicit none
#include "jeveux.h"
#include "asterfort/calcul.h"
#include "asterfort/corich.h"
#include "asterfort/dismoi.h"
#include "asterfort/gcnco2.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mecara.h"
#include "asterfort/megeom.h"
#include "asterfort/reajre.h"
#include "asterfort/utmess.h"
    character(len=8), intent(in) :: model, caraElem
    character(len=24) :: loadNameJv, loadInfoJv, inst, chtn, chti
    character(len=24) :: chlapm, chlapp, veres, mateco
! ----------------------------------------------------------------------
! CALCUL DES VECTEURS ELEMENTAIRES - SECOND MEMBRE DU PROBLEME TRANSPORT
!                                    EN THERMIQUE NON LINEAIRE -
!
! IN  MODELE  : NOM DU MODELE
! IN  CHARGE  : LISTE DES CHARGES
! IN  INFCHA  : INFORMATIONS SUR LES CHARGES
! IN  CARELE  : CHAMP DE CARA_ELEM
! IN  MATE    : CHAMP DE MATERIAU
! IN  INST    : CARTE CONTENANT LA VALEUR DU TEMPS ET AUTRES PARAMETRES
! IN  CHTN    : CHAMP DE TEMPERATURE A L'INSTANT PRECEDENT
! IN  CHTI    : I EME ITERE DU CHAMP DE TEMPERATURE
! IN  CHLAPM  : I EME ITERE DU CHAMP DES FONCTIONS ENTHALPIE
! OUT CHLAPP  : I+1 EME ITERE DU CHAMP DES FONCTIONS ENTHALPIE
! OUT VERES   : VECTEURS ELEMENTAIRES (SECOND MEMBRE)
!
!
!
    character(len=8) :: nomcha, lpain(7), lpaout(4), newnom
    character(len=16) :: option
    character(len=24) :: ligrmo, lchin(7), lchout(4)
    character(len=24) :: chvite, convch, chgeom, chcara(18)
    integer(kind=8) :: iret, jvites
    integer(kind=8) :: icha, ichar, iconv, jchar, jinf, nchar
!
!-----------------------------------------------------------------------
    call jemarq()
    newnom = '.0000000'
    call dismoi('NOM_LIGREL', model, 'MODELE', repk=ligrmo)
    call jeexin(loadNameJv, iret)
    if (iret .ne. 0) then
        call jelira(loadNameJv, 'LONMAX', nchar)
        call jeveuo(loadNameJv, 'L', jchar)
        call jeveuo(loadInfoJv, 'L', jinf)
    else
        nchar = 0
    end if
!
    call megeom(model, chgeom)
    call mecara(caraElem, chcara)
!
    lpaout(1) = 'PVECTTR'
    lpaout(2) = 'PLAGRP '
    lpaout(3) = 'PRESIDU'
    lpain(1) = 'PGEOMER'
    lchin(1) = chgeom
    lpain(2) = 'PMATERC'
    lchin(2) = mateco
    lpain(3) = 'PINSTR'
    lchin(3) = inst
    lpain(4) = 'PTEMPER'
    lchin(4) = chtn
    lpain(5) = 'PTEMPEI'
    lchin(5) = chti
    lpain(6) = 'PLAGRM '
    lchin(6) = chlapm
!
! --- TERME VOLUMIQUE PROVENANT DU COMPORTEMENT
!
    iconv = 0
    chvite = ' '
    option = 'CHAR_THER_TNL'
    lchout(2) = chlapp
    lpain(7) = 'PVITESR'
!
    do ichar = 1, nchar
        nomcha = zk24(jchar+ichar-1) (1:8)
        convch = nomcha//'.CHTH'//'.CONVE'//'.VALE'
        call jeexin(convch, iret)
        if (iret .gt. 0) then
            iconv = iconv+1
            if (iconv .gt. 1) then
                call utmess('F', 'CALCULEL3_72')
            end if
            call jeveuo(convch, 'L', jvites)
            chvite = zk8(jvites)
        end if
    end do
    if (iconv .eq. 0) then
        call utmess('F', 'CALCULEL5_38')
    end if
    lchin(7) = chvite
!
    call gcnco2(newnom)
    lchout(1) = '&&VETRTH.'//newnom(2:8)
    call corich('E', lchout(1), ichin_=-1)
    call gcnco2(newnom)
    lchout(3) = '&&VETRTH.'//newnom(2:8)
    call corich('E', lchout(3), ichin_=-1)
    call calcul('S', option, ligrmo, 7, lchin, &
                lpain, 3, lchout, lpaout, 'V', &
                'OUI')
    call reajre(veres, lchout(3), 'V')
!
! --- TERME SURFACIQUE PROVENANT DES CONDITIONS AUX LIMITES
!
    lpaout(1) = 'PRESIDU'
    lpain(2) = 'PFLUXNL'
    if (nchar .gt. 0) then
        do icha = 1, nchar
            lchin(2) = zk24(jchar+icha-1) (1:8)//'.CHTH.FLUNL.DESC'
            call jeexin(lchin(2), iret)
            if (iret .ne. 0) then
                option = 'CHAR_THER_FLUTNL'
                call gcnco2(newnom)
                lchout(1) = '&&VETRTH.'//newnom(2:8)
                call corich('E', lchout(1), ichin_=-1)
                call calcul('S', option, ligrmo, 5, lchin, &
                            lpain, 1, lchout, lpaout, 'V', &
                            'OUI')
                call reajre(veres, lchout(1), 'V')
            end if
        end do
    end if
!
    call jedema()
end subroutine
