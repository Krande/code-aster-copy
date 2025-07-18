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

subroutine metnth(model, loadNameJv, caraElem, mateco, time, &
                  chtni, metrnl)
!
!
!
!     ARGUMENTS:
!     ----------
    implicit none
#include "jeveux.h"
#include "asterfort/calcul.h"
#include "asterfort/codent.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mecara.h"
#include "asterfort/megeom.h"
#include "asterfort/memare.h"
#include "asterfort/reajre.h"
#include "asterfort/utmess.h"
    character(len=*) :: loadNameJv, mateco
    character(len=8) :: model, caraElem
    character(len=24) :: metrnl, time, chtni
! ----------------------------------------------------------------------
!
!     CALCUL DES MATRICES ELEMENTAIRES DE CONVECTION NATURELLE
!
!     ENTREES:
!
!     LES NOMS QUI SUIVENT SONT LES PREFIXES UTILISATEUR K8:
!        MODELE : NOM DU MODELE
!        LCHAR  : OBJET CONTENANT LA LISTE DES CHARGES
!        MATE   : CHAMP DE MATERIAUX
!        CARA   : CHAMP DE CARAC_ELEM
!        TIME   : CHAMPS DE TEMPSR
!        CHTNI  : IEME ITEREE DU CHAMP DE TEMPERATURE
!        METRNL : NOM DU MATR_ELEM (N RESUELEM) PRODUIT
!
!     SORTIES:
!        METRNL  : EST REMPLI.
!
! ----------------------------------------------------------------------
!
!     VARIABLES LOCALES:
!     ------------------
!
!
    character(len=8) :: nomcha, lpain(6), lpaout(1)
    character(len=8) :: vitess
    character(len=16), parameter :: option = 'RIGI_THER_CONV'
    character(len=24) :: lchin(6), lchout(1), chgeom, chcara(18)
    character(len=24) :: chvite, ligrmo, convch
    integer(kind=8) :: iret, ilires
    integer(kind=8) :: nchar, jchar
!
! DEB-------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: ichar, iconv, jvites
!-----------------------------------------------------------------------
    call jemarq()
!     -- ON VERIFIE LA PRESENCE PARFOIS NECESSAIRE DE CARA_ELEM
    if (model(1:1) .eq. ' ') then
        call utmess('F', 'CALCULEL3_50')
    end if
!
    call jeexin(loadNameJv, iret)
    if (iret .ne. 0) then
        call jelira(loadNameJv, 'LONMAX', nchar)
        call jeveuo(loadNameJv, 'L', jchar)
    else
        nchar = 0
    end if
!
    call megeom(model, chgeom)
    call mecara(caraElem, chcara)
!
    call jeexin(metrnl, iret)
    if (iret .eq. 0) then
        metrnl = '&&METNTH           .RELR'
        call memare('V', metrnl, model(1:8), 'RIGI_THER')
    else
        call jedetr(metrnl)
    end if
!
    chvite = '????'
!
    iconv = 0
!
    lpaout(1) = 'PMATTTR'
    lchout(1) = metrnl(1:8)//'.ME000'
    do ichar = 1, nchar
        nomcha = zk24(jchar+ichar-1) (1:8)
        convch = nomcha//'.CHTH'//'.CONVE'//'.VALE'
        call jeexin(convch, iret)
        if (iret .gt. 0) then
            iconv = iconv+1
            if (iconv .gt. 1) then
                call utmess('F', 'CALCULEL3_72')
            end if
!

            call memare('V', metrnl, model(1:8), option)
!
            call jeveuo(convch, 'L', jvites)
            vitess = zk8(jvites)
            chvite = vitess
            lpain(1) = 'PGEOMER'
            lchin(1) = chgeom
            lpain(2) = 'PMATERC'
            lchin(2) = mateco
            lpain(3) = 'PCACOQU'
            lchin(3) = chcara(7)
            lpain(4) = 'PINSTR'
            lchin(4) = time
            lpain(5) = 'PVITESR'
            lchin(5) = chvite
            lpain(6) = 'PTEMPEI'
            lchin(6) = chtni
!
!
            ligrmo = model(1:8)//'.MODELE'
            ilires = 0
            ilires = ilires+1
            call codent(ilires, 'D0', lchout(1) (12:14))
            call calcul('S', option, ligrmo, 7, lchin, &
                        lpain, 1, lchout, lpaout, 'V', &
                        'OUI')
            call reajre(metrnl, lchout(1), 'V')
!
        end if
    end do
!
    call jedema()
!
end subroutine
