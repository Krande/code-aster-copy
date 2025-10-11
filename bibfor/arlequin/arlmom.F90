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
subroutine arlmom(mailar, modarl)
!
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/nbelem.h"
#include "asterfort/nbgrel.h"
#include "asterfort/typele.h"
#include "asterfort/wkvect.h"
!
!     ARGUMENTS:
!     ----------
    character(len=8) :: mailar, modarl
!
! ----------------------------------------------------------------------
!
! ROUTINE ARLEQUIN
!
! CREATION DU PSEUDO-MODELE
!
! ----------------------------------------------------------------------
!
!
! IN  MAILAR : NOM DU PSEUDO-MAILLAGE
! IN  MODARL : NOM DU PSEUDO-MODELE
!
    character(len=24) :: modmai
    integer(kind=8) :: jmoma
    integer(kind=8) :: jdime
    integer(kind=8) :: nbma
    character(len=19) :: ligarl
    integer(kind=8) :: igrel, iel, ima, nute, nbelgr
    integer(kind=8) :: ialiel, illiel, iaux1
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
! --- INFO SUR LE MAILLAGE
!
    call jeveuo(mailar(1:8)//'.DIME', 'L', jdime)
    nbma = zi(jdime-1+3)
!
! --- CREATION DES SDs DE BASE DE MODELE
!
    call dismoi('NOM_LIGREL', modarl, 'MODELE', repk=ligarl)
    modmai = ligarl//'.TYFE'
    call wkvect(modmai, 'V V I', nbma, jmoma)
!
! --- ACCES AU LIGREL
!
    call jeveuo(ligarl//'.LIEL', 'L', ialiel)
    call jeveuo(jexatr(ligarl//'.LIEL', 'LONCUM'), 'L', illiel)
!
! --- REMPLISSAGE DE LA SD LIGREL//'.TYFE'
!
    do igrel = 1, nbgrel(ligarl)
        nute = typele(ligarl, igrel)
        nbelgr = nbelem(ligarl, igrel)
        iaux1 = ialiel-1+zi(illiel-1+igrel)-1
        do iel = 1, nbelgr
            ima = zi(iaux1+iel)
            if (ima > nbma) then
                ASSERT(.false.)
            else
                zi(jmoma+ima-1) = nute
            end if
        end do
    end do
!
    call jedetr(modmai)
    call jedema()
!
end subroutine
