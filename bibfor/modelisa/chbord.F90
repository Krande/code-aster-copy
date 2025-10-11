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
subroutine chbord(nomo, nbCell, listCellNume, mabord, nbmapr, &
                  nbmabo)
    implicit none
#include "jeveux.h"
#include "asterfort/dismoi.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/teattr.h"
!
    integer(kind=8), pointer :: listCellNume(:)
    integer(kind=8) :: nbCell, mabord(*), nbmapr, nbmabo
    character(len=*) :: nomo
!
!      OPERATEURS :     AFFE_CHAR_MECA ET AFFE_CHAR_MECA_C
!                                      ET AFFE_CHAR_MECA_F
!
! IN  : NOMO   : NOM DU MODELE
! IN  : nbCell : NOMBRE DE MAILLES
! IN  : LISTMA : LISTE DES NUMEROS DE MAILLE
! OUT : MABORD : MABORD(IMA) = 0 , MAILLE "PRINCIPAL"
!                MABORD(IMA) = 1 , MAILLE "BORD"
! OUT : NBMAPR : NOMBRE DE MAILLES "PRINCIPAL"
! OUT : NBMABO : NOMBRE DE MAILLES "BORD"
!
!-----------------------------------------------------------------------
    integer(kind=8) :: iret, nbgrel, igrel, ialiel, nel, itypel, iCell, ier, cellNume, iel
    integer(kind=8) :: traite
    character(len=8) :: modele, dmo, dma
    character(len=16) :: nomte
    character(len=19) :: nolig
!     ------------------------------------------------------------------
!
    modele = nomo
    nbmapr = 0
    nbmabo = 0
    call dismoi('NOM_LIGREL', modele, 'MODELE', repk=nolig)
!
    call jeexin(nolig//'.LIEL', iret)
    if (iret .eq. 0) goto 999
!
    call jelira(nolig//'.LIEL', 'NUTIOC', nbgrel)
    if (nbgrel .le. 0) goto 999
!
    traite = 0
    do igrel = 1, nbgrel
        call jeveuo(jexnum(nolig//'.LIEL', igrel), 'L', ialiel)
        call jelira(jexnum(nolig//'.LIEL', igrel), 'LONMAX', nel)
        itypel = zi(ialiel-1+nel)
        call jenuno(jexnum('&CATA.TE.NOMTE', itypel), nomte)
        do iCell = 1, nbCell
            cellNume = listCellNume(iCell)
            do iel = 1, nel-1
                if (cellNume .eq. zi(ialiel-1+iel)) then
                    traite = traite+1
                    call teattr('S', 'DIM_TOPO_MODELI', dmo, ier, typel=nomte)
                    call teattr('S', 'DIM_TOPO_MAILLE', dma, ier, typel=nomte)
                    if (dmo .eq. dma) then
!                    on a un element principal
                        mabord(iCell) = 0
                        nbmapr = nbmapr+1
                    else
!                    on a un element de bord
                        mabord(iCell) = 1
                        nbmabo = nbmabo+1
                    end if
                    goto 20
                end if
            end do
20          continue
        end do
        if (traite .eq. nbCell) goto 999
    end do
!
999 continue
!
end subroutine
