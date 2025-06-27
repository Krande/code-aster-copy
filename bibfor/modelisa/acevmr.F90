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
subroutine acevmr(nbocc, noma, noemax, noemaf)
    implicit none
#include "jeveux.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/getvem.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
!
    integer(kind=8) :: nbocc, noemax
    character(len=8) :: noma
!
!
! ----------------------------------------------------------------------
!     AFFE_CARA_ELEM
!     VERIFICATION DES DIMENSIONS POUR LES MASSES REPARTIES
! ----------------------------------------------------------------------
! IN  : NBOCC  : NOMBRE D'OCCURENCE
! IN  : NOMA   : NOM DU MAILLAGE
! OUT : NOEMAX : NOMBRE TOTAL DE NOEUDS MAX
! ----------------------------------------------------------------------
    character(len=24) :: magrma, manoma
    character(len=8) :: k8b
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ii, ij, in, inoe
    integer(kind=8) :: ioc, ldgm, ldnm, nb, nbgr, nbgrmx, nbv
    integer(kind=8) :: nm, nn, noema2, noemaf
    character(len=24), pointer :: group_ma(:) => null()
    integer(kind=8), pointer :: parno2(:) => null()
!-----------------------------------------------------------------------
    call jemarq()
    nbgrmx = 0
    magrma = noma//'.GROUPEMA'
    manoma = noma//'.CONNEX'
    do ioc = 1, nbocc
!        --- ON RECUPERE UNE LISTE DE GROUP_MA ---
        call getvem(noma, 'GROUP_MA', 'MASS_AJOU', 'GROUP_MA', ioc, &
                    0, k8b, nbgr)
        nbgr = -nbgr
        nbgrmx = max(nbgrmx, nbgr)
    end do
    AS_ALLOCATE(vk24=group_ma, size=nbgrmx)
    noemax = 0
    noemaf = 0
    do ioc = 1, nbocc
        noema2 = 0
        call getvem(noma, 'GROUP_MA', 'MASS_AJOU', 'GROUP_MA', ioc, &
                    0, k8b, nbgr)
        nbgr = -nbgr
        call getvem(noma, 'GROUP_MA', 'MASS_AJOU', 'GROUP_MA', ioc, &
                    nbgr, group_ma, nbv)
!
!        --- ON ECLATE LES GROUP_MA ---
        do i = 1, nbgr
            call jelira(jexnom(magrma, group_ma(i)), 'LONUTI', nb)
            call jeveuo(jexnom(magrma, group_ma(i)), 'L', ldgm)
            do in = 0, nb-1
                call jelira(jexnum(manoma, zi(ldgm+in)), 'LONMAX', nm)
                call jeveuo(jexnum(manoma, zi(ldgm+in)), 'L', ldnm)
                do nn = 1, nm
                    inoe = zi(ldnm+nn-1)
                    noema2 = max(noema2, inoe)
                end do
            end do
        end do
        noemaf = max(noemaf, noema2)
        AS_ALLOCATE(vi=parno2, size=noema2)
        do i = 1, nbgr
            call jelira(jexnom(magrma, group_ma(i)), 'LONUTI', nb)
            call jeveuo(jexnom(magrma, group_ma(i)), 'L', ldgm)
            do in = 0, nb-1
                call jelira(jexnum(manoma, zi(ldgm+in)), 'LONMAX', nm)
                call jeveuo(jexnum(manoma, zi(ldgm+in)), 'L', ldnm)
                do nn = 1, nm
                    inoe = zi(ldnm+nn-1)
                    parno2(inoe) = parno2(inoe)+1
                end do
            end do
        end do
        ii = 0
        do ij = 1, noema2
            if (parno2(ij) .eq. 0) goto 51
            ii = ii+1
51          continue
        end do
        noema2 = ii
        noemax = noemax+noema2
        AS_DEALLOCATE(vi=parno2)
    end do
    AS_DEALLOCATE(vk24=group_ma)
!
    call jedema()
end subroutine
