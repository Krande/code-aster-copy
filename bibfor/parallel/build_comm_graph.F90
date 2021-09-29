! --------------------------------------------------------------------
! Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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

subroutine build_comm_graph(noma, nomgr, base)
    implicit none
#include "asterf_config.h"
#include "asterf.h"
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetc.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/wkvect.h"
#include "asterfort/asmpi_info.h"
#include "asterc/asmpi_bcast_i.h"
#include "asterc/asmpi_comm.h"
    character(len=8) :: noma
    character(len=*) :: nomgr
    character(len=1) :: base

#ifdef ASTER_HAVE_MPI
!
    integer :: rang, nbproc, jdojoi, nbjoin, iaux, jgraco
    integer :: jmasqu, iproc1, iproc2, nbedge, posit, nmatch
    integer :: jtmp, jordjo, num, iret, jnbjoi
    integer(kind=4) :: iaux4, n4e
    mpi_int :: mrank, msize, mpicou
!
    character(len=8) :: k8bid
!
    call jemarq()
!
    call asmpi_comm('GET', mpicou)
    call asmpi_info(rank = mrank, size = msize)
    rang = to_aster_int(mrank)
    nbproc = to_aster_int(msize)
!
    call wkvect('&&BUCOGR.GRAPH_COMM', 'V V I', nbproc*nbproc, jgraco)
!
    nbjoin = 0
    jdojoi = 0
    call jeexin(noma//'.DOMJOINTS', iret)
    if(iret > 0) then
        call jeveuo(noma//'.DOMJOINTS', 'L', jdojoi)
        call jelira(noma//'.DOMJOINTS', 'LONMAX', nbjoin, k8bid)
    !   CREATION DU GRAPH LOCAL
        do iaux = 1, nbjoin
            zi(jgraco + rang*nbproc + zi(jdojoi + iaux - 1)) = 1
        end do
        nbjoin = nbjoin/2
    endif

    n4e = nbproc
!   ON COMMUNIQUE POUR SAVOIR QUI EST EN RELATION AVEC QUI
    do iaux = 0, nbproc - 1
        iaux4 = iaux
        call asmpi_bcast_i(zi(jgraco+iaux*nbproc), n4e, iaux4, mpicou)
    end do

    nbedge = 0
    do iaux = 1, nbproc*nbproc
       if (zi(jgraco+iaux-1) .eq. 1) nbedge = nbedge+1
    end do
    nbedge = nbedge/2

!   RECHERCHE DES COUPLAGES MAXIMAUX
    call wkvect('&&BUCOGR.MASQUE', 'V V I', nbproc*nbproc, jmasqu)
    call wkvect('&&BUCOGR.TMP', 'V V I', nbproc, jtmp)
    nmatch = 1
60  continue

    do iproc1 = 0, nbproc - 1
        do iproc2 = 0, nbproc - 1
            posit = iproc1*nbproc + iproc2
            if (zi(jgraco + posit) .eq. 1 .and. zi(jtmp + iproc1) .eq. 0 &
                .and. zi(jtmp + iproc2) .eq. 0) then
                zi(jgraco + posit) = 0
                zi(jmasqu + posit) = nmatch
                posit = iproc2*nbproc + iproc1
                zi(jgraco + posit) = 0
                zi(jmasqu + posit) = nmatch
                nbedge = nbedge-1
                zi(jtmp + iproc1) = 1
                zi(jtmp + iproc2) = 1
            endif
        end do
    end do

    nmatch = nmatch + 1
    do iaux = 0, nbproc - 1
        zi(jtmp + iaux) = 0
    end do
    if (nbedge .gt. 0) goto 60

    nmatch = nmatch - 1
    call wkvect('&&BUCOGR.ORDJOI', 'V V I', nmatch, jordjo)
    do iaux = 0, nmatch - 1
        zi(jordjo + iaux) = -1
    end do

    do iaux = 0, nbproc - 1
       num = zi(jmasqu + rang*nbproc + iaux)
       ASSERT(num .le. nmatch)
       if (num .ne. 0) then
           zi(jordjo + num - 1) = iaux
       endif
    end do

    call wkvect(nomgr, base//' V I', 1 + nmatch, jnbjoi)
    zi(jnbjoi) = nmatch
    do iaux = 1, nmatch
        zi(jnbjoi + iaux) = zi(jordjo + iaux - 1)
    enddo
!
    call jedetc('V', '&&BUCOGR', 1)
!
    call jedema()
#endif
!
end subroutine
