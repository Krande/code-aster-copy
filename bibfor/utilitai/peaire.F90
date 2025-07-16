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

subroutine peaire(resu, mesh, nbocc)
    implicit none
#include "jeveux.h"
#include "asterfort/getvem.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/peair1.h"
#include "asterfort/tbajli.h"
#include "asterfort/tbajpa.h"
#include "asterfort/tbcrsd.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    integer(kind=8) :: nbocc
    character(len=*) :: resu
    character(len=8), intent(in) :: mesh
!     OPERATEUR   POST_ELEM
!     TRAITEMENT DU MOT CLE-FACTEUR "AIRE_INTERNE"
!     ------------------------------------------------------------------
!
    integer(kind=8) :: nbparr, ibid, iret, iocc, ng, ngb, jgb, igb, nbb
    integer(kind=8) :: ifm, niv, iadgma
    parameter(nbparr=3)
    real(kind=8) :: valpar(nbparr), aire, long
    character(len=3) :: typarr(nbparr)
    character(len=8) :: k8b
    character(len=16) :: noparr(nbparr)
    character(len=24) :: grpma
    complex(kind=8) :: c16b
!     ------------------------------------------------------------------
    data noparr/'GROUP_MA', 'AIRE', 'LONGUEUR'/
    data typarr/'K24', 'R', 'R'/
!     ------------------------------------------------------------------
!
    call jemarq()
    ibid = 0
    c16b = (0.d0, 0.d0)
!
! --- RECUPERATION DU NIVEAU D'IMPRESSION
    call infniv(ifm, niv)
!
    grpma = mesh//'.GROUPEMA'
!
!     --- CREATION DE LA TABLE ---
!
    call tbcrsd(resu, 'G')
    call tbajpa(resu, nbparr, noparr, typarr)
!
    do iocc = 1, nbocc
        call getvem(mesh, 'GROUP_MA', 'AIRE_INTERNE', 'GROUP_MA_BORD', iocc, &
                    0, k8b, ngb)
        if (ngb .ne. 0) then
            ngb = -ngb
            call wkvect('&&PEAIRE.GROUP_NO', 'V V K24', ngb, jgb)
            call getvem(mesh, 'GROUP_MA', 'AIRE_INTERNE', 'GROUP_MA_BORD', iocc, &
                        ngb, zk24(jgb), ng)
            do igb = 1, ngb
                call jeexin(jexnom(grpma, zk24(jgb+igb-1)), iret)
                if (iret .eq. 0) then
                    call utmess('A', 'UTILITAI3_46', sk=zk24(jgb+igb-1))
                    cycle
                end if
                call jelira(jexnom(grpma, zk24(jgb+igb-1)), 'LONMAX', nbb)
                if (nbb .eq. 0) then
                    call utmess('A', 'UTILITAI3_47', sk=zk24(jgb+igb-1))
                    cycle
                end if
                call jeveuo(jexnom(grpma, zk24(jgb+igb-1)), 'L', iadgma)
!
!              BORD DU TROU : CALCUL DE L'AIRE
!
                call peair1(mesh, nbb, zi(iadgma), aire, long)
                valpar(1) = aire
                valpar(2) = long
                call tbajli(resu, nbparr, noparr, [ibid], valpar, &
                            [c16b], zk24(jgb+igb-1), 0)
            end do
            call jedetr('&&PEAIRE.GROUP_NO')
        end if
    end do
!
    call jedema()
end subroutine
