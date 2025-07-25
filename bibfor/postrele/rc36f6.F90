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
subroutine rc36f6(nbp12, nbp23, nbp13, nbsigr, nbsg1, &
                  nbsg2, nbsg3, sigr, nocc, saltij)
    implicit none
#include "jeveux.h"
#include "asterfort/jelira.h"
#include "asterfort/jeveuo.h"
#include "asterfort/rc36f4.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: nbp12, nbp23, nbp13, nbsigr, nocc(*), nbsg1, nbsg2, nbsg3
    integer(kind=8) :: sigr(*)
    real(kind=8) :: saltij(*)
!
!     CALCUL DU FACTEUR D'USAGE POUR LES SITUATIONS DE PASSAGE
!        SI NOCC(CHEMINS DE PASSAGE) = 0
!        ALORS IL N'EXISTE PLUS DE CHEMIN DE PASSAGE
!
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
    integer(kind=8) :: nbsips, jnpass, i, k, ioc1, nsitup
    character(len=3) :: typass
!     ------------------------------------------------------------------
!
    typass = '1_2'
    if (nbp12 .eq. 0) goto 999
    call jelira('&&RC32SI.PASSAGE_'//typass, 'LONUTI', nbsips)
    call jeveuo('&&RC32SI.PASSAGE_'//typass, 'L', jnpass)
    do i = 1, nbsips
        nsitup = zi(jnpass+i-1)
        do k = 1, nbsigr
            if (sigr(k) .eq. nsitup) then
                ioc1 = k
                goto 14
            end if
        end do
        call utmess('F', 'POSTRCCM_36')
14      continue
        if (nocc(2*(ioc1-1)+1) .ne. 0 .or. nocc(2*(ioc1-1)+2) .ne. 0) goto 999
    end do
    nbp12 = 0
    call rc36f4(typass, nbp12, nbp23, nbp13, nbsigr, &
                nbsg1, nbsg2, nbsg3, saltij)
!
999 continue
    typass = '2_3'
    if (nbp23 .eq. 0) goto 9997
    call jelira('&&RC32SI.PASSAGE_'//typass, 'LONUTI', nbsips)
    call jeveuo('&&RC32SI.PASSAGE_'//typass, 'L', jnpass)
    do i = 1, nbsips
        nsitup = zi(jnpass+i-1)
        do k = 1, nbsigr
            if (sigr(k) .eq. nsitup) then
                ioc1 = k
                goto 24
            end if
        end do
        call utmess('F', 'POSTRCCM_36')
24      continue
        if (nocc(2*(ioc1-1)+1) .ne. 0 .or. nocc(2*(ioc1-1)+2) .ne. 0) goto 9997
    end do
    nbp23 = 0
    call rc36f4(typass, nbp12, nbp23, nbp13, nbsigr, &
                nbsg1, nbsg2, nbsg3, saltij)
!
9997 continue
    typass = '1_3'
    if (nbp13 .eq. 0) goto 9995
    call jelira('&&RC32SI.PASSAGE_'//typass, 'LONUTI', nbsips)
    call jeveuo('&&RC32SI.PASSAGE_'//typass, 'L', jnpass)
    do i = 1, nbsips
        nsitup = zi(jnpass+i-1)
        do k = 1, nbsigr
            if (sigr(k) .eq. nsitup) then
                ioc1 = k
                goto 34
            end if
        end do
        call utmess('F', 'POSTRCCM_36')
34      continue
        if (nocc(2*(ioc1-1)+1) .ne. 0 .or. nocc(2*(ioc1-1)+2) .ne. 0) goto 9995
    end do
    nbp13 = 0
    call rc36f4(typass, nbp12, nbp23, nbp13, nbsigr, &
                nbsg1, nbsg2, nbsg3, saltij)
!
9995 continue
!
end subroutine
