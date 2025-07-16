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
subroutine utcono(mcfac, mocle, iocc, nomail, ndim, &
                  coor, iret)
    implicit none
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
#include "asterfort/utnono.h"
#include "asterfort/char8_to_int.h"
!
    integer(kind=8) :: iocc, ndim, iret
    real(kind=8) :: coor(*)
    character(len=8) :: nomail
    character(len=*) :: mcfac, mocle(3)
!
!
! ----------------------------------------------------------------------
!     ------------------------------------------------------------------
    integer(kind=8) :: n1, n2, n3, numno, i, ier, jcoor
    integer(kind=8) :: vali(2)
    character(len=8) :: k8b, noeud
    character(len=16) :: concep, cmd
    character(len=24) :: coord, nomgrn
    character(len=24) :: valk(3)
!     ------------------------------------------------------------------
    call jemarq()
    iret = 0
!
    call getvr8(mcfac, mocle(1), iocc=iocc, nbval=0, nbret=n1)
    if (n1 .ne. 0) then
        call getvr8(mcfac, mocle(1), iocc=iocc, nbval=ndim, vect=coor, &
                    nbret=n1)
        if (n1 .lt. ndim) then
            call getres(k8b, concep, cmd)
            valk(1) = mcfac
            vali(1) = iocc
            call utmess('F+', 'MODELISA9_23', sk=valk(1), si=vali(1))
            if (ndim .eq. 2) then
                call utmess('F+', 'MODELISA9_24')
            else
                call utmess('F+', 'MODELISA9_25')
            end if
            vali(1) = abs(n1)
            vali(2) = ndim
            valk(1) = mocle(1)
            call utmess('F', 'MODELISA9_26', sk=valk(1), ni=2, vali=vali)
        end if
        iret = 1
        goto 999
    end if
!
    coord = nomail//'.COORDO    .VALE'
    call jeveuo(coord, 'L', jcoor)
!
    call getvtx(mcfac, mocle(2), iocc=iocc, nbval=0, nbret=n2)
    if (n2 .ne. 0) then
        call getvtx(mcfac, mocle(2), iocc=iocc, scal=noeud, nbret=n2)
        numno = char8_to_int(noeud)
        numno = char8_to_int(noeud)
        if (numno .eq. 0) then
            call getres(k8b, concep, cmd)
            valk(1) = mcfac
            valk(2) = mocle(2)
            valk(3) = noeud
            vali(1) = iocc
            call utmess('F', 'MODELISA9_27', nk=3, valk=valk, si=vali(1))
        end if
        do i = 1, ndim
            coor(i) = zr(jcoor+3*(numno-1)+i-1)
        end do
        iret = 1
        goto 999
    end if
!
    call getvtx(mcfac, mocle(3), iocc=iocc, scal=k8b, nbret=n3)
    if (n3 .ne. 0) then
        call getvtx(mcfac, mocle(3), iocc=iocc, scal=nomgrn, nbret=n3)
        call utnono(' ', nomail, 'NOEUD', nomgrn, k8b, &
                    ier)
        if (ier .eq. 10) then
            call getres(k8b, concep, cmd)
            valk(1) = mcfac
            valk(2) = mocle(3)
            valk(3) = nomgrn
            vali(1) = iocc
            call utmess('F', 'MODELISA9_28', nk=3, valk=valk, si=vali(1))
        else if (ier .eq. 1) then
            call getres(k8b, concep, cmd)
            valk(1) = mcfac
            valk(2) = k8b
            vali(1) = iocc
            call utmess('A', 'MODELISA9_29', nk=2, valk=valk, si=vali(1))
        end if
        numno = char8_to_int(k8b)
        do i = 1, ndim
            coor(i) = zr(jcoor+3*(numno-1)+i-1)
        end do
        iret = 1
        goto 999
    end if
!
999 continue
    call jedema()
end subroutine
