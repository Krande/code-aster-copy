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

subroutine op0015()
    implicit none
!     OPERATEUR RESOUDRE
!     ------------------------------------------------------------------
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterfort/assert.h"
#include "asterfort/chpver.h"
#include "asterfort/copisd.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/infmaj.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/resoud.h"
#include "asterfort/titre.h"
!
    integer(kind=8) :: ifm, niv, nb, j1, mxiter, ier
    character(len=8) :: xsol, secmbr, matr, vcine, matf, metres, kvari
    character(len=16) :: concep, nomcmd
    character(len=19) :: solve1, solve2
    complex(kind=8) :: cbid
    real(kind=8) :: eps
    integer(kind=8) :: iret
    cbid = dcmplx(0.d0, 0.d0)
!     ------------------------------------------------------------------
    call jemarq()
!
!
    call infmaj()
    call infniv(ifm, niv)
!
    call getres(xsol, concep, nomcmd)
!
    call getvid('  ', 'MATR', scal=matr, nbret=nb)
    ASSERT(nb .eq. 1)
!
    matf = ' '
    call getvid(' ', 'MATR_PREC', scal=matf, nbret=nb)
!
    call getvid('  ', 'CHAM_NO', scal=secmbr, nbret=nb)
    ASSERT(nb .eq. 1)
    call chpver('F', secmbr, 'NOEU', '*', ier)
!
    vcine = ' '
    call getvid('  ', 'CHAM_CINE', scal=vcine, nbret=nb)
    if (nb .eq. 1) call chpver('F', vcine, 'NOEU', '*', ier)
!
!
!   -- CREATION D'1 SOLVEUR TEMPORAIRE : SOLVE2 (SAUF SI MUMPS)
    if (matf .eq. ' ') then
        call dismoi('SOLVEUR', matr, 'MATR_ASSE', repk=solve1)
    else
        call dismoi('SOLVEUR', matf, 'MATR_ASSE', repk=solve1)
    end if
    call jeveuo(solve1//'.SLVK', 'E', j1)
    metres = zk24(j1-1+1)
    if (metres .ne. 'MUMPS' .and. metres .ne. 'PETSC') then
        solve2 = '&&OP0015.SOLVEUR'
        call copisd('SOLVEUR', 'V', solve1, solve2)
    else
!       -- MUMPS COMME PETSC VERIFIENT QUE LE SOLVEUR LORS DE RESOUD
!          EST LE MEME QUE CELUI DE PRERES. ON EST DONC OBLIGE DE LE
!          MODIFIER
        solve2 = solve1
    end if
!
!     -- MODIFICATION DU SOLVEUR DU FAIT DE CERTAINS MOTS CLES :
    call getvr8(' ', 'RESI_RELA', scal=eps, nbret=nb)
    if (nb .eq. 1) then
        call jeveuo(solve2//'.SLVR', 'E', j1)
        zr(j1-1+2) = eps
    end if
    call getvtx(' ', 'POSTTRAITEMENTS', scal=kvari, nbret=nb)
    if (nb .eq. 1) then
        call jeveuo(solve2//'.SLVK', 'E', j1)
        zk24(j1-1+11) = kvari
    end if
    call getvis(' ', 'NMAX_ITER', scal=mxiter, nbret=nb)
    if (nb .eq. 1) then
        call jeveuo(solve2//'.SLVI', 'E', j1)
        zi(j1-1+2) = mxiter
    end if
    call getvtx(' ', 'ALGORITHME', scal=kvari, nbret=nb)
    if ((nb .eq. 1) .and. (metres .eq. 'PETSC')) then
        call jeveuo(solve2//'.SLVK', 'E', j1)
        zk24(j1-1+6) = kvari
    end if
!
!     -- APPEL A LA ROUTINE RESOUD :
    call resoud(matr, matf, solve2, vcine, 0, &
                secmbr, xsol, 'G', [0.d0], [cbid], &
                ' ', .true._1, 0, iret)
!
    if (metres .ne. 'MUMPS' .and. metres .ne. 'PETSC') then
        call detrsd('SOLVEUR', solve2)
    end if
!
!
    call titre()
    call jedema()
end subroutine
