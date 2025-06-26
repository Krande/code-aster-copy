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

subroutine jenonu(nomlu, numo)
    implicit none
#include "jeveux_private.h"
#include "asterfort/assert.h"
#include "asterfort/jjallc.h"
#include "asterfort/jjcroc.h"
#include "asterfort/jjlide.h"
#include "asterfort/jjvern.h"
#include "asterfort/jxveuo.h"
#include "asterfort/utmess.h"
    character(len=*), intent(in) :: nomlu
    integer(kind=8), intent(out) :: numo
!     ==================================================================
!-----------------------------------------------------------------------
    integer(kind=8) :: iadmex, iadmi, ibacol, ipgcex, jcara, jctab, jdate
    integer(kind=8) :: jhcod, jiadd, jiadm, jlong, jlono, jltyp, jluti
    integer(kind=8) :: jmarq, n
!-----------------------------------------------------------------------
    parameter(n=5)
    common/jiatje/jltyp(n), jlong(n), jdate(n), jiadd(n), jiadm(n),&
     &                 jlono(n), jhcod(n), jcara(n), jluti(n), jmarq(n)
    integer(kind=8) :: iclas, iclaos, iclaco, idatos, idatco, idatoc
    common/iatcje/iclas, iclaos, iclaco, idatos, idatco, idatoc
    integer(kind=8) :: ipgc, kdesma(2), lgd, lgduti, kposma(2), lgp, lgputi
    common/iadmje/ipgc, kdesma, lgd, lgduti, kposma, lgp, lgputi
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
    character(len=32) :: noml32
    integer(kind=8) :: icre, iret, itab(1)
!     ------------------------------------------------------------------
    numo = 0
    ipgcex = ipgc
    ipgc = -2
!
    if (len(nomlu) .ne. 32) then
        call utmess('F', 'JEVEUX_24', sk=nomlu)
    end if
!
    icre = 0
    noml32 = nomlu
    call jjvern(noml32, icre, iret)
!
    if (iret .eq. 0) then
        call utmess('F', 'JEVEUX_23', sk=noml32)
    end if
!
    if (iret .eq. 1) then
!       ----- OBJET DE TYPE REPERTOIRE
        iadmi = iadm(jiadm(iclaos)+2*idatos-1)
        iadmex = iadmi
        if (iadmex .eq. 0) then
            call jxveuo('L', itab, iret, jctab)
        end if
        call jjcroc('        ', icre)
        if (iadmex .eq. 0) then
            call jjlide('JENONU', noml32, iret)
        end if
!
    else if (iret .eq. 2) then
!       ----- REPERTOIRE DE COLLECTION --
        call jjallc(iclaco, idatco, 'L', ibacol)
        call jjcroc(noml32(25:32), icre)
        call jjlide('JENONU', noml32(1:24), iret)
    else
        ASSERT(.false.)
    end if
!
    numo = idatoc
    ipgc = ipgcex
!
end subroutine
