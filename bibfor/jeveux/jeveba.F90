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

subroutine jeveba(nomlu, base)
! person_in_charge: j-pierre.lefebvre at edf.fr
! aslint: disable=W0405,C1002
    use iso_c_binding, only: c_loc, c_ptr, c_f_pointer
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "jeveux_private.h"
#include "asterfort/assert.h"
#include "asterfort/jjvern.h"
#include "asterfort/utmess.h"
!
    character(len=*), intent(in) :: nomlu
    character(len=1), intent(out) :: base
!
!
!   ==================================================================
    integer(kind=8) :: lk1zon, jk1zon, liszon, jiszon
    common/izonje/lk1zon, jk1zon, liszon, jiszon
!-----------------------------------------------------------------------
    integer(kind=8) :: ibacol, iblono, inat, inatb, ixdeso, ixiadd, ixlono
    integer(kind=8) :: jcara, jdate, jdocu, jgenr, jhcod, jiadd, jiadm
    integer(kind=8) :: jlong, jlono, jltyp, jluti, jmarq, jorig, jrnom
    integer(kind=8) :: jtype, lonoi, ltypi, n
!-----------------------------------------------------------------------
    parameter(n=5)
    common/jiatje/jltyp(n), jlong(n), jdate(n), jiadd(n), jiadm(n),&
     &                 jlono(n), jhcod(n), jcara(n), jluti(n), jmarq(n)
!
    common/jkatje/jgenr(n), jtype(n), jdocu(n), jorig(n), jrnom(n)
    integer(kind=8) :: nblmax, nbluti, longbl, kitlec, kitecr, kiadm, iitlec, iitecr
    integer(kind=8) :: nitecr, kmarq
    common/ificje/nblmax(n), nbluti(n), longbl(n),&
     &                 kitlec(n), kitecr(n), kiadm(n),&
     &                 iitlec(n), iitecr(n), nitecr(n), kmarq(n)
    integer(kind=8) :: numatr
    common/idatje/numatr
!     ------------------------------------------------------------------
    integer(kind=8) :: iclas, iclaos, iclaco, idatos, idatco, idatoc
    common/iatcje/iclas, iclaos, iclaco, idatos, idatco, idatoc
    integer(kind=8) :: izr, izc, izl, izk8, izk16, izk24, izk32, izk80
    equivalence(izr, zr), (izc, zc), (izl, zl), (izk8, zk8), (izk16, zk16),&
     &               (izk24, zk24), (izk32, zk32), (izk80, zk80)
! ----------------------------------------------------------------------
    character(len=32) :: noml32
    integer(kind=8) :: icre, iret
!
!   ==================================================================

    noml32 = nomlu
    base = ' '
!
    icre = 0
    call jjvern(noml32, icre, iret)
    if (iret .eq. 0) then
        call utmess('F', 'JEVEUX_26', sk=noml32(1:24))
    end if
    if (iclaos .eq. 1) then
        base = 'G'
    else if (iclaos .eq. 2) then
        base = 'V'
    else
        ASSERT(.false.)
    end if
!
end subroutine
