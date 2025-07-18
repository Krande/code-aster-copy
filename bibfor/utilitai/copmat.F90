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

subroutine copmat(mat_in, mat_out)
    implicit none
#include "jeveux.h"
#include "asterf_types.h"
!
!                              Function
!     _______________________________________________________________________
!    | Extract, inside a temporary work vector, all terms of a given matrix  |
!    |_______________________________________________________________________|
!
! person_in_charge: hassan.berro at edf.fr
!
#include "asterfort/assert.h"
#include "asterfort/jelira.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"

    character(len=8), intent(in) :: mat_in
    real(kind=8), intent(out) :: mat_out(*)

    aster_logical :: lsym
    integer(kind=8) :: neq, jsmhc, nbnonz, nvale, jvale
    integer(kind=8) :: nlong, jval2, kterm, jcoll, iligl
    integer(kind=8) :: i, j

    character(len=1)  :: ktyp
    character(len=14) :: nonu
    character(len=19) :: mat19

    character(len=24), pointer :: refa(:) => null()
    integer(kind=8), pointer           :: smde(:) => null()
    integer(kind=8), pointer           :: smdi(:) => null()
!
! Column-major index storage
#define m_out(i,j) mat_out((j-1)*neq+i)
!
    call jemarq()

    mat19 = mat_in

    call jeveuo(mat19//'.REFA', 'L', vk24=refa)
    nonu = refa(2) (1:14)

    call jeveuo(nonu//'.SMOS.SMDE', 'L', vi=smde)
    neq = smde(1)
    mat_out(1:neq*neq) = 0.d0

    call jeveuo(nonu//'.SMOS.SMDI', 'L', vi=smdi)
    call jeveuo(nonu//'.SMOS.SMHC', 'L', jsmhc)

    nbnonz = smdi(neq)

    call jelira(mat19//'.VALM', 'NMAXOC', nvale)
    if (nvale .eq. 1) then
        lsym = .true.
    else if (nvale .eq. 2) then
        lsym = .false.
    else
        ASSERT(.false.)
    end if

    call jeveuo(jexnum(mat19//'.VALM', 1), 'L', jvale)
    call jelira(jexnum(mat19//'.VALM', 1), 'LONMAX', nlong)
    ASSERT(nlong .eq. nbnonz)
    if (.not. lsym) then
        call jeveuo(jexnum(mat19//'.VALM', 2), 'L', jval2)
        call jelira(jexnum(mat19//'.VALM', 2), 'LONMAX', nlong)
        ASSERT(nlong .eq. nbnonz)
    end if

    call jelira(jexnum(mat19//'.VALM', 1), 'TYPE', cval=ktyp)
    ASSERT(ktyp .eq. 'R')

    jcoll = 1
    do kterm = 1, nbnonz

        if (smdi(jcoll) .lt. kterm) jcoll = jcoll+1
        iligl = zi4(jsmhc-1+kterm)

        i = iligl
        j = jcoll

        m_out(i, j) = zr(jvale-1+kterm)
        if (lsym) then
            m_out(j, i) = zr(jvale-1+kterm)
        else
            if (i .ne. j) then
                m_out(j, i) = zr(jval2-1+kterm)
            end if
        end if
    end do

    call jedema()
end subroutine
