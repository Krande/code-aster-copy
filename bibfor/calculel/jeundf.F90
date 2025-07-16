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
subroutine jeundf(obj, undf0_)
! person_in_charge: jacques.pellet at edf.fr
! A_UTIL
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/ismaem.h"
#include "asterc/r8nnem.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
    character(len=*), intent(in) :: obj
    aster_logical, optional, intent(in) :: undf0_
! ----------------------------------------------------------------------
!     BUT : METTRE A "UNDEF" UN OBJET JEVEUX
!             I   :  ISMAEM()
!             L   :  .FALSE.
!             R   :  R8NNEM()
!             C   :  DCMPLX(R8NNEM(),R8NNEM())
!             K*  : 'XXXXXXXXXXXXXX'
!
!     OBJ   IN/JXVAR  K24 : NOM DE L'OBJET
!     undf0 IN        L   : Initialisation à 0 des R, C et I
!
    character(len=24) :: obj2
    real(kind=8) :: r1undf
    character(len=8) :: typsca, k8df
    character(len=16) :: k16df
    character(len=24) :: k24df
    character(len=32) :: k32df
    character(len=80) :: k80df
    character(len=1) :: xous, type
    complex(kind=8) :: c1undf
    integer(kind=8) :: long, i1undf, ltyp, iad, k
    aster_logical :: undf0
! DEB-------------------------------------------------------------------
!
    call jemarq()
!
!
! - Parameter: is undefined a zero ?
!
    undf0 = ASTER_FALSE
    if (present(undf0_)) then
        undf0 = undf0_
    end if
!
    obj2 = obj
!
    if (undf0) then
        i1undf = 0
        r1undf = 0.d0
    else
        i1undf = ismaem()
        r1undf = r8nnem()
    end if
!
    c1undf = dcmplx(r1undf, r1undf)
    k8df = 'XXXXXXXX'
    k16df = k8df//k8df
    k24df = k16df//k8df
    k32df = k24df//k8df
    k80df = k32df//k32df//k16df
!
!
!     -- DETERMINATION DE TYPSCA :
!     -----------------------------
    call jelira(obj2, 'TYPE', cval=type)
    if (type .eq. 'K') then
        call jelira(obj2, 'LTYP', ltyp)
        if (ltyp .eq. 8) then
            typsca = 'K8'
        else if (ltyp .eq. 16) then
            typsca = 'K16'
        else if (ltyp .eq. 24) then
            typsca = 'K24'
        else if (ltyp .eq. 32) then
            typsca = 'K32'
        else if (ltyp .eq. 80) then
            typsca = 'K80'
        else
            ASSERT(.false.)
        end if
    else
        typsca = type
    end if
!
    call jelira(obj2, 'XOUS', cval=xous)
!     TEST CAS NON PROGRAMME
    ASSERT(xous .ne. 'X')
!
    call jelira(obj2, 'LONMAX', long)
    call jeveuo(obj2, 'E', iad)
!
!
    if (typsca .eq. 'I') then
        do k = 1, long
            zi(iad-1+k) = i1undf
        end do
    else if (typsca .eq. 'L') then
        do k = 1, long
            zl(iad-1+k) = .false.
        end do
    else if (typsca .eq. 'R') then
        do k = 1, long
            zr(iad-1+k) = r1undf
        end do
    else if (typsca .eq. 'C') then
        do k = 1, long
            zc(iad-1+k) = c1undf
        end do
    else if (typsca .eq. 'K8') then
        do k = 1, long
            zk8(iad-1+k) = k8df
        end do
    else if (typsca .eq. 'K16') then
        do k = 1, long
            zk16(iad-1+k) = k16df
        end do
    else if (typsca .eq. 'K24') then
        do k = 1, long
            zk24(iad-1+k) = k24df
        end do
    else if (typsca .eq. 'K32') then
        do k = 1, long
            zk32(iad-1+k) = k32df
        end do
    else if (typsca .eq. 'K80') then
        do k = 1, long
            zk80(iad-1+k) = k80df
        end do
    else
        ASSERT(.false.)
    end if
!
!
    call jedema()
end subroutine
