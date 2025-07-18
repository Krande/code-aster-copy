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
function zerobj(obj)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/jaexin.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
    aster_logical :: zerobj
    character(len=*) :: obj
! person_in_charge: jacques.pellet at edf.fr
! ----------------------------------------------------------------------
!  BUT : DETERMINER SI UN OBJET JEVEUX EST NUL (OU PAS)
!       OBJ     : NOM DE L'OBJET JEVEUX à TESTER
!
!     RESULTAT:
!       ZEROBJ : .TRUE.    SI LES VALEURS DE OBJ SONT TOUTES NULLES
!                .FALSE.   SINON
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
    character(len=1) :: typsca, xous, genr
    character(len=24) :: obj2
    integer(kind=8) :: i, j, n, jval, long, iret, iexi
!
! -DEB------------------------------------------------------------------
!
    call jemarq()
    zerobj = .true.
    obj2 = obj
!
    call jelira(obj2, 'TYPE', cval=typsca)
    ASSERT(typsca .eq. 'R' .or. typsca .eq. 'C')
    call jelira(obj2, 'XOUS', cval=xous)
    call jelira(obj2, 'XOUS', cval=xous)
    call jelira(obj2, 'GENR', cval=genr)
    ASSERT(genr .eq. 'V')
!
!
!     1) CAS DES OBJETS SIMPLES :
!     --------------------------------
    if (xous .eq. 'S') then
        call jeveuo(obj2, 'L', jval)
        call jelira(obj2, 'LONMAX', long)
!
        if (typsca .eq. 'R') then
            do j = 1, long
                if (zr(jval-1+j) .ne. 0.d0) goto 9998
            end do
        else
            do j = 1, long
                if (zc(jval-1+j) .ne. (0.d0, 0.d0)) goto 9998
            end do
        end if
    end if
!
!
!     2) CAS DES COLLECTIONS :
!     --------------------------------
    if (xous .eq. 'X') then
        call jelira(obj2, 'NMAXOC', n)
!
        do i = 1, n
            call jeexin(jexnum(obj2, i), iret)
            if (iret .eq. 0) goto 10
!         -- SI UN OBJET N'A PAS D'ADRESSE DISQUE, C'EST QU'IL EST NUL :
!            (CELA PEUT ARRIVER SI METHODE='GROUP_ELEM')
            call jaexin(jexnum(obj2, i), iexi)
            if (iexi .eq. 0) goto 10
            call jeveuo(jexnum(obj2, i), 'L', jval)
            call jelira(jexnum(obj2, i), 'LONMAX', long)
!
            if (typsca .eq. 'R') then
                do j = 1, long
                    if (zr(jval-1+j) .ne. 0.d0) goto 9998
                end do
            else
                do j = 1, long
                    if (zc(jval-1+j) .ne. (0.d0, 0.d0)) goto 9998
                end do
            end if
10          continue
        end do
    end if
!
!
!
    goto 999
9998 continue
    zerobj = .false.
!
!
999 continue
    call jedema()
end function
