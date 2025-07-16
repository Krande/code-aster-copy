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
subroutine tefrep(option, fieldTypeName, jvForc)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/tecach.h"
#include "asterfort/tecael.h"
#include "asterfort/utmess.h"
!
    character(len=16), intent(in) :: option
    character(len=*), intent(in) :: fieldTypeName
    integer(kind=8), intent(out) :: jvForc
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Check components of force (options CHAR_MECA_F*)
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  fieldTypeName    : type of field
! Out jvForc           : JEVEUX adress of field
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: itab(8), iCmp, iret, iNode, ico
    integer(kind=8) :: iadzi, iazk24, nbValue, jad, nbNode, nbCmp
    character(len=24) :: valk(2)
    character(len=8) :: cellName
!
! --------------------------------------------------------------------------------------------------
!
    call tecach('OON', fieldTypeName, 'L', iret, nval=8, itab=itab)
    ASSERT(iret .eq. 0 .or. iret .eq. 3)
!
    if (iret .eq. 0) then
        jvForc = itab(1)
!
    else if (iret .eq. 3) then
        jvForc = itab(1)
        nbValue = itab(2)
        nbNode = itab(3)
        jad = itab(8)
        nbCmp = nbValue/nbNode
        ASSERT(jvForc .ne. 0)
        ASSERT(nbValue .eq. nbNode*nbCmp)
!
        do iNode = 1, nbNode
            ico = 0
            do iCmp = 1, nbCmp
                if (zl(jad-1+(iNode-1)*nbCmp+iCmp)) then
                    ico = ico+1
                end if
            end do
            if (ico .ne. 0 .and. ico .ne. nbCmp) goto 12

            if (ico .eq. 0) then
                do iCmp = 1, nbCmp
                    zr(jvForc-1+(iNode-1)*nbCmp+iCmp) = 0.d0
                end do
            end if
        end do
        goto 999
12      continue
        call tecael(iadzi, iazk24)
        cellName = zk24(iazk24-1+3) (1:8)
        valk(1) = fieldTypeName
        valk(2) = option
        call utmess('F', 'CALCUL_18', nk=2, valk=valk, si=zi(iadzi))
    end if
!
999 continue
end subroutine
