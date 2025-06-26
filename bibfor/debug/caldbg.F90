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
subroutine caldbg(inout, ncham, lcham, lparam)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/dbgobj.h"
#include "asterfort/exisd.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/uttr24.h"
!
    integer(kind=8) :: ncham
    character(len=*) :: inout
    character(len=19) :: lcham(*)
    character(len=8) :: lparam(*)
!
! --------------------------------------------------------------------------------------------------
!
!     BUT : IMPRIMER SUR UNE LIGNE LA VALEUR
!           D'UNE LISTE DE CHAMPS POUR COMPARER 2 VERSIONS
!
! --------------------------------------------------------------------------------------------------
!
    character(len=19) :: champ
    character(len=24) :: ojb
    character(len=4) :: inou2
    integer(kind=8) :: i, iret, j
    character(len=24), pointer :: fieldParam(:) => null()
    integer(kind=8), pointer :: fieldIndir(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    inou2 = inout

    if (ncham .ne. 0) then
        AS_ALLOCATE(vk24=fieldParam, size=ncham)
        AS_ALLOCATE(vi=fieldIndir, size=ncham)
        do i = 1, ncham
            fieldParam(i) = lparam(i)
        end do
        call uttr24(fieldParam, ncham)
        do i = 1, ncham
            do j = 1, ncham
                if (fieldParam(i) .eq. lparam(j)) then
                    fieldIndir(i) = j
                end if
            end do
        end do
    end if
!
    do i = 1, ncham
        ASSERT(fieldIndir(i) .gt. 0)
        ASSERT(fieldIndir(i) .le. ncham)
        champ = lcham(fieldIndir(i))
        ojb = " "
        call exisd('CARTE', champ, iret)
        if (iret .gt. 0) ojb = champ//'.VALE'
        call exisd('CHAM_GEOM', champ, iret)
        if (iret .gt. 0) ojb = champ//'.VALE'
        call exisd('CHAM_NO', champ, iret)
        if (iret .gt. 0) ojb = champ//'.VALE'
        call exisd('CHAM_ELEM', champ, iret)
        if (iret .gt. 0) ojb = champ//'.CELV'
        call exisd('RESUELEM', champ, iret)
        if (iret .gt. 0) ojb = champ//'.RESL'
        call dbgobj(ojb, 'OUI', 6, '&&CALCUL|'//inou2//'|'//lparam(fieldIndir(i)))
    end do

    AS_DEALLOCATE(vk24=fieldParam)
    AS_DEALLOCATE(vi=fieldIndir)
!
    call jedema()
!
end subroutine
