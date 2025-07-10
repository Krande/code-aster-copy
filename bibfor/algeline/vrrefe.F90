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
subroutine vrrefe(dsName1Z, dsName2Z, ier)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/idensd.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/assert.h"
!
    character(len=*), intent(in) :: dsName1Z, dsName2Z
    integer(kind=8), intent(out) :: ier
!
! --------------------------------------------------------------------------------------------------
!
! Verification que deux concepts ont meme domaine de definition
!     ==> comparaison des ".REFE"
! les concepts doivent etre de type matr_asse_gd, cham_no
! ou cham_elem
!
! Attention : la verification est "minimum". Ce n'est pas parce que
! les objets .REFE sont coherents que l'organisation des champs est identiques.
!
! --------------------------------------------------------------------------------------------------
!
! in  dsName1Z  : ch*19 : nom du 1-er concept
! in  dsName2Z  : ch*19 : nom du 2-nd concept
! out ier       : is   : code retour
!                = 0 pas d'erreur
!                > 0 nombre de descripteurs differents
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: ok
    integer(kind=8) :: ival1, ival2
    character(len=19) :: dsName1, dsName2
    character(len=24) :: refe1, refe2
    aster_logical :: isMatrAsse, isFieldNode, isFieldElem, lgene
    integer(kind=8) :: irefe1, irefe2, iret
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Initializations
    ier = 0
    dsName1 = dsName1Z
    dsName2 = dsName2Z

! - Detect type of object
    refe1 = dsName1//'.REFE'
    call jeexin(refe1, iret)
    isFieldNode = iret .gt. 0
    refe1 = dsName1//'.REFA'
    call jeexin(refe1, iret)
    isMatrAsse = iret .gt. 0
    refe1 = dsName1//'.CELK'
    call jeexin(refe1, iret)
    isFieldElem = iret .gt. 0

! - Access to object
    if (isFieldNode) then
        refe1 = dsName1//'.REFE'
        refe2 = dsName2//'.REFE'
    else if (isMatrAsse) then
        refe1 = dsName1//'.REFA'
        refe2 = dsName2//'.REFA'
    else if (isFieldElem) then
        refe1 = dsName1//'.CELK'
        refe2 = dsName2//'.CELK'
    else
        ASSERT(ASTER_FALSE)
    end if

!   -- recuperation des longueurs des tableaux de reference
    call jelira(refe1, 'LONMAX', ival1)
    call jelira(refe2, 'LONMAX', ival2)
    if (ival1 .ne. ival2) then
        ier = ier+abs(ival1-ival2)
    end if

!   -- recuperation des tableaux d'informations de reference
    call jeveuo(refe1, 'L', irefe1)
    call jeveuo(refe2, 'L', irefe2)

!   -- controle des references
    if (isMatrAsse) then
        lgene = zk24(irefe1-1+10) .eq. 'GENE'
        if (lgene) then
            if (zk24(irefe1-1+2) .ne. zk24(irefe2-1+2)) ier = ier+1
        else
            if (zk24(irefe1-1+1) .ne. zk24(irefe2-1+1)) ier = ier+1
            if (zk24(irefe1-1+2) .ne. zk24(irefe2-1+2)) ier = ier+1
        end if

    else if (isFieldElem) then
        if (zk24(irefe1) .ne. zk24(irefe2)) ier = ier+1
        if (zk24(irefe1+2) .ne. zk24(irefe2+2)) ier = ier+1
        if (zk24(irefe1+3) .ne. zk24(irefe2+3)) ier = ier+1
        if (zk24(irefe1+4) .ne. zk24(irefe2+4)) ier = ier+1
        if (zk24(irefe1+5) .ne. zk24(irefe2+5)) ier = ier+1
        if (zk24(irefe1+6) .ne. zk24(irefe2+6)) ier = ier+1
!       Quelques options metallurgiques produisent
!       des champs que l'on souhaite combiner :
        if (zk24(irefe1+1) .ne. zk24(irefe2+1)) then
            if (zk24(irefe1+1) (1:5) .ne. 'META_') ier = ier+1
            if (zk24(irefe2+1) (1:5) .ne. 'META_') ier = ier+1
        end if
        call jeveuo(dsName1//'.CELD', 'L', irefe1)
        call jeveuo(dsName2//'.CELD', 'L', irefe2)
        if (zi(irefe1) .ne. zi(irefe2)) ier = ier+1

    else if (isFieldNode) then
        if (zk24(irefe1) .ne. zk24(irefe2)) ier = ier+1
        ok = idensd('NUME_EQUA', zk24(irefe1+1), zk24(irefe2+1))
        if (.not. ok) ier = ier+1

    else
        ASSERT(ASTER_FALSE)
    end if
!
    call jedema()
end subroutine
