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
subroutine nmlere(sddisc, action, infz, iterat, valr)
!
! person_in_charge: mickael.abbas at edf.fr
!
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
    character(len=19) :: sddisc
    character(len=1) :: action
    character(len=*) :: infz
    integer(kind=8) :: iterat
    real(kind=8) :: valr(*)
!
! ----------------------------------------------------------------------
!
! ROUTINE *_NON_LINE (STRUCTURES DE DONNES - SD DISCRETISATION)
!
! LECTURE/ECRITURE DANS SD STOCKAGE DES RESIDUS
!
! ----------------------------------------------------------------------
!
! IN  SDDISC : SD DISCRETISATION
! IN  ACTION : 'L' OU 'E'
! IN  ITERAT : NUMERO ITERATION NEWTON
! IN  INFO   : TYPE D'INFO A STOCKER OU A LIRE
! I/O VALR   : REEL   A ECRIRE OU A LIRE
!
! ----------------------------------------------------------------------
!
    character(len=24) :: infore
    integer(kind=8) :: jifre
    character(len=24) :: info
    integer(kind=8) :: iter
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
! --- INITIALISATIONS
!
    infore = sddisc(1:19)//'.IFRE'
    call jeveuo(infore, 'E', jifre)
    info = infz
!
    ASSERT((action .eq. 'E') .or. (action .eq. 'L'))
!
    if (info .eq. 'VRELA') then
        if (action .eq. 'E') then
            zr(jifre+3*iterat+1-1) = valr(1)
        else
            valr(1) = zr(jifre+3*iterat+1-1)
        end if
    else if (info .eq. 'VMAXI') then
        if (action .eq. 'E') then
            zr(jifre+3*iterat+2-1) = valr(1)
        else
            valr(1) = zr(jifre+3*iterat+2-1)
        end if
    else if (info .eq. 'VCHAR') then
        if (action .eq. 'E') then
            zr(jifre+3*iterat+3-1) = valr(1)
        else
            valr(1) = zr(jifre+3*iterat+3-1)
        end if
!
    else if (info .eq. 'VRELA_TOUS') then
        if (action .eq. 'L') then
            do iter = 0, iterat
                valr(iter+1) = zr(jifre+3*iter+1-1)
            end do
        else
            ASSERT(.false.)
        end if
    else if (info .eq. 'VMAXI_TOUS') then
        if (action .eq. 'L') then
            do iter = 0, iterat
                valr(iter+1) = zr(jifre+3*iter+2-1)
            end do
        else
            ASSERT(.false.)
        end if
    else if (info .eq. 'VCHAR_TOUS') then
        if (action .eq. 'L') then
            do iter = 0, iterat
                valr(iter+1) = zr(jifre+3*iter+3-1)
            end do
        else
            ASSERT(.false.)
        end if
    else
        ASSERT(.false.)
    end if
!
    call jedema()
!
end subroutine
