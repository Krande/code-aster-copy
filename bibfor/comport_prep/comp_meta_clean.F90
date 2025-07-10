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
subroutine comp_meta_clean(comporMetaInfo)
!
    use Metallurgy_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
!
    character(len=19), intent(in) :: comporMetaInfo
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of comportment (metallurgy)
!
! Delete informations about internal variables
!
! --------------------------------------------------------------------------------------------------
!
! In  comporMetaInfo    : name of object for information about internal variables and comportement
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call jedetr(comporMetaInfo(1:19)//'.ZONE')
    call jedetr(comporMetaInfo(1:19)//'.RELA')
    call jedetr(comporMetaInfo(1:19)//'.VARI')
    call jedetr(comporMetaInfo(1:19)//'.INFO')
    call jedema()
!
end subroutine
