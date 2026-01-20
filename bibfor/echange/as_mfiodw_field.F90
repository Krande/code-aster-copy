! --------------------------------------------------------------------
! Copyright (C) 1991 - 2026 - EDF - www.code-aster.org
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

subroutine as_mfiodw_field(fid, cha, desc, cret)
!
!
    implicit none
#include "asterf_types.h"
#include "asterf.h"
#include "asterfort/utmess.h"
#include "asterfort/as_mfinvr.h"
#include "med/mfioex.h"
#include "med/mfiodw.h"
    med_idt :: fid
    aster_int :: cret
    character(len=*) :: cha, desc
#ifndef ASTER_HAVE_MED
    call utmess('F', 'FERMETUR_2')
#else
!
#if !ASTER_MED_SAME_INT_IDT
    med_idt :: fidm
    med_int :: cret4, oexist4, class4
    med_int :: maj4, mini4, rel4
    fidm = to_med_idt(fid)
    ! class4 = 1 <=> field type
    class4 = 1_4
    !
    call as_mfinvr(fidm, maj4, mini4, rel4, cret4)
    if ((cret4 .eq. 0) .and. (maj4 .eq. 4 .and. mini4 .ge. 2 .or. maj4 .ge. 5)) then
        ! On verifie que le champ existe bien avant d'appeler mfiodw
        call mfioex(fidm, class4, cha, oexist4, cret4)
        if (oexist4 .eq. 1) then
#ifdef ASTER_HAVE_MED_MFIODW
            call mfiodw(fidm, class4, cha, desc, cret4)
#else
            cret4 = 0
#endif
            cret = cret4
        else
            cret = -1
        end if
    else
        cret = 0
    end if

#else
    aster_int :: oexist, class
    aster_int :: maj, mini, rel
    ! class = 1 <=> field type
    class = 1
    !
    call as_mfinvr(fid, maj, mini, rel, cret)
    if ((cret .eq. 0) .and. (maj .eq. 4 .and. mini .ge. 2 .or. maj .ge. 5)) then
        ! On verifie que le champ existe bien avant d'appeler mfiodw
        call mfioex(fid, class, cha, oexist, cret)
        if (oexist .eq. 1) then
#ifdef ASTER_HAVE_MED_MFIODW
            call mfiodw(fid, class, cha, desc, cret)
#else
            cret = 0
#endif
        else
            cret = -1
        end if
    else
        cret = 0
    end if
#endif
!
#endif
end subroutine
