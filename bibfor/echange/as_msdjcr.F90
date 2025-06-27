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

subroutine as_msdjcr(fid, lmname, jname, des, dom, rmname, cret)
! person_in_charge: nicolas.sellenet at edf.fr
!
!
    implicit none
#include "asterf_config.h"
#include "asterf_types.h"
#include "asterfort/utmess.h"
#include "med/msdjcr.h"
    med_idt :: fid
    aster_int :: dom, cret
    character(len=*) :: lmname, jname, des, rmname

#ifndef ASTER_HAVE_MED
    call utmess('F', 'FERMETUR_2')
#else

#ifdef ASTER_DEBUG_MED
    write (6, *) '=== as_msdjcr fid=', fid
    write (6, *) '=== as_msdjcr lmname=', lmname
    write (6, *) '=== as_msdjcr jname=', jname
    write (6, *) '=== as_msdjcr des=', des
    write (6, *) '=== as_msdjcr dom=', dom
    write (6, *) '=== as_msdjcr rmname=', rmname
#endif

#if !ASTER_MED_SAME_INT_IDT
    med_idt :: fid4
    med_int :: dom4, cret4
!
    fid4 = to_med_idt(fid)
    dom4 = to_med_int(dom)
    call msdjcr(fid4, lmname, jname, des, dom4, rmname, cret4)
    cret = to_aster_int(cret4)
#else
    call msdjcr(fid, lmname, jname, des, dom, rmname, cret)
#endif

#ifdef ASTER_DEBUG_MED
    write (6, *) '=== as_msdjcr cret=', cret
#endif

#endif
end subroutine
