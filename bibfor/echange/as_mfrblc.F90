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

subroutine as_mfrblc(fid, nent, nvent, ncent, cs, &
                     swm, stm, pname, start, stride, &
                     count, bsize, lbsize, flt, cret)
! person_in_charge: nicolas.sellenet at edf.fr
!
    implicit none
#include "asterf_types.h"
#include "asterf.h"
#include "asterfort/utmess.h"
#include "med/mfrblc.h"
    med_idt :: fid
    aster_int :: nent, nvent, ncent, cs, swm, stm, start
    aster_int :: stride, count, bsize, lbsize, flt, cret
    character(len=*) :: pname
#ifndef ASTER_HAVE_MED
    call utmess('F', 'FERMETUR_2')
#else
!
#if !ASTER_MED_SAME_INT_IDT
    med_idt :: fidm
    med_int :: nent4, nvent4, ncent4, cs4, swm4, stm4, start4
    med_int :: stride4, count4, bsize4, lbsize4, cret4
    fidm = to_med_idt(fid)
    nent4 = to_med_idt(nent)
    nvent4 = to_med_idt(nvent)
    ncent4 = to_med_idt(ncent)
    cs4 = to_med_idt(cs)
    swm4 = to_med_idt(swm)
    stm4 = to_med_idt(stm)
    start4 = to_med_idt(start)
    stride4 = to_med_idt(stride)
    count4 = to_med_idt(count)
    bsize4 = to_med_idt(bsize)
    lbsize4 = to_med_idt(lbsize)
    call mfrblc(fidm, nent4, nvent4, ncent4, cs4, &
                swm4, stm4, pname, start4, stride4, &
                count4, bsize4, lbsize4, flt, cret4)
    cret = cret4
#else
    call mfrblc(fid, nent, nvent, ncent, cs, &
                swm, stm, pname, start, stride, &
                count, bsize, lbsize, flt, cret)
#endif
!
#endif
end subroutine
