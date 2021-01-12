! --------------------------------------------------------------------
! Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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

subroutine as_mfafai(fid, maa, ind, fam, num,&
                     gro, cret)
! person_in_charge: nicolas.sellenet at edf.fr
!
!
    implicit none
#include "asterf_types.h"
#include "asterf.h"
#include "asterfort/utmess.h"
#include "med/mfafai.h"
    med_idt :: fid
    aster_int :: num, cret, ind
    character(len=*) :: maa, fam, gro(*)
#ifndef ASTER_HAVE_MED
    call utmess('F', 'FERMETUR_2')
#else
!
#if !ASTER_MED_SAME_INT_IDT
    med_idt :: fidm
    med_int :: num4
    med_int :: cret4, ind4
!      ASSERT(LEN(MAA).EQ.32)
!      ASSERT(LEN(FAM).EQ.32)
!      ASSERT(LEN(GRO(1)).EQ.80)
    fidm = to_med_idt(fid)
    ind4 = ind
    call mfafai(fidm, maa, ind4, fam, num4,&
                gro, cret4)
    num = num4
    cret = cret4
#else
    call mfafai(fid, maa, ind, fam, num,&
                gro, cret)
#endif
!
#endif
end subroutine
