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

subroutine as_mfdcsi(fid, cha, ind, numdt, numo, &
                     dt, cret)
! person_in_charge: nicolas.sellenet at edf.fr
!
!
    implicit none
#include "asterf_types.h"
#include "asterf.h"
#include "asterfort/utmess.h"
#include "med/mfdcsi.h"
    med_idt :: fid
    aster_int :: ind, numdt, numo, cret
    character(len=*) :: cha
    real(kind=8) :: dt
#ifndef ASTER_HAVE_MED
    call utmess('F', 'FERMETUR_2')
#else
!
#if !ASTER_MED_SAME_INT_IDT
    med_idt :: fidm
    med_int :: ind4, numdt4, numo4, cret4
    fidm = to_med_idt(fid)
    ind4 = ind
    call mfdcsi(fidm, cha, ind4, numdt4, numo4, &
                dt, cret4)
    numdt = numdt4
    numo = numo4
    cret = cret4
#else
    call mfdcsi(fid, cha, ind, numdt, numo, &
                dt, cret)
#endif
!
#endif
end subroutine
