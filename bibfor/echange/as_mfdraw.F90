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

subroutine as_mfdraw(fid, cha, filter, val, locname, &
                     typent, typgeo, numdt, dt, numo, &
                     cret)
! person_in_charge: nicolas.sellenet at edf.fr
!
    implicit none
#include "asterf_types.h"
#include "asterf.h"
#include "asterfort/utmess.h"
#include "med/mfdraw.h"
    character(len=*) :: cha, locname
    med_idt :: fid
    aster_int :: typent, typgeo, cret
    aster_int :: numdt, numo, filter
    real(kind=8) :: dt
    real(kind=8) :: val(*)
#ifndef ASTER_HAVE_MED
    call utmess('F', 'FERMETUR_2')
#else
!
#if !ASTER_MED_SAME_INT_IDT
    med_idt :: fidm
    med_int :: typen4, typge4, cret4
    med_int :: numdt4, numo4
    fidm = to_med_idt(fid)
    typen4 = typent
    typge4 = typgeo
    numdt4 = numdt
    numo4 = numo
    call mfdraw(fidm, cha, numdt4, numo4, dt, &
                typen4, typge4, locname, filter, val, &
                cret4)
    cret = cret4
#else
    call mfdraw(fid, cha, numdt, numo, dt, &
                typent, typgeo, locname, filter, val, &
                cret)
#endif
!
#endif
end subroutine
