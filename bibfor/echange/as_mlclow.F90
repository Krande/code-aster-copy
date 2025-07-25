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

subroutine as_mlclow(fid, typgeo, refcoo, modeco, ngauss, &
                     gscoo, wg, locname, ndim, nomasu, &
                     cret)
! person_in_charge: nicolas.sellenet at edf.fr
!
!
    implicit none
#include "asterf_types.h"
#include "asterf.h"
#include "asterfort/utmess.h"
#include "med/mlclow.h"
    med_idt :: fid
    aster_int :: typgeo, modeco, ngauss, cret, ndim
    real(kind=8) :: refcoo(*), gscoo(*), wg(*)
    character(len=*) :: locname, nomasu
#ifndef ASTER_HAVE_MED
    call utmess('F', 'FERMETUR_2')
#else
!
#if !ASTER_MED_SAME_INT_IDT
    med_idt :: fidm
    med_int :: typge4, mode_4, ngaus4, cret4, ndim4
    fidm = to_med_idt(fid)
    typge4 = typgeo
    mode_4 = modeco
    ngaus4 = ngauss
    ndim4 = ndim
    call mlclow(fidm, locname, typge4, ndim4, refcoo, &
                mode_4, ngaus4, gscoo, wg, '', &
                nomasu, cret4)
    cret = cret4
#else
    call mlclow(fid, locname, typgeo, ndim, refcoo, &
                modeco, ngauss, gscoo, wg, '', &
                nomasu, cret)
#endif
!
#endif
end subroutine
