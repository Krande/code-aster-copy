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

subroutine as_mmhaaw(fid, maa, fam, n, filter, &
                     typent, typgeo, cret)
! person_in_charge: nicolas.sellenet at edf.fr
!
    implicit none
#include "asterf_types.h"
#include "asterf.h"
#include "asterfort/conv_int.h"
#include "asterfort/utmess.h"
#include "med/mmhaaw.h"
    med_idt :: fid
    aster_int :: filter
    aster_int :: fam(*), n, typent, typgeo, cret, mdnont, mdnoit, dtype
    character(len=*) :: maa
#ifndef ASTER_HAVE_MED
    call utmess('F', 'FERMETUR_2')
#else
!
#if  !ASTER_MED_SAME_INT_IDT
    med_idt :: fidm
    med_int :: typen4, typge4, cret4
    med_int :: mdnon4, mdnoi4, dtyp4
    med_int, allocatable :: fam4(:)
    mdnont = -1
    mdnoit = -1
    dtyp4 = 4
    fidm = to_med_idt(fid)
    allocate (fam4(n))
    call conv_int('ast->med', n, vi_ast=fam, vi_med=fam4)
    typen4 = typent
    typge4 = typgeo
    mdnon4 = mdnont
    mdnoi4 = mdnoit
    call mmhaaw(fidm, maa, dtyp4, mdnon4, mdnoi4, &
                typen4, typge4, filter, fam4, cret4)
    cret = cret4
    deallocate (fam4)
#else
    mdnont = -1
    mdnoit = -1
    dtype = 4
    call mmhaaw(fid, maa, dtype, mdnont, mdnoit, &
                typent, typgeo, filter, fam, cret)
#endif
!
#endif
end subroutine
