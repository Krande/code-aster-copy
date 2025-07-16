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

subroutine as_mmhyaw(fid, maa, conn, csize, typent, &
                     typgeo, typcon, filter, cret)
! person_in_charge: nicolas.sellenet at edf.fr
!     L'ARGUMENT CSIZE N'EST PAS DANS L'API MED
!
!
    implicit none
#include "asterf_types.h"
#include "asterf.h"
#include "asterfort/conv_int.h"
#include "asterfort/utmess.h"
#include "med/mmhyaw.h"
    character(len=*) :: maa
    med_idt :: fid
    aster_int :: conn(*), csize, typent, typgeo, typcon, cret
    aster_int :: mdnont, mdnoit, filter
    real(kind=8) :: mdnodt
#ifndef ASTER_HAVE_MED
    call utmess('F', 'FERMETUR_2')
#else
!
#if !ASTER_MED_SAME_INT_IDT
    med_int, allocatable :: conn4(:)
    med_idt :: fidm
    med_int :: typen4, typge4, typco4, cret4
    med_int :: mdnon4, mdnoi4
    mdnont = -1
    mdnoit = -1
    mdnodt = -1.d0
    fidm = to_med_idt(fid)
    allocate (conn4(csize))
    call conv_int('ast->med', csize, vi_ast=conn, vi_med=conn4)
    typen4 = typent
    typge4 = typgeo
    typco4 = typcon
    mdnon4 = mdnont
    mdnoi4 = mdnoit
    call mmhyaw(fidm, maa, mdnon4, mdnoi4, mdnodt, &
                typen4, typge4, typco4, filter, conn4, &
                cret4)
    cret = cret4
    deallocate (conn4)
#else
    mdnont = -1
    mdnoit = -1
    mdnodt = -1.d0
    call mmhyaw(fid, maa, mdnont, mdnoit, mdnodt, &
                typent, typgeo, typcon, filter, conn, &
                cret)
#endif
!
#endif
end subroutine
