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

subroutine as_mmhnme(fid, maa, quoi, typent, typgeo, &
                     typcon, n, cret)
! person_in_charge: nicolas.sellenet at edf.fr
!
!
    implicit none
#include "asterf_types.h"
#include "asterf.h"
#include "asterfort/utmess.h"
#include "med/mmhnme.h"
    character(len=*) :: maa
    med_idt :: fid
    aster_int :: typent, typgeo, cret, typcon, n, quoi, mdnont, mdnoit
    aster_int :: chtseq, chttra
#ifndef ASTER_HAVE_MED
    call utmess('F', 'FERMETUR_2')
#else
!
#if !ASTER_MED_SAME_INT_IDT
    med_idt :: fidm
    med_int :: typen4, typge4, cret4, typco4, n4, quoi4
    med_int :: mdnon4, mdnoi4, chtse4, chttr4
    mdnont = -1
    mdnoit = -1
    fidm = to_med_idt(fid)
    typen4 = typent
    typge4 = typgeo
    typco4 = typcon
    quoi4 = quoi
    mdnon4 = mdnont
    mdnoi4 = mdnoit
    call mmhnme(fidm, maa, mdnon4, mdnoi4, typen4, &
                typge4, quoi4, typco4, chtse4, chttr4, &
                n4, cret4)
    n = n4
    cret = cret4
#else
    mdnont = -1
    mdnoit = -1
    call mmhnme(fid, maa, mdnont, mdnoit, typent, &
                typgeo, quoi, typcon, chtseq, chttra, &
                n, cret)
#endif
!
#endif
end subroutine
