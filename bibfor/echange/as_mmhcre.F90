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

subroutine as_mmhcre(fid, nom, dim, type, desc, &
                     descdt, typrep, nocomp, unit, cret)
! person_in_charge: nicolas.sellenet at edf.fr
!
!
    implicit none
#include "asterf_types.h"
#include "asterf.h"
#include "asterfort/utmess.h"
#include "med/mmhcre.h"
    character(len=*) :: nom
    character(len=*) :: desc, descdt
    character(len=16) :: nocomp(3), unit(3)
    med_idt :: fid
    aster_int :: dim, type, cret, stunde, typrep
#ifndef ASTER_HAVE_MED
    call utmess('F', 'FERMETUR_2')
#else
!
#if !ASTER_MED_SAME_INT_IDT
    med_idt :: fidm
    med_int :: dim4, type4, cret4, stund4, typre4
    stunde = 1
    fidm = to_med_idt(fid)
    dim4 = dim
    type4 = type
    stund4 = stunde
    typre4 = typrep
    call mmhcre(fidm, nom, dim4, dim4, type4, &
                desc, descdt, stund4, typre4, nocomp, &
                unit, cret4)
    cret = cret4
#else
    stunde = 1
    call mmhcre(fid, nom, dim, dim, type, &
                desc, descdt, stunde, typrep, nocomp, &
                unit, cret)
#endif
!
#endif
end subroutine
