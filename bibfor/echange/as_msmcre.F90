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

subroutine as_msmcre(fid, nom, dim, desc, typrep, &
                     nocomp, unit, cret)
! person_in_charge: nicolas.sellenet at edf.fr
!
!
    implicit none
#include "asterf_types.h"
#include "asterf.h"
#include "asterfort/utmess.h"
#include "med/msmcre.h"
    character(len=*) :: nom
    character(len=*) :: desc
    character(len=16) :: nocomp(3), unit(3)
    med_idt :: fid
    aster_int :: dim, cret, typrep
#ifndef ASTER_HAVE_MED
    call utmess('F', 'FERMETUR_2')
#else
!
#if !ASTER_MED_SAME_INT_IDT
    med_idt :: fidm
    med_int :: dim4, cret4, typre4
    fidm = to_med_idt(fid)
    dim4 = to_med_int(dim)
    typre4 = to_med_int(typrep)
    call msmcre(fidm, nom, dim4, dim4, desc, &
                typre4, nocomp, unit, cret4)
    cret = cret4
#else
    call msmcre(fid, nom, dim, dim, desc, &
                typrep, nocomp, unit, cret)
#endif
!
#endif
end subroutine
