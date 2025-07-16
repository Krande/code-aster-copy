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
!
subroutine as_mpfope(fid, nom, acces, comm, cret)
!
    implicit none
#include "asterf.h"
#include "asterf_types.h"
#include "asterfort/utmess.h"
#include "med/mpfope.h"
    med_idt, intent(inout) :: fid
    character(len=*), intent(in) :: nom
    aster_int, intent(in) :: acces
    aster_int, intent(in) :: comm
    aster_int, intent(out) :: cret
#ifndef ASTER_HAVE_MED
    call utmess('F', 'FERMETUR_2')
#else
!
#if !ASTER_MED_SAME_INT_IDT
    med_idt :: fidm
    med_int :: acces4, comm4, info4, cret4
#endif
    cret = 0
    if (cret .eq. 0) then
#ifdef ASTER_HAVE_MPI
#if !ASTER_MED_SAME_INT_IDT
        acces4 = to_med_int(acces)
        comm4 = to_med_int(comm)
        info4 = 0
        call mpfope(fidm, nom, acces4, comm4, info4, &
                    cret4)
        fid = to_med_idt(fidm)
        cret = to_aster_int(cret4)
#else
        call mpfope(fid, nom, acces, comm, 0, &
                    cret)
#endif
#endif
    end if
!
#endif
end subroutine
