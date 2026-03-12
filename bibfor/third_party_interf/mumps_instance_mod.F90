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
! Module providing MUMPS instance storage with correct struct layout.
!
! This module MUST be compiled WITHOUT integer-size promotion flags
! (e.g. -fdefault-integer-8, /integer-size:64) so that bare INTEGER
! in the MUMPS struct headers resolves to INTEGER(4), matching the
! LP64 MUMPS library ABI.
!
! Callers compiled WITH integer-size promotion can safely USE this
! module: Fortran module import preserves the type layout from the
! module's compilation (SEQUENCE types have fixed layout).
!
module mumps_instance_mod
#include "asterf.h"
#ifdef ASTER_HAVE_MUMPS
    implicit none
!
! --- MUMPS struct type definitions (LP64: INTEGER = 4 bytes here)
#   include "smumps_struc.h"
#   include "cmumps_struc.h"
#   include "dmumps_struc.h"
#   include "zmumps_struc.h"
!
! --- Maximum number of simultaneous MUMPS instances
    integer(kind=8), parameter :: nmxins = 50
!
! --- Instance metadata (replaces COMMON /mumpsh/)
    character(len=1)  :: roucs(nmxins), precs(nmxins)
    character(len=4)  :: etams(nmxins)
    character(len=14) :: nonus(nmxins)
    character(len=19) :: nomats(nmxins), nosols(nmxins)
!
! --- MUMPS struct arrays (replaces COMMON /mumpss/, /mumpsc/, /mumpsd/, /mumpsz/)
    type(smumps_struc), target :: smps(nmxins)
    type(cmumps_struc), target :: cmps(nmxins)
    type(dmumps_struc), target :: dmps(nmxins)
    type(zmumps_struc), target :: zmps(nmxins)
!
! --- Interfaces to MUMPS library subroutines
    interface
        subroutine smumps(smpsk)
            import :: smumps_struc
            type(smumps_struc) :: smpsk
        end subroutine smumps
        subroutine cmumps(cmpsk)
            import :: cmumps_struc
            type(cmumps_struc) :: cmpsk
        end subroutine cmumps
        subroutine dmumps(dmpsk)
            import :: dmumps_struc
            type(dmumps_struc) :: dmpsk
        end subroutine dmumps
        subroutine zmumps(zmpsk)
            import :: zmumps_struc
            type(zmumps_struc) :: zmpsk
        end subroutine zmumps
    end interface
#else
    implicit none
    integer(kind=8), parameter :: nmxins = 50
#endif
end module mumps_instance_mod
