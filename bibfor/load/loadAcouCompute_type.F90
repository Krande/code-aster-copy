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
module loadAcouCompute_type
!
    implicit none
!
#include "asterf_types.h"
#include "LoadTypes_type.h"
! ==================================================================================================
!
! Global variables - General
!
! Definition of parameters for all AFFE_CHAR_Acou loads
!
! ==================================================================================================

! - Keyword to define load in AFFE_CHAR_ACOU
    character(len=24), parameter :: acouLoadKeyword(LOAD_NEUA_NBTYPE) = (/ &
                                    'VITE_FACE               '/)

! - Name of field to define load in AFFE_CHAR_ACOU
    character(len=6), parameter :: acouLoadField(LOAD_NEUA_NBTYPE) = (/ &
                                   '.VFACE'/)
!
end module loadAcouCompute_type
