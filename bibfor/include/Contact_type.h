! --------------------------------------------------------------------
! Copyright (C) 1991 - 2026 - EDF - www.code-aster.org
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
! ==================================================================================================
!
! Predefined types for contact (DEFI_CONTACT)
!
! ==================================================================================================
#define CONT_FORM_UNDEF 0
#define CONT_FORM_DISC  1
#define CONT_FORM_CONT  2
#define CONT_FORM_UNIL  4
#define CONT_FORM_LAC   5

! ==================================================================================================
!
! Algorithm
!
! ==================================================================================================
#define ALGO_FIXE       0
#define ALGO_NEWT       1
#define NB_DATA_CYCL    75

! Total number of coordinates for intersection points (old one, to suppress)
#define SIZE_MAX_INTE_SL 16
