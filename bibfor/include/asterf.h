! --------------------------------------------------------------------
! Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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
! person_in_charge: mathieu.courtois@edf.fr
!
#ifndef ASTERF_H_
#define ASTERF_H_

#include "asterf_config.h"

! The fortran preprocessor is compiler dependent
#if  ASTER_STRINGIFY_USE_QUOTES == 1
#   define ASTER_TO_STRING(name)  "name"
#elif ASTER_STRINGIFY_USE_OPERATOR == 1
#   define ASTER_TO_STRING(name)  #name
#else
#   define ASTER_TO_STRING(name)  "?"
#endif

! Exceptions identifiers - keep consistency with astercxx.h
#define ASTER_ERROR 1
#define ASTER_CONVERGENCE_ERROR 2
#define ASTER_INTEGRATION_ERROR 3
#define ASTER_SOLVER_ERROR 4
#define ASTER_CONTACT_ERROR 5
#define ASTER_TIMELIMIT_ERROR 6

#endif
