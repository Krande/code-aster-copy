! --------------------------------------------------------------------
! Copyright (C) 1991 - 2024 - EDF R&D - www.code-aster.org
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
#ifndef ASTERF_DEBUG_H_
#define ASTERF_DEBUG_H_

#include "asterf_config.h"
!
! Here are defined some flags to add debugging informations.
!
! If the flag is defined, a function prints informations on stdout starting
! with the MARKER (to make grep easy).
! If the flag is not defined, the function must be empty macro.
!
! to add all traces, define ASTER_DEBUG_ALL
! to debug MPI (add to asterc_debug.h too), #define ASTER_DEBUG_MPI
!

#ifdef ASTER_DEBUG_ALL
#   define ASTER_DEBUG_ALLOCATE
#   define ASTER_DEBUG_MPI
#   define ASTER_DEBUG_LOC
#endif

! all prints should start with the same marker
#define MARKER "DEBUG: "

! trace AS_ALLOCATE / AS_DEALLOCATE
#ifdef ASTER_DEBUG_ALLOCATE
#   define DEBUG_ALLOCATE(a, b, c) print *, MARKER, a, ':', b, c
#else
#   define DEBUG_ALLOCATE(a, b, c) continue
#endif

! trace MPI communications
#ifdef ASTER_DEBUG_MPI
#   define DEBUG_MPI(a, b, c) print *, MARKER, a, ':', b, c
#else
#   define DEBUG_MPI(a, b, c) continue
#endif

! print localization
#ifdef ASTER_DEBUG_LOC
#   define DEBUG_LOC(label, a, b) write(6,"(1X,A,A,'@',A,':',I4)") MARKER, label, a, b
#else
#   define DEBUG_LOC(label, a, b) continue
#endif

! print hho
! #define ASTER_DEBUG_TIMER
! print localization
#ifdef ASTER_DEBUG_TIMER
#   define DEBUG_TIME(a, b) print *, MARKER, a, ':', b, " sec"
#else
#   define DEBUG_TIME(a, b) continue
#endif

#define DEBUG_TIMER(a) call cpu_time(a)


#endif
