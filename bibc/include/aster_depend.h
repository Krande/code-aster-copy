/* -------------------------------------------------------------------- */
/* Copyright (C) 1991 - 2025 - EDF R&D - www.code-aster.org             */
/* This file is part of code_aster.                                     */
/*                                                                      */
/* code_aster is free software: you can redistribute it and/or modify   */
/* it under the terms of the GNU General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or    */
/* (at your option) any later version.                                  */
/*                                                                      */
/* code_aster is distributed in the hope that it will be useful,        */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of       */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        */
/* GNU General Public License for more details.                         */
/*                                                                      */
/* You should have received a copy of the GNU General Public License    */
/* along with code_aster.  If not, see <http://www.gnu.org/licenses/>.  */
/* -------------------------------------------------------------------- */

#ifndef ASTER_DEPEND_H_
#define ASTER_DEPEND_H_

/*
 * Supported "platforms" are:
 * - ASTER_PLATFORM_LINUX
 * - ASTER_PLATFORM_DARWIN
 * - ASTER_PLATFORM_FREEBSD
 * - ASTER_PLATFORM_MINGW
 *
 * ASTER_PLATFORM_POSIX means LINUX or DARWIN or FREEBSD.
 *
 */

#include "asterc_config.h"

#ifdef ASTER_PLATFORM_LINUX64
#define ASTER_PLATFORM_LINUX
#endif

#ifdef ASTER_PLATFORM_DARWIN64
#define ASTER_PLATFORM_DARWIN
#endif

#if defined ASTER_PLATFORM_MINGW32 || defined ASTER_PLATFORM_MINGW64 || defined __MINGW32__
#ifndef ASTER_PLATFORM_MINGW
#define ASTER_PLATFORM_MINGW
#endif
#endif

#if ( defined ASTER_PLATFORM_FREEBSD64 ) || ( defined __FreeBSD__ )
#define ASTER_PLATFORM_FREEBSD
#endif

/* test required value */
#if ( !defined ASTER_PLATFORM_POSIX ) && ( !defined ASTER_PLATFORM_MINGW )
#error ERROR ASTER_PLATFORM_POSIX or ASTER_PLATFORM_MINGW is required
#endif
#if ( defined ASTER_PLATFORM_POSIX ) && ( defined ASTER_PLATFORM_MINGW )
#error ERROR only one of ASTER_PLATFORM_POSIX or ASTER_PLATFORM_MINGW, not both
#endif

/* MS Windows platforms */
#ifdef ASTER_PLATFORM_MINGW

/* win64 - use LLP64 model */
#ifdef ASTER_HAVE_64_BITS
#define ASTER_STRLEN_AT_END
#define ASTER_INT_SIZE 8
#define ASTER_REAL8_SIZE 8
#define ASTER_HAVE_LONG_LONG 1
#endif

/* stdcall must be defined explicitly because it does not seem required anywhere */
#define ASTER_STRLEN_AT_END

#else
/* Linux & Unix platforms */
#define ASTER_STRLEN_AT_END

/* end platforms type */
#endif

#ifdef ASTER_HAVE_64_BITS
#define INTEGER_NB_CHIFFRES_SIGNIFICATIFS 19
#define REAL_NB_CHIFFRES_SIGNIFICATIFS 16
#else
#define INTEGER_NB_CHIFFRES_SIGNIFICATIFS 9
#define REAL_NB_CHIFFRES_SIGNIFICATIFS 16
#endif

#define STRING_SIZE ASTER_C_STRING_SIZE
typedef ASTER_C_FORTRAN_INT4 ASTERINTEGER4;
typedef ASTER_C_FORTRAN_INT ASTERINTEGER;
typedef ASTER_C_FORTRAN_REAL8 ASTERDOUBLE;
typedef ASTER_C_FORTRAN_REAL4 ASTERFLOAT;
typedef ASTER_C_FORTRAN_LOGICAL ASTERLOGICAL;

/* flags d'optimisation */
/* taille de bloc dans MULT_FRONT */
#ifdef ASTER_HAVE_64_BITS
#define ASTER_MULT_FRONT_BLOCK_SIZE__ 96
#else
#define ASTER_MULT_FRONT_BLOCK_SIZE__ 32
#endif

#ifndef ASTER_MULT_FRONT_BLOCK_SIZE
#define ASTER_MULT_FRONT_BLOCK_SIZE ASTER_MULT_FRONT_BLOCK_SIZE__
#endif

/* Comportement par défaut des FPE dans matfpe pour les blas/lapack */
/* On non GNU/Linux systems, FPE are always enabled */
#if defined ASTER_PLATFORM_LINUX || defined ASTER_PLATFORM_MINGW
#ifndef ASTER_HAVE_SUPPORT_FPE
#define ASTER_HAVE_SUPPORT_FPE
#endif

#else
#undef ASTER_HAVE_SUPPORT_FPE
#endif

#endif
