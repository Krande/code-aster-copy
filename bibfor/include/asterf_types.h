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

#ifndef ASTERF_TYPES_H_
#define ASTERF_TYPES_H_
!
! Definition of types used by aster
!
#include "asterf.h"
!
#define aster_int integer(kind=ASTER_INT_SIZE)
#define to_aster_int(a) int(a, ASTER_INT_SIZE)
!
#define aster_logical logical(kind=ASTER_LOGICAL_SIZE)
#define to_aster_logical(a) logical(a, ASTER_LOGICAL_SIZE)
!
! convert C "bool" stored as a long (1:true, 0,false) as a logical
#define int_to_logical(a) transfer(a, .true.)
!
#define ASTER_TRUE to_aster_logical(.true.)
#define ASTER_FALSE to_aster_logical(.false.)
!
#ifndef ASTER_HAVE_HDF5
#define ASTER_HDF5_HID_SIZE 4
#endif
#define hid_t integer(kind=ASTER_HDF5_HID_SIZE)
#define to_hid_t(a) int(a, ASTER_HDF5_HID_SIZE)
!
#ifndef ASTER_HAVE_MED
#define ASTER_MED_INT_SIZE 4
#define ASTER_MED_IDT_SIZE 4
#endif
#define med_int integer(kind=ASTER_MED_INT_SIZE)
#define to_med_int(a) int(a, ASTER_MED_INT_SIZE)
#define med_idt integer(kind=ASTER_MED_IDT_SIZE)
#define to_med_idt(a) int(a, ASTER_MED_IDT_SIZE)
#define ASTER_MED_SAME_INT_IDT (ASTER_MED_INT_SIZE == ASTER_INT_SIZE \
                                && ASTER_MED_IDT_SIZE == ASTER_INT_SIZE)
!
#ifndef ASTER_HAVE_MPI
#   define ASTER_MPI_INT_SIZE 4
#endif
#define mpi_int integer(kind=ASTER_MPI_INT_SIZE)
#define mpi_bool logical(kind=ASTER_MPI_INT_SIZE)
#define to_mpi_int(a) int(a, ASTER_MPI_INT_SIZE)
!
#ifndef ASTER_HAVE_MUMPS
#   define ASTER_MUMPS_INT_SIZE 4
#endif
#define mumps_int integer(kind=ASTER_MUMPS_INT_SIZE)
#define to_mumps_int(a) int(a, ASTER_MUMPS_INT_SIZE)
!
#define blas_int integer(kind=ASTER_BLAS_INT_SIZE)
#define to_blas_int(a) int(a, ASTER_BLAS_INT_SIZE)
!
#ifdef ASTER_PETSC_64BIT_INDICES
#   define ASTER_PETSC_INT_SIZE 8
#else
#   define ASTER_PETSC_INT_SIZE 4
#endif
#define petsc_int integer(kind=ASTER_PETSC_INT_SIZE)
#define to_petsc_int(a) int(a, ASTER_PETSC_INT_SIZE)
!
#endif
