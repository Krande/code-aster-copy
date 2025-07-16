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

subroutine entete()
! person_in_charge: mathieu.courtois at edf.fr
!
!
    use superv_module, only: asthread_getmax
    implicit none
! ----------------------------------------------------------------------
!     ECRITURE DE L'ENTETE
! ----------------------------------------------------------------------
#include "asterc/lihdfv.h"
#include "asterc/limedv.h"
#include "asterc/liscov.h"
#include "asterc/matfpe.h"
#include "asterc/prhead.h"
#include "asterf_types.h"
#include "asterf.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/utmess.h"
    mpi_int :: rank, size
    integer(kind=8) :: vali(3)
!
#ifdef ASTER_HAVE_OPENMP
    integer(kind=8) :: maxThreads
#endif
! ----------------------------------------------------------------------
! --- INFORMATIONS GLOBALES
    call utmess('I', 'SUPERVIS2_98')
    call prhead(1)
! --- CONFIGURATION MPI
    call asmpi_info(rank=rank, size=size)
#ifdef ASTER_HAVE_MPI
    vali(1) = to_aster_int(rank)
    vali(2) = to_aster_int(size)
    call utmess('I', 'SUPERVIS2_11', ni=2, vali=vali)
#else
    call utmess('I', 'SUPERVIS2_12')
#endif
! --- CONFIGURATION OPENMP
#ifdef ASTER_HAVE_OPENMP
    maxThreads = asthread_getmax()
    call utmess('I', 'SUPERVIS2_13', si=maxThreads)
#endif
! --- LIBRARIES HDF5 ET MED
#ifdef ASTER_HAVE_HDF5
!   Disable fpe that may be raised by H5Tinit_float in hdf5<1.14.4
!   Must be done before any other hdf5 use.
    call matfpe(-1)
    call lihdfv(vali(1), vali(2), vali(3))
    call utmess('I', 'SUPERVIS2_14', ni=3, vali=vali)
    call matfpe(+1)
#else
    call utmess('I', 'SUPERVIS2_15')
#endif
#ifdef ASTER_HAVE_MED
    call limedv(vali(1), vali(2), vali(3))
    call utmess('I', 'SUPERVIS2_16', ni=3, vali=vali)
#else
    call utmess('I', 'SUPERVIS2_17')
#endif
#ifdef ASTER_HAVE_MFRONT
#   define vers0 ASTER_MFRONT_VERSION
    call utmess('I', 'SUPERVIS2_27', sk=vers0)
#else
    call utmess('I', 'SUPERVIS2_28')
#endif
! --- LIBRARIES SOLVEURS
! for backward compatibility
#ifdef ASTER_HAVE_MUMPS
#   ifndef ASTER_MUMPS_VERSION
#       define ASTER_MUMPS_VERSION MUMPS_VERSION
#   endif
!   to avoid C1510, use vers1
#   define vers1 ASTER_MUMPS_VERSION
    call utmess('I', 'SUPERVIS2_18', sk=vers1)
#else
    call utmess('I', 'SUPERVIS2_19')
#endif
#ifdef ASTER_HAVE_PETSC
!   to avoid C1510, use vers2
#   define vers2 ASTER_PETSC_VERSION
    call utmess('I', 'SUPERVIS2_25', sk=vers2)
#else
    call utmess('I', 'SUPERVIS2_26')
#endif
#ifdef ASTER_HAVE_SCOTCH
    call liscov(vali(1), vali(2), vali(3))
    call utmess('I', 'SUPERVIS2_20', ni=3, vali=vali)
#else
    call utmess('I', 'SUPERVIS2_21')
#endif
!     SAUT DE LIGNE
    call utmess('I', 'VIDE_1')
end subroutine
