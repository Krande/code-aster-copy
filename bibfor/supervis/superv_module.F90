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

!> This module manages the global values stored during the execution.
!
module superv_module
!
! person_in_charge: mathieu.courtois@edf.fr

! warning on dummy argument (W0104) may occur because of ifdef

    use calcul_module, only: calcul_init
    implicit none
    private

#include "asterf_types.h"
#include "threading_interfaces.h"
#include "asterc/gtopti.h"
#include "asterc/asmpi_comm.h"
#include "asterfort/assert.h"
#include "asterfort/check_aster_allocate.h"
#include "asterfort/deleteCachedObjects.h"
#include "asterfort/foint0.h"
#include "asterfort/jermxd.h"
#include "asterfort/utgtme.h"
#include "asterfort/utmess.h"
#include "asterfort/utptme.h"

    logical :: first = .true.
    integer(kind=8) :: initMaxThreads = 0

    public :: superv_before, superv_after
    public :: asthread_getmax, asthread_setnum, asthread_blasset, asthread_getnum

contains

!>  Initialize the values or reinitialize them between before executing an operator
!
!>  @todo Remove treatments from execop.
    subroutine superv_before()
        mpi_int :: world, current
        integer(kind=8) :: maxThreads, iret
        real(kind=8) :: rval(6), vx(3), v0
        character(len=8) :: k8tab(6)

!   Check MPI communicators: must be equal between operators
        call asmpi_comm('GET_WORLD', world)
        call asmpi_comm('GET', current)
        ASSERT(world == current)
!   OpenMP variables
        if (first) then
            first = .false.
            call gtopti('numthreads', maxThreads, iret)
            initMaxThreads = maxThreads
        end if
        call asthread_setnum(initMaxThreads, blas_max=1)
!   Memory allocation
!       Adjust Jeveux parameters
        k8tab(1) = 'LIMIT_JV'
        k8tab(2) = 'MEM_TOTA'
        k8tab(3) = 'VMSIZE'
        k8tab(4) = 'CMAX_JV'
        k8tab(5) = 'RLQ_MEM'
        k8tab(6) = 'COUR_JV'
        call utgtme(6, k8tab, rval, iret)
        if (rval(3) .gt. 0 .and. rval(3)-rval(6) .lt. rval(5)) then
!           the remaining memory decreased: adjust it
            call utptme('RLQ_MEM ', rval(3)-rval(6), iret)
        end if
        if (rval(2)-rval(5) .ge. 0) then
            v0 = rval(1)
            if ((rval(2)-rval(5)) .gt. v0) then
!               reduce memory limit
                call jermxd((rval(2)-rval(5))*1024*1024, iret)
                if (iret .eq. 0) then
                    k8tab(1) = 'RLQ_MEM'
                    k8tab(2) = 'LIMIT_JV'
                    call utgtme(2, k8tab, rval, iret)
                    if (abs(rval(2)-v0) .gt. v0*0.1d0) then
                        vx(1) = rval(1)
                        vx(2) = rval(3)
                        vx(3) = rval(2)-v0
                        call utmess('I', 'JEVEUX1_73', nr=3, valr=vx)
                    end if
                end if
            end if
        end if
!       Reinit calcul mark in case of exception
        call calcul_init()
!       Reinitialize counter for as_[de]allocate
        call check_aster_allocate(stage=0)
    end subroutine superv_before

!>  Initialize the values or reinitialize them between before executing an operator
!
!>  This subroutine is called after the execution of each operator to clean
!>  the memory of temporary objects (volatile), matrix...
!
!>  @param[in] exception tell if an exception/error will be raised
    subroutine superv_after(exception)
        logical, optional :: exception
        integer(kind=8) :: stage
        stage = 1
        if (present(exception)) then
!           Do not add another error message if an error has been raised
            stage = 2
        end if
!   Memory allocation
!       Check for not deallocated vectors
        call check_aster_allocate(stage)
!
!       Reset commons used for function interpolation
        call foint0()
!       Delete cached and temporary Jeveux objects
        call deleteCachedObjects()
    end subroutine superv_after

!>  Return the current maximum number of available threads
!
!>  @return current number of threads
    function asthread_getmax()
        implicit none
        integer(kind=8) :: asthread_getmax
#ifdef ASTER_HAVE_OPENMP
        asthread_getmax = omp_get_max_threads()
#else
        asthread_getmax = 1
#endif
    end function asthread_getmax

!>  Set the maximum number of threads for OpenMP and Blas
!
!>  @param[in] nbThreads new maximum number of threads
    subroutine asthread_setnum(nbThreads, blas_max)
        implicit none
        integer(kind=8), intent(in) :: nbThreads
        integer(kind=8), intent(in), optional :: blas_max
#ifdef ASTER_HAVE_OPENMP
        call omp_set_num_threads(nbThreads)
#endif
        if (present(blas_max)) then
            if (blas_max .eq. 1) then
                call asthread_blasset(initMaxThreads)
            end if
        end if
    end subroutine asthread_setnum

!>  Set the maximum number of threads for Blas functions
!
!>  @param[in] nbThreads new maximum number of threads for Blas
    subroutine asthread_blasset(nbThreads)
        implicit none
        integer(kind=8), intent(in) :: nbThreads
#ifdef ASTER_HAVE_OPENMP
#ifdef ASTER_HAVE_OPENBLAS
!       no effect, conflicts with numpy init. run_aster sets OPENBLAS_NUM_THREADS=1
        call openblas_set_num_threads(nbThreads)
# endif
#ifdef ASTER_HAVE_MKL
        call mkl_set_num_threads(nbThreads)
# endif
#endif
    end subroutine asthread_blasset

!>  Return the current thread id
!
!>  @return the current thread id
    function asthread_getnum()
        implicit none
        integer(kind=8) :: asthread_getnum
#ifdef ASTER_HAVE_OPENMP
        asthread_getnum = omp_get_thread_num()
#else
        asthread_getnum = 0
#endif
    end function asthread_getnum

end module superv_module
