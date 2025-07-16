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

subroutine asmpi_comm_logical(op, nbval, scl, vl)
!
    implicit none
#include "asterf_types.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/asmpi_comm_vect.h"
#include "asterfort/assert.h"
!
    aster_logical, intent(inout), optional :: scl, vl(*)
    integer(kind=8), intent(in), optional :: nbval
    character(len=*), intent(in) :: op
!
!
!
!  FONCTION REALISEE : SUR-COUCHE MPI
!
!  FAIRE UN ECHANGE SUR UN LOGICAL
!
! POUR UN VECTEUR ON APPLIQUE PAR COMPOSANTE
!
! Arguments d'appels
! in optmpi :
!       /'MPI_LAND' == True if all mpi value are true else false
!       /'MPI_LOR'  == True if at least one mpi value is true else false
!
#ifdef ASTER_HAVE_MPI
!
    integer(kind=4) :: i
    integer(kind=4), pointer :: vi4(:) => null()
!
    if (present(scl)) then
        if (scl) then
            i = 1
        else
            i = 0
        end if
!
        if (op == "MPI_LAND") then
            call asmpi_comm_vect('MPI_MIN', 'S', sci4=i)
        elseif (op == "MPI_LOR") then
            call asmpi_comm_vect('MPI_MAX', 'S', sci4=i)
        else
            ASSERT(ASTER_FALSE)
        end if
!
        if (i == 0) then
            scl = ASTER_FALSE
        else
            scl = ASTER_TRUE
        end if
    else
        ASSERT(present(vl))
        ASSERT(present(nbval))
        AS_ALLOCATE(vi4=vi4, size=nbval)
        do i = 1, nbval
            if (vl(i)) then
                vi4(i) = 1
            else
                vi4(i) = 0
            end if
        end do
!
        if (op == "MPI_LAND") then
            call asmpi_comm_vect('MPI_MIN', 'S', nbval=nbval, vi4=vi4)
        elseif (op == "MPI_LOR") then
            call asmpi_comm_vect('MPI_MAX', 'S', nbval=nbval, vi4=vi4)
        else
            ASSERT(ASTER_FALSE)
        end if
!
        do i = 1, nbval
            if (vi4(i) == 0) then
                vl(i) = ASTER_FALSE
            else
                vl(i) = ASTER_TRUE
            end if
        end do
!
        AS_DEALLOCATE(vi4=vi4)
    end if
#else
    aster_logical :: tmp
    character(len=8) :: kbid

    tmp = present(nbval)
    tmp = present(scl)
    tmp = present(vl)
    kbid = op
#endif
!
end subroutine
