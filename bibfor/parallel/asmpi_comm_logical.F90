! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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

subroutine asmpi_comm_logical(op, svl)
!
implicit none
#include "asterf_types.h"
#include "asterfort/asmpi_comm_vect.h"
#include "asterfort/assert.h"
!
    aster_logical, intent(inout) :: svl
    character(len=*), intent(in) :: op
!
!
!
!  FONCTION REALISEE : SUR-COUCHE MPI
!
!  FAIRE UN ECHANGE SUR UN LOGICAL
!
! Arguments d'appels
! in optmpi :
!       /'MPI_LAND' == True if all mpi value are true else false
!       /'MPI_LOR'  == True if at least one mpi value is true else false
!
#ifdef ASTER_HAVE_MPI
!
    integer(kind=4) :: i
!
    if(svl) then
        i = 1
    else
        i = 0
    end if
!
    if( op == "MPI_LAND" ) then
        call asmpi_comm_vect('MPI_MIN', 'I', sci4=i)
        if(i == 0) then
            svl = ASTER_FALSE
        else
            svl = ASTER_TRUE
        end if
    elseif( op == "MPI_LOR" ) then
        call asmpi_comm_vect('MPI_MAX', 'I', sci4=i)
        if(i == 0) then
            svl = ASTER_FALSE
        else
            svl = ASTER_TRUE
        end if
    else
        ASSERT(ASTER_FALSE)
    end if
#else
    aster_logical :: tmp
    character(len=8) :: kbid

    tmp = svl
    kbid = op
#endif
!
end subroutine
