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
subroutine crcoch_getloads(listLoad, nbLoad, nb_ondp, v_ondp)
!
    use listLoad_module
    use listLoad_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/getvid.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "LoadTypes_type.h"
!
    character(len=19), intent(in) :: listLoad
    integer(kind=8), intent(out) :: nb_ondp, nbLoad
    character(len=8), pointer :: v_ondp(:)
!
! --------------------------------------------------------------------------------------------------
!
! CREA_RESU /CONV_CHAR
!
! Get loads
!
! --------------------------------------------------------------------------------------------------
!
! In  listLoad       : name of datastructure for list of loads
! Out nbLoad         : number of loads
! Out nb_ondp        : number of ONDE_PLANE loads
! PTR v_ondp         : pointer to list of ONDE_PLANE loads
!
! --------------------------------------------------------------------------------------------------
!
    character(len=4), parameter :: phenom = "MECA"
    character(len=1), parameter :: jvBase = "V"
    character(len=16), parameter :: loadKeyword = 'CONV_CHAR'
    integer(kind=8) :: iocc, nocc
    character(len=8), pointer :: loadFromUser(:) => null()
    type(ListLoad_Prep) :: listLoadPrep
    aster_logical, parameter :: kineExcl = ASTER_TRUE, diriExcl = ASTER_TRUE
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Output parameters
    nb_ondp = 0
    nbLoad = 0

! - Get loads for loads datastructure
    iocc = 1
    call getvid(loadKeyword, 'CHARGE', iocc=iocc, nbval=0, nbret=nocc)
    nbLoad = -nocc
    ASSERT(nbLoad .ge. 1)
    AS_ALLOCATE(vk8=loadFromUser, size=nbLoad)
    call getvid(loadKeyword, 'CHARGE', iocc=iocc, nbval=nbLoad, vect=loadFromUser)

! - Create list of loads datastructure
    listLoadPrep%model = " "
    call creaListLoadFromList(phenom, listLoadPrep, &
                              listLoad, jvBase, &
                              nbLoad, loadFromUser, &
                              kineExcl, diriExcl)

! - Get loads from type
    call detectMecaNeumLoad(listLoad, LOAD_NEUM_PWAVE, &
                            nb_ondp, v_ondp)

! - Clean
    AS_DEALLOCATE(vk8=loadFromUser)
!
    call jedema()
end subroutine
